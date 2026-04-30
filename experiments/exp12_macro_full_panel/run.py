"""exp12 — Macro-AUC over the full 6-feature panel (catalog + waveform).

Combines exp07 catalog features with exp11 waveform features per region, then
applies the same pre-reg gates (CI95-lower > 0.5 AND |effect| > 0.05 AND
permutation z > 3 AND p < Bonferroni-corrected α) to the now-complete
6-feature panel.

Pre-reg v1 SHA:  a4f1c6fa093ecea8bd00077fb8dbcab008c70eaa
PRA-2 SHA:       05a4b0f4f7d26b076fc5169c5cda493e9f343652

Outputs:
    full_panel_feature_summary.csv     all (region, window_kind, t_start, t_end)
                                        rows merged with both catalog + waveform
                                        features
    full_panel_macro_auc_table.csv     macro AUC, CI, z, p, gates per
                                        (feature, null_kind)
    full_panel_macro_plot.png          6-feature × 2-null bar chart with
                                        cross-region scatter
    summary.json
"""
from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.region_pipeline import (  # noqa: E402
    PipelineParams,
    auc_with_significance,
    per_feature_auc,
)

EXP_DIR = Path(__file__).resolve().parent
CATALOG_CSV = ROOT / "experiments" / "exp07_macro_pra2" / "feature_summary.csv"
WAVEFORM_DIR = ROOT / "experiments" / "exp11_full_waveform_features"

# All six features per pre-reg v1 §4. Catalog half is unchanged; waveform half
# joins via the shared (region, window_kind, t_start, t_end) keys.
CATALOG_FEATURES = ["b", "b_drift_per_day", "benioff_total_log10",
                    "benioff_curv", "n_above_mc"]
WAVEFORM_FEATURES = ["wf_spectral_slope", "wf_waveform_entropy",
                     "wf_hht_imf1_if_median_hz"]
ALL_FEATURES = CATALOG_FEATURES + WAVEFORM_FEATURES


def _ts() -> str:
    return dt.datetime.now().strftime("%H:%M:%S")


def cross_region_bootstrap(per_region_aucs, n_boot=1000, seed=42):
    rng = np.random.default_rng(seed)
    region_names = list(per_region_aucs.keys())
    n_regions = len(region_names)
    if n_regions == 0:
        return {}
    feature_keys = list(next(iter(per_region_aucs.values())).keys())
    out = {}
    for fk in feature_keys:
        region_aucs = np.array([per_region_aucs[r][fk]["auc"] for r in region_names])
        macro = float(np.nanmean(region_aucs))
        boot_macros = np.empty(n_boot)
        for i in range(n_boot):
            idx = rng.integers(0, n_regions, size=n_regions)
            boot_macros[i] = np.nanmean(region_aucs[idx])
        ci_lo, ci_hi = np.nanpercentile(boot_macros, [2.5, 97.5])
        out[fk] = {
            "macro_auc": macro, "ci_lo": float(ci_lo), "ci_hi": float(ci_hi),
            "per_region_aucs": region_aucs.tolist(),
        }
    return out


def pooled_macro_significance(merged_df, regions, feature_col, null_kind,
                              n_perm=1000, seed=42):
    rng = np.random.default_rng(seed)
    obs_aucs = []
    region_data = []
    for r in regions:
        sub = merged_df[merged_df["region"] == r]
        p = sub.loc[sub["window_kind"] == "precursor", feature_col].to_numpy(dtype=float)
        q = sub.loc[sub["window_kind"] == null_kind, feature_col].to_numpy(dtype=float)
        p = p[np.isfinite(p)]; q = q[np.isfinite(q)]
        region_data.append((p, q))
        obs_aucs.append(per_feature_auc(p, q))
    obs_macro = float(np.nanmean(obs_aucs))
    perm_macros = np.empty(n_perm)
    for i in range(n_perm):
        per_region_perm = []
        for p, q in region_data:
            if len(p) < 2 or len(q) < 2:
                per_region_perm.append(np.nan)
                continue
            pooled = np.concatenate([p, q])
            rng.shuffle(pooled)
            per_region_perm.append(per_feature_auc(pooled[:len(p)], pooled[len(p):]))
        perm_macros[i] = float(np.nanmean(per_region_perm))
    perm_mean = float(np.nanmean(perm_macros))
    perm_std = float(np.nanstd(perm_macros))
    z = (obs_macro - perm_mean) / perm_std if perm_std > 0 else float("nan")
    from scipy.stats import norm
    p_val = 2 * (1 - norm.cdf(abs(z)))
    return {
        "macro_auc": obs_macro,
        "perm_mean": perm_mean,
        "perm_std": perm_std,
        "z": float(z),
        "p_two_sided": float(p_val),
        "per_region_obs_aucs": obs_aucs,
    }


def main() -> int:
    print(f"[{_ts()}] [exp12] start — full 6-feature panel macro-AUC", flush=True)

    # === Load catalog features (exp07) ===
    print(f"[{_ts()}] [exp12] loading catalog features from exp07", flush=True)
    cat_df = pd.read_csv(CATALOG_CSV)
    cat_df["t_start"] = pd.to_datetime(cat_df["t_start"], utc=True, format="ISO8601")
    cat_df["t_end"] = pd.to_datetime(cat_df["t_end"], utc=True, format="ISO8601")
    print(f"[{_ts()}] [exp12]   {len(cat_df)} rows from {sorted(cat_df['region'].unique())}",
          flush=True)

    # === Load per-region waveform features (exp11) ===
    print(f"[{_ts()}] [exp12] loading waveform features from exp11", flush=True)
    wf_dfs = []
    for csv in sorted(WAVEFORM_DIR.glob("*_waveform_features.csv")):
        wf = pd.read_csv(csv)
        wf["t_start"] = pd.to_datetime(wf["t_start"], utc=True, format="ISO8601")
        wf["t_end"] = pd.to_datetime(wf["t_end"], utc=True, format="ISO8601")
        wf_dfs.append(wf)
        print(f"[{_ts()}] [exp12]   {csv.name}: {len(wf)} rows", flush=True)
    if not wf_dfs:
        print(f"[{_ts()}] [exp12] no waveform CSVs in {WAVEFORM_DIR}; aborting", flush=True)
        return 2
    wf_df = pd.concat(wf_dfs, ignore_index=True)

    # === Merge on (region, window_kind, t_start, t_end) ===
    keys = ["region", "window_kind", "t_start", "t_end"]
    merged = cat_df.merge(
        wf_df[keys + ["wf_n_successful", "wf_n_snapshots"] + WAVEFORM_FEATURES],
        on=keys, how="left",
    )
    print(f"[{_ts()}] [exp12] merged: {len(merged)} rows; "
          f"{int(merged['wf_n_successful'].notna().sum())} have waveform features",
          flush=True)
    merged.to_csv(EXP_DIR / "full_panel_feature_summary.csv", index=False)

    # === Determine qualifying regions per pre-reg + PRA-2 §3 (N_pre >= 8) ===
    qual_regions = []
    for r in sorted(merged["region"].unique()):
        n_pre = int(((merged["region"] == r) & (merged["window_kind"] == "precursor")).sum())
        if n_pre >= 8:
            qual_regions.append(r)
    print(f"[{_ts()}] [exp12] qualifying regions (N_pre>=8): {qual_regions}", flush=True)

    # === Compute per-region per-feature AUC ===
    params = PipelineParams()
    per_region_aucs: dict[str, dict] = {}
    for r in qual_regions:
        sub = merged[merged["region"] == r]
        per_region_aucs[r] = {}
        for col in ALL_FEATURES:
            for null_kind in ("null_A", "null_B"):
                p = sub.loc[sub["window_kind"] == "precursor", col].to_numpy(dtype=float)
                q = sub.loc[sub["window_kind"] == null_kind, col].to_numpy(dtype=float)
                res = auc_with_significance(p, q, n_boot=params.n_boot,
                                             n_perm=params.n_perm,
                                             seed=params.random_seed)
                per_region_aucs[r][f"{col}__{null_kind}"] = res

    # === Cross-region bootstrap + permutation z + Bonferroni ===
    boot = cross_region_bootstrap(per_region_aucs, n_boot=1000)

    perm_table = {}
    for col in ALL_FEATURES:
        for null_kind in ("null_A", "null_B"):
            perm_table[f"{col}__{null_kind}"] = pooled_macro_significance(
                merged, qual_regions, feature_col=col, null_kind=null_kind,
                n_perm=params.n_perm, seed=42,
            )

    # PRA-2 §3 Bonferroni: 0.05 / (n_features × 2 nulls × 2 test_regions ×
    # n_qualifying). Test regions deferred to honest-null framing here, so we
    # use the same denominator as exp07 = 5 features × 2 nulls × 2 test ×
    # 4 qualifying = 80, then scale by features ratio. Concretely use the
    # full 6 feature panel: 6 * 2 * 2 * len(qual_regions).
    bonferroni_n = len(ALL_FEATURES) * 2 * 2 * max(len(qual_regions), 1)
    bonferroni_alpha = 0.05 / bonferroni_n

    # === Headline ===
    rows = []
    print(f"\n[{_ts()}] [headline] full 6-feature panel macro AUC, "
          f"Bonferroni α = 0.05/{bonferroni_n} = {bonferroni_alpha:.2e}", flush=True)
    print(f"{'feature':<26s} {'null':<7s} {'macro':>6s} {'CI95':>14s} "
          f"{'z':>6s} {'p':>10s}  {'3σ':>3s} {'Bonf':>4s}", flush=True)
    for col in ALL_FEATURES:
        for null_kind in ("null_A", "null_B"):
            fk = f"{col}__{null_kind}"
            b = boot[fk]; pm = perm_table[fk]
            ci_lo_above = b["ci_lo"] > 0.5
            ci_hi_below = b["ci_hi"] < 0.5
            effect = abs(b["macro_auc"] - 0.5)
            sigma3 = (ci_lo_above or ci_hi_below) and effect > 0.05 and abs(pm["z"]) > 3
            bonf = pm["p_two_sided"] < bonferroni_alpha
            print(f"{col:<26s} {null_kind:<7s} {b['macro_auc']:>6.3f} "
                  f"[{b['ci_lo']:.2f},{b['ci_hi']:.2f}] {pm['z']:+6.2f} "
                  f"{pm['p_two_sided']:>10.2e}  "
                  f"{'Y' if sigma3 else 'n':>3s} {'Y' if bonf else 'n':>4s}", flush=True)
            rows.append({
                "feature": col, "null": null_kind,
                "is_waveform": col in WAVEFORM_FEATURES,
                "macro_auc": b["macro_auc"],
                "ci_lo": b["ci_lo"], "ci_hi": b["ci_hi"],
                "perm_z": pm["z"], "perm_p": pm["p_two_sided"],
                "passes_3sigma": bool(sigma3),
                "passes_bonferroni": bool(bonf),
                "per_region_aucs": b["per_region_aucs"],
            })
    macro_df = pd.DataFrame(rows)
    macro_df.to_csv(EXP_DIR / "full_panel_macro_auc_table.csv", index=False)

    # === Plot ===
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5), dpi=120, sharey=True)
    x = np.arange(len(ALL_FEATURES))
    for ax_i, null_kind in enumerate(("null_A", "null_B")):
        ax = axes[ax_i]
        macros = [boot[f"{c}__{null_kind}"]["macro_auc"] for c in ALL_FEATURES]
        ci_los = [boot[f"{c}__{null_kind}"]["ci_lo"] for c in ALL_FEATURES]
        ci_his = [boot[f"{c}__{null_kind}"]["ci_hi"] for c in ALL_FEATURES]
        err_lo = np.array(macros) - np.array(ci_los)
        err_hi = np.array(ci_his) - np.array(macros)
        colors = ["#7d8aa6" if c in CATALOG_FEATURES else "#a06f3a"
                  for c in ALL_FEATURES]
        ax.bar(x, macros, width=0.7, color=colors, alpha=0.5,
               yerr=[err_lo, err_hi], capsize=5, edgecolor="black")
        for r_idx, region in enumerate(qual_regions):
            for ci, c in enumerate(ALL_FEATURES):
                auc = per_region_aucs[region][f"{c}__{null_kind}"]["auc"]
                ax.scatter(ci + 0.18 * (r_idx - len(qual_regions) / 2 + 0.5) /
                           max(1, len(qual_regions)),
                           auc, s=20, alpha=0.85,
                           color=plt.cm.tab10(r_idx % 10),
                           edgecolor="white", linewidth=0.5,
                           label=region if ci == 0 else None)
        ax.axhline(0.5, color="black", ls="--", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([c.replace("_", "\n") for c in ALL_FEATURES],
                           fontsize=7.5)
        ax.set_ylim(0.20, 0.85)
        ax.set_ylabel("AUC")
        ax.set_title(f"vs {null_kind}  (n_qual={len(qual_regions)})")
        ax.grid(alpha=0.3, axis="y")
        # Annotate catalog vs waveform
        ax.axvspan(-0.5, 4.5, alpha=0.04, color="#7d8aa6")
        ax.axvspan(4.5, len(ALL_FEATURES) - 0.5, alpha=0.04, color="#a06f3a")
        ax.text(2, 0.83, "catalog", ha="center", fontsize=8, color="#5d6a86")
        ax.text(6, 0.83, "waveform", ha="center", fontsize=8, color="#7d4f1a")
        if ax_i == 1 and qual_regions:
            ax.legend(fontsize=8, loc="lower right", title="region")
    fig.suptitle("exp12: Full 6-feature panel macro-AUC under PRA-2", fontsize=11)
    fig.tight_layout()
    fig.savefig(EXP_DIR / "full_panel_macro_plot.png")
    plt.close(fig)
    print(f"[{_ts()}] [plot] full_panel_macro_plot.png", flush=True)

    n_3 = sum(1 for r in rows if r["passes_3sigma"])
    n_b = sum(1 for r in rows if r["passes_bonferroni"])
    n_3_wf = sum(1 for r in rows if r["passes_3sigma"] and r["is_waveform"])
    n_b_wf = sum(1 for r in rows if r["passes_bonferroni"] and r["is_waveform"])
    print(f"\n[{_ts()}] [exp12] HEADLINE", flush=True)
    print(f"[{_ts()}] [exp12]   total: {n_3}/{len(rows)} pass 3σ, "
          f"{n_b}/{len(rows)} pass Bonferroni", flush=True)
    print(f"[{_ts()}] [exp12]   waveform half: {n_3_wf} pass 3σ, "
          f"{n_b_wf} pass Bonferroni", flush=True)

    summary = {
        "experiment": "exp12_macro_full_panel",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "pre_reg_v1_sha": "a4f1c6fa093ecea8bd00077fb8dbcab008c70eaa",
        "pra2_sha": "05a4b0f4f7d26b076fc5169c5cda493e9f343652",
        "qualifying_regions": qual_regions,
        "n_features_total": len(ALL_FEATURES),
        "n_catalog_features": len(CATALOG_FEATURES),
        "n_waveform_features": len(WAVEFORM_FEATURES),
        "bonferroni_n_tests": bonferroni_n,
        "bonferroni_alpha": bonferroni_alpha,
        "macro_auc_table": rows,
        "n_pass_3sigma": n_3,
        "n_pass_bonferroni": n_b,
        "n_pass_3sigma_waveform": n_3_wf,
        "n_pass_bonferroni_waveform": n_b_wf,
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[{_ts()}] [persist] summary.json", flush=True)
    print(f"[{_ts()}] [exp12] done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
