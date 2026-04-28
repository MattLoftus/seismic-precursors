"""exp06_cross_regional_macro: macro-AUC pool across 6 training regions.

Applies the locked pre-reg protocol (a4f1c6f) to each of the 6 training
regions (California, Cascadia, Japan, Chile, Turkey, Italy), then aggregates:

    1. Macro-AUC across regions (unweighted mean of per-region AUCs)
    2. Pooled AUC (concatenate all regions' precursor and null windows)
    3. Cross-region bootstrap CI (resample REGIONS with replacement,
       recompute macro-AUC each time)
    4. Per-feature direction-consistency sign test (binomial: how many
       of 6 regions show AUC<0.5 vs >0.5?)
    5. Pre-reg gates: macro-AUC bootstrap CI95-lower > 0.5  AND
       |effect| > 0.05  AND  permutation z > 3, with Bonferroni α = 0.05/36

Test regions (Mexico + Alaska) are HELD OUT and not run here.

Same declared deviations as exp05: catalog M_min=2.5 (preview; pre-reg=1.5)
and Null C skipped (catalog-only).

Outputs:
    catalog_<region>.csv (per region; gitignored)
    region_results.json
    feature_summary.csv  (all regions × all windows)
    macro_auc_table.csv
    macro_auc_plot.png
    sign_test_plot.png
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
    FEATURE_COLS,
    PipelineParams,
    per_feature_auc,
    run_region_pipeline,
)
from src.regions import TRAINING_REGIONS  # noqa: E402

START = dt.datetime(2000, 1, 1)
END = dt.datetime(2025, 1, 1)
EXP_DIR = Path(__file__).resolve().parent

PARAMS = PipelineParams()  # all defaults from pre-reg


def _ts() -> str:
    return dt.datetime.now().strftime("%H:%M:%S")


def cross_region_bootstrap(per_region_aucs: dict[str, dict[str, dict]],
                           n_boot: int = 1000,
                           seed: int = 42) -> dict:
    """Resample REGIONS with replacement; recompute macro-AUC for each draw.

    per_region_aucs: {region_name: {feature_key: {auc, z, ...}}}
    Returns: {feature_key: {macro_auc, ci_lo, ci_hi, n_regions}}
    """
    rng = np.random.default_rng(seed)
    region_names = list(per_region_aucs.keys())
    n_regions = len(region_names)
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
            "macro_auc": macro,
            "ci_lo": float(ci_lo),
            "ci_hi": float(ci_hi),
            "per_region_aucs": region_aucs.tolist(),
            "n_regions_finite": int(np.isfinite(region_aucs).sum()),
        }
    return out


def sign_test(per_region_aucs: np.ndarray) -> dict:
    """Binomial sign test: how many of N regions have AUC ≠ 0.5?

    Reports: n_above (AUC > 0.5), n_below (AUC < 0.5), n_total finite,
    p-values for "all above" or "all below" via two-sided binomial.
    """
    finite = per_region_aucs[np.isfinite(per_region_aucs)]
    n = len(finite)
    if n < 2:
        return {"n": n, "n_above": 0, "n_below": 0, "p_two_sided": float("nan")}
    n_above = int((finite > 0.5).sum())
    n_below = int((finite < 0.5).sum())
    # Two-sided binomial p: P(|X - n/2| >= |observed - n/2|) under p=0.5
    from scipy.stats import binomtest
    extreme = max(n_above, n_below)
    p = binomtest(extreme, n, p=0.5, alternative="two-sided").pvalue
    return {"n": n, "n_above": n_above, "n_below": n_below, "p_two_sided": float(p)}


def pooled_macro_significance(per_region_results, feature_col: str, null_kind: str,
                              n_perm: int = 1000, seed: int = 42) -> dict:
    """Permutation z-score for the macro-AUC null distribution.

    Build the null distribution by independently shuffling labels within each
    region, computing per-region AUC, then averaging. Repeat n_perm times.
    The observed macro-AUC's z-score is (obs - perm_mean) / perm_std.
    """
    rng = np.random.default_rng(seed)
    obs_aucs = []
    region_data = []  # list of (precursor_arr, null_arr) per region
    for r in per_region_results:
        feat_df = r.feat_df
        p = feat_df.loc[feat_df["window_kind"] == "precursor", feature_col].to_numpy(dtype=float)
        q = feat_df.loc[feat_df["window_kind"] == null_kind, feature_col].to_numpy(dtype=float)
        p = p[np.isfinite(p)]
        q = q[np.isfinite(q)]
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
    print(f"[{_ts()}] [exp06] start — 6 training regions", flush=True)
    print(f"[{_ts()}] [exp06] regions: {[r.name for r in TRAINING_REGIONS]}", flush=True)

    region_results = []
    for region in TRAINING_REGIONS:
        catalog_cache = EXP_DIR / f"catalog_{region.name}.csv"
        try:
            res = run_region_pipeline(
                region, START, END, PARAMS,
                catalog_cache=catalog_cache,
                log=lambda m: print(f"[{_ts()}] {m}", flush=True),
            )
        except Exception as e:
            print(f"[{_ts()}] [exp06] {region.name} FAILED: {type(e).__name__}: {e}", flush=True)
            continue
        region_results.append(res)
        print(f"[{_ts()}] [exp06] {region.name} done: "
              f"N_pre={res.n_precursor} N_A={res.n_null_a} N_B={res.n_null_b}", flush=True)

    print(f"\n[{_ts()}] [exp06] aggregating across {len(region_results)} regions", flush=True)

    # Persist all region feature dataframes concatenated
    all_feats = pd.concat([r.feat_df for r in region_results], ignore_index=True)
    all_feats.to_csv(EXP_DIR / "feature_summary.csv", index=False)
    print(f"[{_ts()}] [persist] feature_summary.csv ({len(all_feats)} rows across regions)",
          flush=True)

    # Per-region AUC table
    per_region_aucs = {r.region_name: r.auc_table for r in region_results}

    # Cross-region bootstrap (resample regions)
    print(f"[{_ts()}] [aggregate] cross-region bootstrap...", flush=True)
    boot_table = cross_region_bootstrap(per_region_aucs, n_boot=1000, seed=42)

    # Permutation z for macro-AUC null
    print(f"[{_ts()}] [aggregate] permutation z for macro-AUC...", flush=True)
    perm_table = {}
    for col in FEATURE_COLS:
        for null_kind in ("null_A", "null_B"):
            key = f"{col}__{null_kind}"
            perm_table[key] = pooled_macro_significance(
                region_results, feature_col=col, null_kind=null_kind, n_perm=PARAMS.n_perm, seed=42
            )

    # Sign test
    sign_table = {}
    for fk in boot_table:
        per_region = np.array(boot_table[fk]["per_region_aucs"])
        sign_table[fk] = sign_test(per_region)

    # Headline gates per pre-reg §6.3 + §6.4
    rows = []
    bonferroni_alpha = 0.05 / 36
    print(f"\n[{_ts()}] [headline] pre-reg gates: 3σ = (CI_lo>0.5 AND |eff|>0.05 AND z>3); "
          f"Bonferroni α = {bonferroni_alpha:.2e}", flush=True)
    print(f"[{_ts()}] [headline] {'feature':<22s} {'null':<7s} {'macro':>6s} {'CI95':>14s} "
          f"{'z':>6s} {'p':>10s} {'sign':>10s} {'3σ':>4s} {'Bonf':>5s}", flush=True)
    for col in FEATURE_COLS:
        for null_kind in ("null_A", "null_B"):
            fk = f"{col}__{null_kind}"
            b = boot_table[fk]
            pm = perm_table[fk]
            st = sign_table[fk]
            ci_lo_above = b["ci_lo"] > 0.5
            ci_hi_below = b["ci_hi"] < 0.5
            effect = abs(b["macro_auc"] - 0.5)
            sigma3 = (ci_lo_above or ci_hi_below) and effect > 0.05 and abs(pm["z"]) > 3
            bonf = pm["p_two_sided"] < bonferroni_alpha
            print(f"[{_ts()}] [headline] {col:<22s} {null_kind:<7s} "
                  f"{b['macro_auc']:>6.3f} [{b['ci_lo']:.2f},{b['ci_hi']:.2f}] "
                  f"{pm['z']:+6.2f} {pm['p_two_sided']:>10.2e} "
                  f"{st['n_above']}/{st['n_below']}/{st['n']:>2d} "
                  f"  {'Y' if sigma3 else 'n':>4s} {'Y' if bonf else 'n':>5s}", flush=True)
            rows.append({
                "feature": col, "null": null_kind,
                "macro_auc": b["macro_auc"], "ci_lo": b["ci_lo"], "ci_hi": b["ci_hi"],
                "perm_z": pm["z"], "perm_p": pm["p_two_sided"],
                "n_above_half": st["n_above"], "n_below_half": st["n_below"], "n_regions": st["n"],
                "sign_p_two_sided": st["p_two_sided"],
                "passes_3sigma": bool(sigma3),
                "passes_bonferroni": bool(bonf),
                "per_region_aucs": b["per_region_aucs"],
            })
    macro_df = pd.DataFrame(rows)
    macro_df.to_csv(EXP_DIR / "macro_auc_table.csv", index=False)
    print(f"\n[{_ts()}] [persist] macro_auc_table.csv", flush=True)

    # Plot: macro AUC bars per feature × null, with per-region scatter
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), dpi=120, sharey=True)
    width = 0.7
    x = np.arange(len(FEATURE_COLS))
    for ax_i, null_kind in enumerate(("null_A", "null_B")):
        ax = axes[ax_i]
        macros = [boot_table[f"{c}__{null_kind}"]["macro_auc"] for c in FEATURE_COLS]
        ci_los = [boot_table[f"{c}__{null_kind}"]["ci_lo"] for c in FEATURE_COLS]
        ci_his = [boot_table[f"{c}__{null_kind}"]["ci_hi"] for c in FEATURE_COLS]
        err_lo = np.array(macros) - np.array(ci_los)
        err_hi = np.array(ci_his) - np.array(macros)
        ax.bar(x, macros, width=width, color="#7d8aa6", alpha=0.5,
               yerr=[err_lo, err_hi], capsize=6, edgecolor="black")
        for region_idx, r in enumerate(region_results):
            for ci, c in enumerate(FEATURE_COLS):
                auc = r.auc_table[f"{c}__{null_kind}"]["auc"]
                ax.scatter(ci + 0.18 * (region_idx - len(region_results) / 2 + 0.5) / len(region_results),
                           auc, s=22, alpha=0.85,
                           color=plt.cm.tab10(region_idx % 10),
                           edgecolor="white", linewidth=0.5,
                           label=r.region_name if ci == 0 else None)
        ax.axhline(0.5, color="black", ls="--", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([c.replace("_", "\n") for c in FEATURE_COLS], fontsize=9)
        ax.set_ylim(0.20, 0.80)
        ax.set_ylabel("AUC")
        ax.set_title(f"vs {null_kind}")
        ax.grid(alpha=0.3, axis="y")
        if ax_i == 1:
            ax.legend(fontsize=8, loc="upper right", title="region", ncol=1)
    fig.suptitle("exp06: macro-AUC across 6 training regions (bars: cross-region bootstrap CI95;\n"
                 "dots: per-region AUCs)", fontsize=10)
    fig.tight_layout()
    fig.savefig(EXP_DIR / "macro_auc_plot.png")
    plt.close(fig)
    print(f"[{_ts()}] [plot] macro_auc_plot.png", flush=True)

    # Plot: sign test (per-feature direction consistency)
    fig, ax = plt.subplots(figsize=(9, 5), dpi=120)
    feat_labels = []
    n_above = []
    n_below = []
    for col in FEATURE_COLS:
        for null_kind in ("null_A", "null_B"):
            fk = f"{col}__{null_kind}"
            st = sign_table[fk]
            feat_labels.append(f"{col}\nvs {null_kind[-1]}")
            n_above.append(st["n_above"])
            n_below.append(-st["n_below"])  # negative so "below" goes left
    x = np.arange(len(feat_labels))
    ax.barh(x, n_above, color="#1f5fa6", alpha=0.7, label="AUC > 0.5 (precursor higher)")
    ax.barh(x, n_below, color="#c0392b", alpha=0.7, label="AUC < 0.5 (precursor lower)")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(x)
    ax.set_yticklabels(feat_labels, fontsize=8)
    ax.set_xlabel("# regions (N=6 training)")
    ax.set_title("Sign-consistency across regions: how often does direction align?")
    ax.set_xlim(-6.5, 6.5)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3, axis="x")
    fig.tight_layout()
    fig.savefig(EXP_DIR / "sign_test_plot.png")
    plt.close(fig)
    print(f"[{_ts()}] [plot] sign_test_plot.png", flush=True)

    # Per-region results JSON
    region_results_json = []
    for r in region_results:
        region_results_json.append({
            "region": r.region_name,
            "mc": r.mc, "mc_diag": r.mc_diag,
            "n_total": r.n_total, "n_background": r.n_background,
            "n_targets_all": r.n_targets_all, "n_targets_background": r.n_targets_background,
            "n_precursor": r.n_precursor, "n_precursor_rejected": r.n_precursor_rejected,
            "n_null_a": r.n_null_a, "n_null_b": r.n_null_b,
            "auc_table": {k: {kk: float(vv) if isinstance(vv, (np.floating, float, int)) else vv
                              for kk, vv in v.items() if kk != "perm_mean" and kk != "perm_std"}
                          for k, v in r.auc_table.items()},
        })

    summary = {
        "experiment": "exp06_cross_regional_macro",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "preview_deviations": [
            "Catalog M_min=2.5 (pre-reg=1.5; computational practicality)",
            "Null C skipped (catalog-only, no waveforms)",
        ],
        "params": {"window_days": PARAMS.window_days, "target_M_min": PARAMS.target_m_min,
                   "catalog_M_min": PARAMS.catalog_m_min,
                   "zbz_log_eta_threshold": PARAMS.zbz_log_eta_threshold,
                   "n_null_a": PARAMS.n_null_a, "n_null_b": PARAMS.n_null_b,
                   "n_boot": PARAMS.n_boot, "n_perm": PARAMS.n_perm},
        "n_regions": len(region_results),
        "region_results": region_results_json,
        "macro_auc_table": rows,
        "headline_pass_3sigma": [r for r in rows if r["passes_3sigma"]],
        "headline_pass_bonferroni": [r for r in rows if r["passes_bonferroni"]],
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[{_ts()}] [persist] summary.json", flush=True)

    n_3sigma = sum(1 for r in rows if r["passes_3sigma"])
    n_bonf = sum(1 for r in rows if r["passes_bonferroni"])
    print(f"\n[{_ts()}] [exp06] HEADLINE: {n_3sigma}/10 features pass 3σ, "
          f"{n_bonf}/10 pass Bonferroni-corrected", flush=True)
    print(f"[{_ts()}] [exp06] done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
