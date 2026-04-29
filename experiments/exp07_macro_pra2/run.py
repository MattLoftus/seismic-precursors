"""exp07_macro_pra2: macro-AUC pool under Pre-Registration Amendment v2.

Re-runs the cross-regional macro-AUC test with three amendments to pre-reg v1:

    Amendment 1 — catalog source switches from ANSS ComCat to ISC global bulletin.
    Amendment 2 — overlap-rejection forbidden zone becomes [t', t'+60] (drop pre-event side).
    Amendment 3 — per-region precursor minimum N >= 8 KEPT windows; failing regions drop from macro.

Pre-reg v1 SHA:  a4f1c6fa093ecea8bd00077fb8dbcab008c70eaa
PRA-2 SHA:       05a4b0f4f7d26b076fc5169c5cda493e9f343652
Surfacing exp06: 3a896b53e2eb96fc5acf6f0a3aa0d50b39167f73

Outputs:
    catalog_<region>.csv (per region; gitignored)
    pra2_macro_auc_table.csv
    pra2_macro_auc_plot.png
    pra2_sign_test_plot.png
    pra2_comparison_v1_vs_v2.png
    summary.json
    run.log
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

# === PRA-2 protocol parameters ===
START = dt.datetime(2000, 1, 1)
END = dt.datetime(2025, 1, 1)
EXP_DIR = Path(__file__).resolve().parent

PARAMS = PipelineParams(
    catalog_source="ISC",                      # PRA-2 amendment 1
    null_buffer_days_before=0,                 # PRA-2 amendment 2
    null_buffer_days_after=60,
    min_kept_precursor_per_region=8,           # PRA-2 amendment 3
    catalog_m_min=2.5,                         # preview deviation; pre-reg=1.5
)

# Reference: exp06's per-region results under v1, for the v1-vs-v2 comparison panel
EXP06_NPRE = {
    "California": 25, "Cascadia": 4, "Japan": 0, "Chile": 0, "Turkey": 2, "Italy": 19,
}


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
            "n_regions_finite": int(np.isfinite(region_aucs).sum()),
        }
    return out


def sign_test(per_region_aucs):
    finite = per_region_aucs[np.isfinite(per_region_aucs)]
    n = len(finite)
    if n < 2:
        return {"n": n, "n_above": 0, "n_below": 0, "p_two_sided": float("nan")}
    n_above = int((finite > 0.5).sum())
    n_below = int((finite < 0.5).sum())
    from scipy.stats import binomtest
    extreme = max(n_above, n_below)
    p = binomtest(extreme, n, p=0.5, alternative="two-sided").pvalue
    return {"n": n, "n_above": n_above, "n_below": n_below, "p_two_sided": float(p)}


def pooled_macro_significance(qual_results, feature_col, null_kind, n_perm, seed):
    rng = np.random.default_rng(seed)
    obs_aucs = []
    region_data = []
    for r in qual_results:
        feat_df = r.feat_df
        p = feat_df.loc[feat_df["window_kind"] == "precursor", feature_col].to_numpy(dtype=float)
        q = feat_df.loc[feat_df["window_kind"] == null_kind, feature_col].to_numpy(dtype=float)
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
        "macro_auc": obs_macro, "perm_mean": perm_mean, "perm_std": perm_std,
        "z": float(z), "p_two_sided": float(p_val), "per_region_obs_aucs": obs_aucs,
    }


def main() -> int:
    print(f"[{_ts()}] [exp07] start — PRA-2 amended protocol", flush=True)
    print(f"[{_ts()}] [exp07] catalog_source = {PARAMS.catalog_source}", flush=True)
    print(f"[{_ts()}] [exp07] overlap forbidden zone = "
          f"[t', t'+{PARAMS.null_buffer_days_after}d]  (pre-event side dropped)", flush=True)
    print(f"[{_ts()}] [exp07] min_kept_precursor_per_region = "
          f"{PARAMS.min_kept_precursor_per_region}", flush=True)

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
            print(f"[{_ts()}] [exp07] {region.name} FAILED: {type(e).__name__}: {e}",
                  flush=True)
            continue
        region_results.append(res)
        print(f"[{_ts()}] [exp07] {region.name} done: N_pre={res.n_precursor} "
              f"N_A={res.n_null_a} N_B={res.n_null_b}", flush=True)

    # Apply Amendment 3: drop regions with N_pre < threshold
    qual_results = [r for r in region_results
                    if r.n_precursor >= PARAMS.min_kept_precursor_per_region]
    dropped = [(r.region_name, r.n_precursor) for r in region_results
               if r.n_precursor < PARAMS.min_kept_precursor_per_region]
    print(f"\n[{_ts()}] [exp07] qualifying regions (N_pre>={PARAMS.min_kept_precursor_per_region}): "
          f"{[r.region_name for r in qual_results]}", flush=True)
    if dropped:
        print(f"[{_ts()}] [exp07] dropped: {dropped}", flush=True)

    if len(qual_results) < 3:
        print(f"[{_ts()}] [exp07] WARNING: fewer than 3 qualifying regions; "
              f"macro pool too narrow for cross-regional claim per PRA-2 §3", flush=True)

    # Persist per-region feature data
    if qual_results:
        all_feats = pd.concat([r.feat_df for r in qual_results], ignore_index=True)
        all_feats.to_csv(EXP_DIR / "feature_summary.csv", index=False)
        print(f"[{_ts()}] [persist] feature_summary.csv ({len(all_feats)} rows)", flush=True)

    # Macro stats on qualifying regions only
    per_region_aucs = {r.region_name: r.auc_table for r in qual_results}
    print(f"[{_ts()}] [aggregate] cross-region bootstrap on {len(qual_results)} regions...",
          flush=True)
    boot_table = cross_region_bootstrap(per_region_aucs, n_boot=1000, seed=42)

    perm_table = {}
    for col in FEATURE_COLS:
        for null_kind in ("null_A", "null_B"):
            key = f"{col}__{null_kind}"
            perm_table[key] = pooled_macro_significance(
                qual_results, feature_col=col, null_kind=null_kind,
                n_perm=PARAMS.n_perm, seed=42,
            )

    sign_table = {}
    for fk in boot_table:
        per_region = np.array(boot_table[fk]["per_region_aucs"])
        sign_table[fk] = sign_test(per_region)

    # Apply PRA-2 §3 Bonferroni: alpha = 0.05 / (5 features × 2 nulls × 2 test_regions × n_qual)
    n_qual = len(qual_results)
    bonferroni_n = 5 * 2 * 2 * max(n_qual, 1)
    bonferroni_alpha = 0.05 / bonferroni_n

    rows = []
    print(f"\n[{_ts()}] [headline] PRA-2 gates: 3σ = (CI_lo>0.5 AND |eff|>0.05 AND z>3); "
          f"Bonferroni α = 0.05/{bonferroni_n} = {bonferroni_alpha:.2e}", flush=True)
    print(f"[{_ts()}] [headline] {'feature':<22s} {'null':<7s} "
          f"{'macro':>6s} {'CI95':>14s} {'z':>6s} {'p':>10s} "
          f"{'sign':>9s} {'3σ':>4s} {'Bonf':>5s}", flush=True)
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
                  f"{'Y' if sigma3 else 'n':>4s} {'Y' if bonf else 'n':>5s}", flush=True)
            rows.append({
                "feature": col, "null": null_kind,
                "macro_auc": b["macro_auc"], "ci_lo": b["ci_lo"], "ci_hi": b["ci_hi"],
                "perm_z": pm["z"], "perm_p": pm["p_two_sided"],
                "n_above_half": st["n_above"], "n_below_half": st["n_below"],
                "n_regions": st["n"], "sign_p_two_sided": st["p_two_sided"],
                "passes_3sigma": bool(sigma3), "passes_bonferroni": bool(bonf),
                "per_region_aucs": b["per_region_aucs"],
            })
    macro_df = pd.DataFrame(rows)
    macro_df.to_csv(EXP_DIR / "pra2_macro_auc_table.csv", index=False)

    # Plot: macro AUC bars w/ per-region scatter
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), dpi=120, sharey=True)
    x = np.arange(len(FEATURE_COLS))
    for ax_i, null_kind in enumerate(("null_A", "null_B")):
        ax = axes[ax_i]
        macros = [boot_table[f"{c}__{null_kind}"]["macro_auc"] for c in FEATURE_COLS]
        ci_los = [boot_table[f"{c}__{null_kind}"]["ci_lo"] for c in FEATURE_COLS]
        ci_his = [boot_table[f"{c}__{null_kind}"]["ci_hi"] for c in FEATURE_COLS]
        err_lo = np.array(macros) - np.array(ci_los)
        err_hi = np.array(ci_his) - np.array(macros)
        ax.bar(x, macros, width=0.7, color="#7d8aa6", alpha=0.5,
               yerr=[err_lo, err_hi], capsize=6, edgecolor="black")
        for region_idx, r in enumerate(qual_results):
            for ci, c in enumerate(FEATURE_COLS):
                auc = r.auc_table[f"{c}__{null_kind}"]["auc"]
                ax.scatter(ci + 0.18 * (region_idx - len(qual_results) / 2 + 0.5) /
                           max(1, len(qual_results)),
                           auc, s=25, alpha=0.85,
                           color=plt.cm.tab10(region_idx % 10),
                           edgecolor="white", linewidth=0.5,
                           label=r.region_name if ci == 0 else None)
        ax.axhline(0.5, color="black", ls="--", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([c.replace("_", "\n") for c in FEATURE_COLS], fontsize=9)
        ax.set_ylim(0.20, 0.80)
        ax.set_ylabel("AUC")
        ax.set_title(f"vs {null_kind} (n_qual={len(qual_results)})")
        ax.grid(alpha=0.3, axis="y")
        if ax_i == 1 and qual_results:
            ax.legend(fontsize=8, loc="upper right", title="region")
    fig.suptitle(f"exp07: macro-AUC under PRA-2 (catalog={PARAMS.catalog_source}, "
                 f"overlap [t',+{PARAMS.null_buffer_days_after}d])", fontsize=10)
    fig.tight_layout()
    fig.savefig(EXP_DIR / "pra2_macro_auc_plot.png")
    plt.close(fig)
    print(f"[{_ts()}] [plot] pra2_macro_auc_plot.png", flush=True)

    # v1 vs v2 comparison panel
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=120)
    region_names = [r.name for r in TRAINING_REGIONS]
    v1_n = [EXP06_NPRE.get(n, 0) for n in region_names]
    v2_n = []
    for n in region_names:
        match = next((r for r in region_results if r.region_name == n), None)
        v2_n.append(match.n_precursor if match else 0)
    x = np.arange(len(region_names))
    width = 0.38
    ax.bar(x - width / 2, v1_n, width, label="v1 (exp06)", color="#c0392b", alpha=0.7)
    ax.bar(x + width / 2, v2_n, width, label="v2 (PRA-2)", color="#1f5fa6", alpha=0.7)
    ax.axhline(PARAMS.min_kept_precursor_per_region, color="black", ls="--", lw=0.8,
               label=f"PRA-2 N>={PARAMS.min_kept_precursor_per_region} threshold")
    ax.set_xticks(x); ax.set_xticklabels(region_names, fontsize=9)
    ax.set_ylabel("kept precursor windows")
    ax.set_title("v1 (ComCat + 90-day overlap) vs v2 (ISC + 60-day overlap)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(EXP_DIR / "pra2_comparison_v1_vs_v2.png")
    plt.close(fig)
    print(f"[{_ts()}] [plot] pra2_comparison_v1_vs_v2.png", flush=True)

    summary = {
        "experiment": "exp07_macro_pra2",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "pre_reg_v1_sha": "a4f1c6fa093ecea8bd00077fb8dbcab008c70eaa",
        "pra2_sha": "05a4b0f4f7d26b076fc5169c5cda493e9f343652",
        "params": {
            "catalog_source": PARAMS.catalog_source,
            "catalog_m_min": PARAMS.catalog_m_min,
            "null_buffer_days_before": PARAMS.null_buffer_days_before,
            "null_buffer_days_after": PARAMS.null_buffer_days_after,
            "min_kept_precursor_per_region": PARAMS.min_kept_precursor_per_region,
        },
        "all_regions_attempted": [r.region_name for r in region_results],
        "qualifying_regions": [r.region_name for r in qual_results],
        "dropped_regions_below_min": dropped,
        "v1_v2_precursor_count_comparison": dict(zip(region_names, zip(v1_n, v2_n))),
        "macro_auc_table": rows,
        "headline_pass_3sigma": [r for r in rows if r["passes_3sigma"]],
        "headline_pass_bonferroni": [r for r in rows if r["passes_bonferroni"]],
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[{_ts()}] [persist] summary.json", flush=True)

    n_3 = sum(1 for r in rows if r["passes_3sigma"])
    n_b = sum(1 for r in rows if r["passes_bonferroni"])
    print(f"\n[{_ts()}] [exp07] HEADLINE: {n_3}/10 features pass 3σ, "
          f"{n_b}/10 pass Bonferroni-corrected", flush=True)
    print(f"[{_ts()}] [exp07] qualifying regions: {len(qual_results)}/{len(region_results)} "
          f"({100*len(qual_results)/max(1,len(region_results)):.0f}%)", flush=True)
    print(f"[{_ts()}] [exp07] done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
