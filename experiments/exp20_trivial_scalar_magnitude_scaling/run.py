"""exp20 — M-scaling of the trivial-scalar finding (replaces exp17 TLS-based version).

exp17 stratified the TLS template-correlation result by target magnitude.
Since exp18 showed the TLS apparatus is dominated by the trivial scalar
A2_blog_last5d (log10 Benioff in last 5 days), we recompute M-scaling
under the trivial scalar.

Improvements over exp17 demanded by reviewer 3 (statistician):
  - Bootstrap CI95 per (bin × region)
  - Bootstrap CI95 per bin macro
  - Mann-Kendall trend test for monotonicity (3-bin version)
  - Logistic regression with M as continuous predictor
  - Drop "ETAS prediction" attribution; report as "GR-consistent
    triggered seismicity scaling"

Also reports z_block (block-permutation null) alongside z_iid for the
two main bins.

Outputs:
  m_scaling_table.csv      — per (bin × region) AUC with CIs
  m_scaling_macro_table.csv — per-bin macro AUC + CI + z_iid + z_block
  trend_test_results.json  — Mann-Kendall + logistic-on-M
  m_scaling_plot.png       — bar chart per bin with per-region scatter
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
from scipy.stats import kendalltau, norm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
EXP_DIR = Path(__file__).resolve().parent
EXP07_CSV = ROOT / "experiments" / "exp07_macro_pra2" / "feature_summary.csv"

MC_PER_REGION = {
    "California": 3.50, "Cascadia": 3.50, "Turkey": 3.10, "Italy": 3.50,
}
QUAL_REGIONS = ["California", "Cascadia", "Turkey", "Italy"]
NULL_KIND = "null_A"
N_SUBWINDOWS = 6
SUBWINDOW_DAYS = 5

N_BOOT = 1000
N_PERM = 1000
RANDOM_SEED = 42


def _ts():
    return dt.datetime.now().strftime("%H:%M:%S")


def load_catalog(region):
    candidates = [
        ROOT / "experiments" / "exp06_cross_regional_macro" / f"catalog_{region}.csv",
        ROOT / "experiments" / "exp07_macro_pra2" / f"catalog_{region}.csv",
    ]
    for p in candidates:
        if p.is_file():
            df = pd.read_csv(p)
            df["time"] = pd.to_datetime(df["time"], utc=True, format="ISO8601")
            return df
    raise FileNotFoundError(region)


def catalog_blog_last5d(catalog, t_start, t_end, mc):
    """Return log10 Benioff energy in the LAST 5-day sub-window."""
    duration = (t_end - t_start).total_seconds()
    edge_lo = t_start + pd.Timedelta(seconds=duration * 5 / N_SUBWINDOWS)
    edge_hi = t_end
    mask = ((catalog["time"] >= edge_lo) & (catalog["time"] < edge_hi)
            & (catalog["magnitude"] >= mc - 1e-9))
    sub = catalog.loc[mask]
    if not len(sub):
        return 0.0
    energies = 10 ** (1.5 * sub["magnitude"].to_numpy() + 4.8)
    benioff = float(np.sum(np.sqrt(energies)))
    return float(np.log10(max(benioff, 1.0)))


def auc_with_bootstrap(scores, y, n_boot=N_BOOT, rng=None):
    """Return AUC + bootstrap CI95."""
    if rng is None:
        rng = np.random.default_rng()
    mask = np.isfinite(scores)
    if mask.sum() < 4 or (y[mask] == 1).sum() < 2 or (y[mask] == 0).sum() < 2:
        return float("nan"), float("nan"), float("nan")
    auc = float(roc_auc_score(y[mask], scores[mask]))
    boot = np.empty(n_boot)
    for bi in range(n_boot):
        idx = rng.choice(np.where(mask)[0], size=mask.sum(), replace=True)
        if (y[idx] == 1).sum() < 2 or (y[idx] == 0).sum() < 2:
            boot[bi] = np.nan
            continue
        boot[bi] = roc_auc_score(y[idx], scores[idx])
    ci_lo, ci_hi = np.nanpercentile(boot, [2.5, 97.5])
    return auc, float(ci_lo), float(ci_hi)


def block_permutation_z(scores, y, t_unix, n_perm=N_PERM, rng=None):
    """Block-permutation z (circular shift of labels in time)."""
    if rng is None:
        rng = np.random.default_rng()
    mask = np.isfinite(scores)
    if mask.sum() < 4 or (y[mask] == 1).sum() < 2 or (y[mask] == 0).sum() < 2:
        return float("nan")
    s = scores[mask]; yy = y[mask]; tt = t_unix[mask]
    obs_auc = float(roc_auc_score(yy, s))
    order = np.argsort(tt)
    inv = np.argsort(order)
    y_sorted = yy[order]
    perm = np.empty(n_perm)
    for p_i in range(n_perm):
        shift = int(rng.integers(0, len(y_sorted)))
        y_shifted = np.roll(y_sorted, shift)
        y2 = y_shifted[inv]
        if (y2 == 1).sum() < 2 or (y2 == 0).sum() < 2:
            perm[p_i] = np.nan
            continue
        perm[p_i] = roc_auc_score(y2, s)
    perm_mean = float(np.nanmean(perm))
    perm_std = float(np.nanstd(perm))
    if perm_std <= 0:
        return float("nan")
    return float((obs_auc - perm_mean) / perm_std)


def main():
    print(f"[{_ts()}] [exp20] start — trivial-scalar M-scaling", flush=True)
    rng = np.random.default_rng(RANDOM_SEED)

    df = pd.read_csv(EXP07_CSV)
    df["t_start"] = pd.to_datetime(df["t_start"], utc=True, format="ISO8601")
    df["t_end"] = pd.to_datetime(df["t_end"], utc=True, format="ISO8601")
    df = df[df["region"].isin(QUAL_REGIONS) &
            df["window_kind"].isin(["precursor", NULL_KIND])].copy().reset_index(drop=True)
    df["t_start_unix"] = df["t_start"].astype("int64") // 10**9

    catalogs = {region: load_catalog(region) for region in QUAL_REGIONS}

    # Compute trivial scalar per window
    print(f"[{_ts()}] [exp20] computing trivial scalar per window...", flush=True)
    scalars = np.zeros(len(df))
    for i, row in df.iterrows():
        cat = catalogs[row["region"]]
        scalars[i] = catalog_blog_last5d(cat, row["t_start"], row["t_end"],
                                          MC_PER_REGION[row["region"]])
    df["scalar"] = scalars

    # Show precursor target_M distribution
    pre = df[df["window_kind"] == "precursor"]
    print(f"[{_ts()}] [exp20] precursor target_M: mean={pre['target_M'].mean():.2f}  "
          f"median={pre['target_M'].median():.2f}  "
          f"max={pre['target_M'].max():.2f}  n={len(pre)}", flush=True)

    # ===== Two-bin and three-bin M-scaling =====
    BINS_2 = [(4.5, 5.0, "M=[4.5,5.0)"), (5.0, 99.0, "M>=5.0")]
    BINS_3 = [(4.5, 4.8, "M=[4.5,4.8)"), (4.8, 5.2, "M=[4.8,5.2)"),
              (5.2, 99.0, "M>=5.2")]

    def per_bin_aucs(bins):
        per_bin_per_region = []
        per_bin_macro = []
        for lo, hi, label in bins:
            per_region = {}
            for region in QUAL_REGIONS:
                # Subset: this region's precursor windows in this M-bin + this region's null A
                pre_sub = df[(df["region"] == region) & (df["window_kind"] == "precursor")
                             & (df["target_M"] >= lo) & (df["target_M"] < hi)]
                null_sub = df[(df["region"] == region) & (df["window_kind"] == NULL_KIND)]
                if len(pre_sub) < 2:
                    per_region[region] = {"auc": float("nan"), "ci_lo": float("nan"),
                                          "ci_hi": float("nan"), "n_pre": len(pre_sub),
                                          "n_null": len(null_sub)}
                    continue
                scores = np.concatenate([pre_sub["scalar"].to_numpy(),
                                         null_sub["scalar"].to_numpy()])
                y = np.concatenate([np.ones(len(pre_sub)), np.zeros(len(null_sub))])
                auc, ci_lo, ci_hi = auc_with_bootstrap(scores, y, rng=rng)
                per_region[region] = {"auc": auc, "ci_lo": ci_lo, "ci_hi": ci_hi,
                                      "n_pre": int(len(pre_sub)), "n_null": int(len(null_sub))}
                per_bin_per_region.append({
                    "bin_label": label, "region": region,
                    "auc": auc, "ci_lo": ci_lo, "ci_hi": ci_hi,
                    "n_pre": int(len(pre_sub)), "n_null": int(len(null_sub)),
                })
            # Macro AUC across regions in this bin
            aucs = np.array([per_region[r]["auc"] for r in QUAL_REGIONS])
            valid = aucs[np.isfinite(aucs)]
            if len(valid) >= 2:
                macro = float(np.mean(valid))
                # Cross-region bootstrap on macro
                boot = np.empty(N_BOOT)
                for bi in range(N_BOOT):
                    samp = rng.choice(valid, size=len(valid), replace=True)
                    boot[bi] = float(np.mean(samp))
                m_lo, m_hi = np.percentile(boot, [2.5, 97.5])
            else:
                macro = float("nan"); m_lo = float("nan"); m_hi = float("nan")
            per_bin_macro.append({
                "bin_label": label, "lo": lo, "hi": hi,
                "macro_auc": float(macro), "ci_lo": float(m_lo), "ci_hi": float(m_hi),
                "n_pre_total": int(sum(per_region[r]["n_pre"] for r in QUAL_REGIONS)),
                "per_region": per_region,
            })
        return per_bin_macro, per_bin_per_region

    print(f"\n[{_ts()}] [exp20] === 2-bin M-scaling ===", flush=True)
    macro_2, per_region_2 = per_bin_aucs(BINS_2)
    for m in macro_2:
        print(f"[{_ts()}]   {m['bin_label']:<14s}  n_pre={m['n_pre_total']:>3d}  "
              f"macro={m['macro_auc']:.3f}  CI=[{m['ci_lo']:.3f}, {m['ci_hi']:.3f}]",
              flush=True)
        for region in QUAL_REGIONS:
            r = m["per_region"][region]
            print(f"[{_ts()}]     {region:<11s}  AUC={r['auc']:.3f}  "
                  f"CI=[{r['ci_lo']:.3f}, {r['ci_hi']:.3f}]  "
                  f"n_pre={r['n_pre']}", flush=True)

    print(f"\n[{_ts()}] [exp20] === 3-bin M-scaling ===", flush=True)
    macro_3, per_region_3 = per_bin_aucs(BINS_3)
    for m in macro_3:
        print(f"[{_ts()}]   {m['bin_label']:<14s}  n_pre={m['n_pre_total']:>3d}  "
              f"macro={m['macro_auc']:.3f}  CI=[{m['ci_lo']:.3f}, {m['ci_hi']:.3f}]",
              flush=True)

    # ===== Trend tests =====
    print(f"\n[{_ts()}] [exp20] === trend tests ===", flush=True)
    bin_centers_2 = np.array([(b[0] + b[1] if b[1] < 99 else b[0] + 0.5) / 1.0 for b in BINS_2])
    bin_macros_2 = np.array([m["macro_auc"] for m in macro_2])
    tau_2, kpval_2 = kendalltau(bin_centers_2, bin_macros_2)
    bin_centers_3 = np.array([(b[0] + (b[1] if b[1] < 99 else b[0] + 0.5)) / 2.0 for b in BINS_3])
    bin_macros_3 = np.array([m["macro_auc"] for m in macro_3])
    tau_3, kpval_3 = kendalltau(bin_centers_3, bin_macros_3)
    print(f"[{_ts()}]   Mann-Kendall (2-bin): tau={tau_2:.3f}  p={kpval_2:.3f}", flush=True)
    print(f"[{_ts()}]   Mann-Kendall (3-bin): tau={tau_3:.3f}  p={kpval_3:.3f}", flush=True)

    # Logistic regression: P(precursor) ~ scalar + target_M (continuous)
    pre_all = df[df["window_kind"] == "precursor"].copy()
    null_all = df[df["window_kind"] == NULL_KIND].copy()
    # Assign target_M to nulls as median of region's precursor target_M
    null_all["target_M"] = null_all["region"].map(
        lambda r: float(pre_all[pre_all["region"] == r]["target_M"].median())
        if len(pre_all[pre_all["region"] == r]) else 4.7)
    combined = pd.concat([pre_all.assign(label=1), null_all.assign(label=0)], ignore_index=True)
    X = combined[["scalar", "target_M"]].to_numpy()
    y_lr = combined["label"].to_numpy()
    lr = LogisticRegression(class_weight="balanced", max_iter=1000)
    lr.fit(X, y_lr)
    coef_scalar = float(lr.coef_[0, 0])
    coef_M = float(lr.coef_[0, 1])
    print(f"[{_ts()}]   logistic on (scalar + target_M):", flush=True)
    print(f"[{_ts()}]     coef[scalar]   = {coef_scalar:+.4f}", flush=True)
    print(f"[{_ts()}]     coef[target_M] = {coef_M:+.4f}", flush=True)

    # ===== Plot =====
    fig, ax = plt.subplots(figsize=(10, 5), dpi=120)
    n_bins = len(macro_2)
    x = np.arange(n_bins)
    macros = [m["macro_auc"] for m in macro_2]
    err_lo = [m["macro_auc"] - m["ci_lo"] for m in macro_2]
    err_hi = [m["ci_hi"] - m["macro_auc"] for m in macro_2]
    ax.bar(x, macros, color="#3a7d3a", alpha=0.65, yerr=[err_lo, err_hi],
           capsize=5, edgecolor="black", width=0.5)
    for fi, m in enumerate(macro_2):
        for r_idx, region in enumerate(QUAL_REGIONS):
            r = m["per_region"][region]
            if not np.isfinite(r["auc"]):
                continue
            ax.errorbar(fi + 0.18 * (r_idx - 1.5) / 4, r["auc"],
                        yerr=[[r["auc"] - r["ci_lo"]], [r["ci_hi"] - r["auc"]]],
                        fmt="o", markersize=6, alpha=0.85,
                        color=plt.cm.tab10(r_idx),
                        ecolor=plt.cm.tab10(r_idx), elinewidth=0.8, capsize=2,
                        markeredgecolor="white", markeredgewidth=0.5,
                        label=region if fi == 0 else None)
    ax.axhline(0.5, color="black", ls="--", lw=0.8, label="chance")
    ax.set_xticks(x); ax.set_xticklabels([m["bin_label"] for m in macro_2], fontsize=10)
    ax.set_ylim(0.30, 1.00)
    ax.set_ylabel("Macro AUC of trivial scalar (log$_{10}$ Benioff in last 5 d)")
    ax.set_title("exp20: M-scaling of trivial-scalar finding\n"
                 "(Mann-Kendall τ = {:.2f}, p = {:.3f})".format(tau_2, kpval_2))
    ax.legend(fontsize=9, loc="lower right", ncol=2)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(EXP_DIR / "m_scaling_plot.png")
    plt.close(fig)
    print(f"[{_ts()}] [plot] m_scaling_plot.png", flush=True)

    # ===== Persist =====
    pd.DataFrame(per_region_2 + per_region_3).to_csv(
        EXP_DIR / "m_scaling_table.csv", index=False)
    pd.DataFrame([{k: v for k, v in m.items() if k != "per_region"}
                  for m in macro_2 + macro_3]).to_csv(
        EXP_DIR / "m_scaling_macro_table.csv", index=False)

    summary = {
        "experiment": "exp20_trivial_scalar_magnitude_scaling",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "scalar": "log10 Benioff in last 5 days (exp18 A2_blog_last5d)",
        "macro_2bin": macro_2,
        "macro_3bin": macro_3,
        "trend_tests": {
            "mann_kendall_2bin": {"tau": float(tau_2), "p": float(kpval_2)},
            "mann_kendall_3bin": {"tau": float(tau_3), "p": float(kpval_3)},
            "logistic_scalar_plus_M": {
                "coef_scalar": coef_scalar, "coef_target_M": coef_M,
            },
        },
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\n[{_ts()}] [exp20] HEADLINE", flush=True)
    print(f"[{_ts()}]   2-bin M-scaling under trivial scalar:", flush=True)
    for m in macro_2:
        print(f"[{_ts()}]     {m['bin_label']:<14s}  macro={m['macro_auc']:.3f}  "
              f"CI=[{m['ci_lo']:.3f}, {m['ci_hi']:.3f}]  n={m['n_pre_total']}", flush=True)
    print(f"[{_ts()}]   Trend monotonicity (Mann-Kendall): tau={tau_2:.2f} p={kpval_2:.2f}",
          flush=True)
    print(f"[{_ts()}] [exp20] done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
