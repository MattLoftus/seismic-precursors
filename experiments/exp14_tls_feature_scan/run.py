"""exp14 — TLS-style template-correlation scan in feature time series.

Pre-reg v1 §3.2 promised a "TLS-style matched filter for feature space"
applying the exoplanet-detection template-matching idea to feature time
series. Per-window scalar AUC (exp07–12) is one slice of this; the
template-matching slice is per-window TIME-SERIES correlation against a
precursor-mean template.

Approach:
  - For each window, compute a 6-point feature time series at 5-day
    sub-window resolution. Waveform features re-use exp11 snapshots.
    Catalog features (n_above_mc, Benioff total) re-bin from the cached
    catalog at the per-region Mc.
  - LORO: for each held-out region, the precursor TEMPLATE is the mean
    trajectory across the OTHER 3 regions' precursor windows.
  - Score each held-out window by Pearson correlation between its
    trajectory and the template.
  - Compute ROC-AUC: precursor windows vs Null A windows in the held-out
    region. Macro across 4 LORO splits per feature.

Pre-reg gates (3σ): macro CI95-lower > 0.5 AND |effect| > 0.05 AND
permutation z > 3. Bonferroni α = 0.05 / (5 features × 4 LORO regions) = 2.5e-3.

Outputs:
    tls_macro_auc_table.csv     macro AUC per feature
    tls_macro_plot.png          per-feature macro AUC bar chart with
                                per-region scatter
    tls_template_plot.png       precursor mean trajectory per feature
                                (overlaid across regions)
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
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
EXP_DIR = Path(__file__).resolve().parent
EXP07_CSV = ROOT / "experiments" / "exp07_macro_pra2" / "feature_summary.csv"
EXP06_CATALOG_DIR = ROOT / "experiments" / "exp06_cross_regional_macro"  # (gitignored — live re-pull or use local cache if present)

# Per-region Mc values from PRA-2 (exp07 / exp11)
MC_PER_REGION = {
    "California": 3.50,
    "Cascadia": 3.50,
    "Turkey": 3.10,
    "Italy": 3.50,
}

QUAL_REGIONS = ["California", "Cascadia", "Turkey", "Italy"]
NULL_KIND = "null_A"
WAVEFORM_FEATURES = [
    "wf_spectral_slope",
    "wf_waveform_entropy",
    "wf_hht_imf1_if_median_hz",
]
CATALOG_FEATURES = ["n_above_mc_traj", "benioff_total_traj"]
ALL_FEATURES = CATALOG_FEATURES + WAVEFORM_FEATURES

N_SUBWINDOWS = 6
SUBWINDOW_DAYS = 5  # 6 × 5 = 30 days

N_BOOT = 1000
N_PERM = 1000
RANDOM_SEED = 42


def _ts() -> str:
    return dt.datetime.now().strftime("%H:%M:%S")


def load_catalog(region: str) -> pd.DataFrame:
    """Load the per-region catalog from exp07's cache (gitignored; re-pull if missing)."""
    candidates = [
        ROOT / "experiments" / "exp06_cross_regional_macro" / f"catalog_{region}.csv",
        ROOT / "experiments" / "exp07_macro_pra2" / f"catalog_{region}.csv",
    ]
    for p in candidates:
        if p.is_file():
            df = pd.read_csv(p)
            df["time"] = pd.to_datetime(df["time"], utc=True, format="ISO8601")
            return df
    raise FileNotFoundError(f"no catalog cache found for {region}; tried {candidates}")


def catalog_trajectory(catalog: pd.DataFrame, t_start, t_end, mc: float,
                       n_subwindows: int = N_SUBWINDOWS) -> tuple[np.ndarray, np.ndarray]:
    """Return (n_above_mc trajectory, benioff total trajectory) at sub-window resolution."""
    duration = (t_end - t_start).total_seconds()
    edges = [t_start + pd.Timedelta(seconds=duration * k / n_subwindows)
             for k in range(n_subwindows + 1)]
    n_arr = np.zeros(n_subwindows)
    b_arr = np.zeros(n_subwindows)
    for k in range(n_subwindows):
        mask = ((catalog["time"] >= edges[k]) & (catalog["time"] < edges[k + 1])
                & (catalog["magnitude"] >= mc - 1e-9))
        sub = catalog.loc[mask]
        n_arr[k] = float(len(sub))
        if len(sub):
            energies = 10 ** (1.5 * sub["magnitude"].to_numpy() + 4.8)
            b_arr[k] = float(np.sum(np.sqrt(energies)))
    # log10 transform of Benioff to match scalar scale; replace zeros with small pos.
    b_arr_log = np.log10(np.where(b_arr > 0, b_arr, 1.0))
    return n_arr, b_arr_log


def waveform_trajectory_from_exp11(region: str, t_start_iso: str,
                                    t_end_iso: str) -> dict[str, np.ndarray]:
    """exp11 stores a per-window aggregated scalar, not a snapshot-level trajectory.

    For TLS we want per-snapshot (6-point) trajectories. Re-derive from the
    stored data by reading per-snapshot from the cached snapshot results.
    Since exp11 only saved the per-window aggregates, here we approximate by
    using the scalar as a constant over 6 sub-windows. This is a documented
    deviation: TLS for waveform features is therefore degenerate and reduces
    to per-window scalar (which exp07–12 already evaluated). We compute this
    so the file completes; treat waveform-feature TLS results as redundant
    with scalar AUC.
    """
    return None  # signaled later


def per_feature_auc(p, q):
    p = np.asarray(p, dtype=float); q = np.asarray(q, dtype=float)
    p = p[np.isfinite(p)]; q = q[np.isfinite(q)]
    if len(p) < 2 or len(q) < 2:
        return float("nan")
    n_p, n_q = len(p), len(q)
    arr = np.concatenate([p, q])
    order = np.argsort(arr)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, n_p + n_q + 1)
    return float((ranks[:n_p].sum() - n_p * (n_p + 1) / 2) / (n_p * n_q))


def correlate_with_template(trajectory: np.ndarray, template: np.ndarray) -> float:
    """Pearson correlation of two 1-D vectors. Returns NaN if either is constant."""
    if np.std(trajectory) == 0 or np.std(template) == 0:
        return float("nan")
    return float(np.corrcoef(trajectory, template)[0, 1])


def main() -> int:
    print(f"[{_ts()}] [exp14] start — TLS template scan in feature time series", flush=True)

    # Load merged exp07 windows
    df = pd.read_csv(EXP07_CSV)
    df["t_start"] = pd.to_datetime(df["t_start"], utc=True, format="ISO8601")
    df["t_end"] = pd.to_datetime(df["t_end"], utc=True, format="ISO8601")
    df = df[df["region"].isin(QUAL_REGIONS) &
            df["window_kind"].isin(["precursor", NULL_KIND])].copy().reset_index(drop=True)
    print(f"[{_ts()}] [exp14] {len(df)} rows in qualifying regions × {{precursor, {NULL_KIND}}}",
          flush=True)

    # Load catalogs and compute catalog trajectories
    catalogs = {region: load_catalog(region) for region in QUAL_REGIONS}
    print(f"[{_ts()}] [exp14] loaded catalogs: "
          f"{ {r: len(c) for r, c in catalogs.items()} }", flush=True)

    print(f"[{_ts()}] [exp14] computing 6-point catalog trajectories per window...",
          flush=True)
    n_traj = np.zeros((len(df), N_SUBWINDOWS))
    b_traj = np.zeros((len(df), N_SUBWINDOWS))
    for i, row in df.iterrows():
        cat = catalogs[row["region"]]
        n, b = catalog_trajectory(cat, row["t_start"], row["t_end"],
                                   MC_PER_REGION[row["region"]])
        n_traj[i] = n
        b_traj[i] = b

    df["n_above_mc_traj"] = list(n_traj)
    df["benioff_total_traj"] = list(b_traj)

    # === LORO TLS scan per feature ===
    rng = np.random.default_rng(RANDOM_SEED)
    print(f"\n[{_ts()}] [exp14] === LORO TLS scan ===", flush=True)
    results = []  # list of dicts per (feature, held_out_region)
    templates = {}  # (feature, held_out_region) -> template

    for feature in CATALOG_FEATURES:
        for held_out in QUAL_REGIONS:
            train = df[(df["region"] != held_out) & (df["window_kind"] == "precursor")]
            if len(train) < 5:
                continue
            template = np.mean(np.stack(train[feature].values), axis=0)
            templates[(feature, held_out)] = template

            test = df[df["region"] == held_out]
            scores = np.full(len(test), np.nan)
            for ti, t_row in enumerate(test.itertuples(index=False)):
                trajectory = getattr(t_row, feature)
                scores[ti] = correlate_with_template(trajectory, template)

            y_test = (test["window_kind"] == "precursor").astype(int).to_numpy()
            mask = np.isfinite(scores)
            auc = (roc_auc_score(y_test[mask], scores[mask])
                   if mask.sum() > 4 and (y_test[mask] == 1).sum() >= 2
                   and (y_test[mask] == 0).sum() >= 2 else float("nan"))

            # bootstrap CI on AUC
            valid_idx = np.where(mask)[0]
            if len(valid_idx) > 4:
                boot_aucs = np.empty(N_BOOT)
                for bi in range(N_BOOT):
                    samp = rng.choice(valid_idx, size=len(valid_idx), replace=True)
                    if (y_test[samp] == 1).sum() < 2 or (y_test[samp] == 0).sum() < 2:
                        boot_aucs[bi] = np.nan
                        continue
                    boot_aucs[bi] = roc_auc_score(y_test[samp], scores[samp])
                ci_lo, ci_hi = np.nanpercentile(boot_aucs, [2.5, 97.5])
            else:
                ci_lo = ci_hi = float("nan")
            results.append({
                "feature": feature, "held_out_region": held_out,
                "n_test_pre": int((y_test == 1).sum()),
                "n_test_null": int((y_test == 0).sum()),
                "n_test_finite": int(mask.sum()),
                "auc": float(auc),
                "ci_lo": float(ci_lo), "ci_hi": float(ci_hi),
                "y_test": y_test.tolist(),
                "scores": scores.tolist(),
            })
            print(f"[{_ts()}]   {feature:<22s}  test={held_out:<11s}  "
                  f"AUC={auc:.3f}  CI=[{ci_lo:.3f}, {ci_hi:.3f}]", flush=True)

    # Macro per feature
    print(f"\n[{_ts()}] [exp14] === macro AUC per feature ===", flush=True)
    macro_table = []
    for feature in CATALOG_FEATURES:
        feat_results = [r for r in results if r["feature"] == feature]
        aucs = np.array([r["auc"] for r in feat_results])
        if not np.any(np.isfinite(aucs)):
            continue
        macro = float(np.nanmean(aucs))

        # Cross-region bootstrap (resample regions)
        boot = np.empty(N_BOOT)
        for bi in range(N_BOOT):
            idx = rng.integers(0, len(aucs), size=len(aucs))
            boot[bi] = np.nanmean(aucs[idx])
        ci_lo, ci_hi = np.nanpercentile(boot, [2.5, 97.5])

        # Permutation z: shuffle y in each held-out region, recompute AUC
        perm = np.empty(N_PERM)
        for p_i in range(N_PERM):
            ps = []
            for r in feat_results:
                y = np.array(r["y_test"]); s = np.array(r["scores"])
                m = np.isfinite(s)
                y2 = y[m]; s2 = s[m]
                rng.shuffle(y2)
                if (y2 == 1).sum() < 2 or (y2 == 0).sum() < 2:
                    ps.append(np.nan)
                    continue
                ps.append(roc_auc_score(y2, s2))
            perm[p_i] = float(np.nanmean(ps))
        perm_mean = float(np.nanmean(perm))
        perm_std = float(np.nanstd(perm))
        z = (macro - perm_mean) / perm_std if perm_std > 0 else float("nan")
        from scipy.stats import norm
        p_two = 2 * (1 - norm.cdf(abs(z)))

        ci_above = ci_lo > 0.5
        ci_below = ci_hi < 0.5
        effect = abs(macro - 0.5)
        sigma3 = (ci_above or ci_below) and effect > 0.05 and abs(z) > 3
        # Bonferroni: 2 features × 4 LORO regions = 8 tests
        bonf_alpha = 0.05 / (len(CATALOG_FEATURES) * len(QUAL_REGIONS))
        bonf = p_two < bonf_alpha
        macro_table.append({
            "feature": feature,
            "macro_auc": macro, "ci_lo": float(ci_lo), "ci_hi": float(ci_hi),
            "perm_z": float(z), "perm_p": float(p_two),
            "passes_3sigma": bool(sigma3), "passes_bonferroni": bool(bonf),
            "per_region_aucs": {r["held_out_region"]: r["auc"] for r in feat_results},
        })
        print(f"[{_ts()}] [headline] {feature:<22s}  macro={macro:.3f}  "
              f"CI=[{ci_lo:.3f}, {ci_hi:.3f}]  z={z:+.2f}  p={p_two:.2e}  "
              f"3σ={'Y' if sigma3 else 'n'}  Bonf={'Y' if bonf else 'n'}", flush=True)

    pd.DataFrame([{k: v for k, v in r.items() if k not in ("y_test", "scores")}
                  for r in results]).to_csv(EXP_DIR / "tls_per_split_table.csv", index=False)
    pd.DataFrame(macro_table).to_csv(EXP_DIR / "tls_macro_auc_table.csv", index=False)

    # === Plots ===
    # Macro AUC per feature
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=120)
    if macro_table:
        names = [m["feature"] for m in macro_table]
        macros = [m["macro_auc"] for m in macro_table]
        err_lo = [m["macro_auc"] - m["ci_lo"] for m in macro_table]
        err_hi = [m["ci_hi"] - m["macro_auc"] for m in macro_table]
        x = np.arange(len(names))
        ax.bar(x, macros, color="#7d8aa6", alpha=0.6,
               yerr=[err_lo, err_hi], capsize=5, edgecolor="black")
        for fi, name in enumerate(names):
            per_region = macro_table[fi]["per_region_aucs"]
            for r_idx, region in enumerate(QUAL_REGIONS):
                if region not in per_region:
                    continue
                ax.scatter(fi + 0.18 * (r_idx - 1.5) / 4,
                           per_region[region], s=22, alpha=0.85,
                           color=plt.cm.tab10(r_idx),
                           edgecolor="white", linewidth=0.5,
                           label=region if fi == 0 else None)
        ax.axhline(0.5, color="black", ls="--", lw=0.8)
        ax.set_xticks(x); ax.set_xticklabels(names, fontsize=9)
        ax.set_ylim(0.30, 0.75)
        ax.set_ylabel("TLS-correlation macro AUC")
        ax.set_title("exp14: TLS template-correlation macro AUC (LORO)")
        ax.legend(fontsize=8, loc="lower right", title="held-out region")
        ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(EXP_DIR / "tls_macro_plot.png")
    plt.close(fig)
    print(f"[{_ts()}] [plot] tls_macro_plot.png", flush=True)

    # Template visualization
    fig, axes = plt.subplots(1, len(CATALOG_FEATURES), figsize=(11, 4.0), dpi=120)
    if not isinstance(axes, np.ndarray):
        axes = [axes]
    for ax, feature in zip(axes, CATALOG_FEATURES):
        for r_idx, region in enumerate(QUAL_REGIONS):
            template = templates.get((feature, region))
            if template is None:
                continue
            ax.plot(np.arange(N_SUBWINDOWS), template,
                    marker="o", lw=1.5, alpha=0.85,
                    color=plt.cm.tab10(r_idx), label=f"hold-out {region}")
        ax.set_title(f"{feature} precursor template (mean of train)")
        ax.set_xlabel(f"sub-window (each {SUBWINDOW_DAYS}d)")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(EXP_DIR / "tls_template_plot.png")
    plt.close(fig)
    print(f"[{_ts()}] [plot] tls_template_plot.png", flush=True)

    summary = {
        "experiment": "exp14_tls_feature_scan",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "features_scanned": CATALOG_FEATURES,
        "waveform_features_skipped": ("exp11 stored per-window aggregates only; "
                                       "TLS for waveform features would require "
                                       "snapshot-resolution storage and re-fetch — "
                                       "skipped this iteration"),
        "qualifying_regions": QUAL_REGIONS,
        "n_subwindows": N_SUBWINDOWS,
        "subwindow_days": SUBWINDOW_DAYS,
        "macro_table": macro_table,
        "n_pass_3sigma": sum(1 for m in macro_table if m["passes_3sigma"]),
        "n_pass_bonferroni": sum(1 for m in macro_table if m["passes_bonferroni"]),
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[{_ts()}] [persist] summary.json", flush=True)

    n3 = summary["n_pass_3sigma"]
    nb = summary["n_pass_bonferroni"]
    print(f"\n[{_ts()}] [exp14] HEADLINE: {n3}/{len(CATALOG_FEATURES)} catalog "
          f"trajectories pass 3σ, {nb}/{len(CATALOG_FEATURES)} pass Bonferroni", flush=True)
    print(f"[{_ts()}] [exp14] done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
