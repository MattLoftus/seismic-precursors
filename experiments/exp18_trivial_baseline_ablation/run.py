"""exp18 — Trivial-baseline ablation for the TLS Benioff finding (peer-review response).

All three Round-E2 cold-read peer reviewers (domain, methods, statistics)
flagged the same concern: the TLS Benioff trajectory scan
(macro AUC=0.704, z=+8.01) may be no more than a slope detector that
trivially picks up the late-window event uptick caused by foreshocks.
Since Benioff strain is monotonically non-decreasing and the precursor
template has a positive slope by construction, a Pearson correlation
against it is essentially "is there a final-week increment?"

This experiment runs the trivial-baseline ablation:
  Baseline A1: count of events in the LAST 5-day sub-window (n[5])
  Baseline A2: log10 cumulative Benioff in the LAST sub-window (b[5])
  Baseline B1: count of events in the WHOLE window (sum n[1..6])
  Baseline B2: total log10 Benioff over the whole window (b[6])
  Baseline C1: increment in count over last sub-window (n[6] - n[5])
  Baseline C2: increment in log10 Benioff over last sub-window (b[6]-b[5])

Each baseline is a SCALAR per window. AUC is computed per LORO held-out
region against {precursor vs Null A} labels, then macro-averaged across
the 4 LORO splits. No template, no Pearson, no degrees of freedom.

If Baseline A1 (or any of these) reaches macro AUC ~0.70, the TLS
machinery in exp14 is ornamental and the paper's central positive is
"there are more events in the last 5 days of a precursor window," which
is the Trugman & Ross 2019 finding cross-regionally. If all baselines
are ≲0.60 while TLS is 0.70, the TLS template adds genuine information.

Outputs:
  baseline_macro_table.csv   — per-baseline macro AUC + per-region AUCs
  baseline_per_split_table.csv — per-baseline × per-LORO-region AUC + n
  baseline_comparison_plot.png — bar chart of all baselines vs TLS=0.704
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

MC_PER_REGION = {
    "California": 3.50,
    "Cascadia": 3.50,
    "Turkey": 3.10,
    "Italy": 3.50,
}
QUAL_REGIONS = ["California", "Cascadia", "Turkey", "Italy"]
NULL_KIND = "null_A"
N_SUBWINDOWS = 6
SUBWINDOW_DAYS = 5

N_BOOT = 1000
N_PERM = 1000
RANDOM_SEED = 42

TLS_REFERENCE_AUC = 0.704  # exp14 macro for benioff_total_traj


def _ts() -> str:
    return dt.datetime.now().strftime("%H:%M:%S")


def load_catalog(region: str) -> pd.DataFrame:
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
    """Return (n_above_mc trajectory, log10 Benioff trajectory) per sub-window."""
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
    b_arr_log = np.log10(np.where(b_arr > 0, b_arr, 1.0))
    return n_arr, b_arr_log


def per_feature_auc(scores: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(scores)
    if mask.sum() < 4 or (y[mask] == 1).sum() < 2 or (y[mask] == 0).sum() < 2:
        return float("nan")
    return float(roc_auc_score(y[mask], scores[mask]))


def macro_with_inference(per_split_aucs: list[dict], rng) -> dict:
    """Compute macro across 4 LORO splits with bootstrap CI95 + permutation z."""
    aucs = np.array([r["auc"] for r in per_split_aucs])
    if not np.any(np.isfinite(aucs)):
        return {"macro": float("nan"), "ci_lo": float("nan"), "ci_hi": float("nan"),
                "perm_z": float("nan"), "perm_p": float("nan"),
                "per_region_aucs": {r["held_out_region"]: r["auc"] for r in per_split_aucs}}

    macro = float(np.nanmean(aucs))

    boot = np.empty(N_BOOT)
    for bi in range(N_BOOT):
        idx = rng.integers(0, len(aucs), size=len(aucs))
        boot[bi] = np.nanmean(aucs[idx])
    ci_lo, ci_hi = np.nanpercentile(boot, [2.5, 97.5])

    # Permutation null: shuffle labels within each held-out region's test set
    perm = np.empty(N_PERM)
    for p_i in range(N_PERM):
        ps = []
        for r in per_split_aucs:
            y = np.array(r["y_test"]); s = np.array(r["scores"])
            m = np.isfinite(s)
            y2 = y[m].copy(); s2 = s[m]
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
    p_two = float(2 * (1 - norm.cdf(abs(z))))
    return {
        "macro": macro, "ci_lo": float(ci_lo), "ci_hi": float(ci_hi),
        "perm_z": float(z), "perm_p": p_two,
        "per_region_aucs": {r["held_out_region"]: r["auc"] for r in per_split_aucs},
    }


def main() -> int:
    print(f"[{_ts()}] [exp18] start — trivial-baseline ablation for TLS Benioff", flush=True)
    rng = np.random.default_rng(RANDOM_SEED)

    df = pd.read_csv(EXP07_CSV)
    df["t_start"] = pd.to_datetime(df["t_start"], utc=True, format="ISO8601")
    df["t_end"] = pd.to_datetime(df["t_end"], utc=True, format="ISO8601")
    df = df[df["region"].isin(QUAL_REGIONS) &
            df["window_kind"].isin(["precursor", NULL_KIND])].copy().reset_index(drop=True)
    print(f"[{_ts()}] [exp18] {len(df)} windows in qualifying regions × {{precursor, {NULL_KIND}}}",
          flush=True)
    print(f"[{_ts()}] [exp18] per-region: " +
          ", ".join(f"{r}: {(df['region']==r).sum()}" for r in QUAL_REGIONS), flush=True)

    catalogs = {region: load_catalog(region) for region in QUAL_REGIONS}
    print(f"[{_ts()}] [exp18] loaded catalogs: " +
          ", ".join(f"{r}={len(c)}" for r, c in catalogs.items()), flush=True)

    # Compute trajectories once
    print(f"[{_ts()}] [exp18] computing 6-point trajectories per window...", flush=True)
    n_traj = np.zeros((len(df), N_SUBWINDOWS))
    b_traj = np.zeros((len(df), N_SUBWINDOWS))
    for i, row in df.iterrows():
        cat = catalogs[row["region"]]
        n, b = catalog_trajectory(cat, row["t_start"], row["t_end"],
                                   MC_PER_REGION[row["region"]])
        n_traj[i] = n
        b_traj[i] = b
        if (i + 1) % 50 == 0:
            print(f"[{_ts()}]   trajectories: {i+1}/{len(df)} done", flush=True)
    print(f"[{_ts()}] [exp18] trajectories complete", flush=True)

    # === Define trivial baselines ===
    # Each baseline maps a window's 6-point trajectory to a single scalar.
    # Index 5 = LAST sub-window (days 25-30). Index 0 = FIRST sub-window (days 0-5).
    baselines = {
        "A1_n_last5d": n_traj[:, 5],
        "A2_blog_last5d": b_traj[:, 5],
        "B1_n_total": n_traj.sum(axis=1),
        "B2_blog_total": b_traj[:, -1],  # last value of cumulative log = ~total
        "C1_n_increment_last": n_traj[:, 5] - n_traj[:, 4],
        "C2_blog_increment_last": b_traj[:, 5] - b_traj[:, 4],
    }
    print(f"[{_ts()}] [exp18] {len(baselines)} baselines defined: {list(baselines.keys())}",
          flush=True)

    # === LORO per-baseline AUC ===
    print(f"\n[{_ts()}] [exp18] === per-baseline LORO scan ===", flush=True)
    all_per_split = []
    for bname, bvec in baselines.items():
        for held_out in QUAL_REGIONS:
            test_mask = df["region"] == held_out
            test_idx = np.where(test_mask)[0]
            if len(test_idx) < 4:
                continue
            scores = bvec[test_idx]
            y_test = (df.iloc[test_idx]["window_kind"] == "precursor").astype(int).to_numpy()
            auc = per_feature_auc(scores, y_test)
            all_per_split.append({
                "baseline": bname, "held_out_region": held_out,
                "n_test_pre": int((y_test == 1).sum()),
                "n_test_null": int((y_test == 0).sum()),
                "auc": float(auc),
                "y_test": y_test.tolist(),
                "scores": scores.tolist(),
            })
            print(f"[{_ts()}]   {bname:<26s}  test={held_out:<11s}  AUC={auc:.3f}  "
                  f"(n_pre={int((y_test==1).sum())}, n_null={int((y_test==0).sum())})", flush=True)

    # === Macro per baseline with full inference ===
    print(f"\n[{_ts()}] [exp18] === macro AUC per baseline ===", flush=True)
    macro_table = []
    for bname in baselines:
        feat_results = [r for r in all_per_split if r["baseline"] == bname]
        m = macro_with_inference(feat_results, rng)
        m["baseline"] = bname
        macro_table.append(m)
        print(f"[{_ts()}] [headline] {bname:<26s}  macro={m['macro']:.3f}  "
              f"CI=[{m['ci_lo']:.3f}, {m['ci_hi']:.3f}]  z={m['perm_z']:+.2f}  "
              f"p={m['perm_p']:.2e}", flush=True)

    pd.DataFrame([{k: v for k, v in r.items() if k not in ("y_test", "scores")}
                  for r in all_per_split]).to_csv(
        EXP_DIR / "baseline_per_split_table.csv", index=False)
    pd.DataFrame([{k: v for k, v in m.items() if k != "per_region_aucs"}
                  for m in macro_table]).to_csv(
        EXP_DIR / "baseline_macro_table.csv", index=False)

    # === Comparison plot vs TLS reference (exp14 = 0.704) ===
    fig, ax = plt.subplots(figsize=(10, 5), dpi=120)
    names = [m["baseline"] for m in macro_table] + ["TLS_benioff_template (exp14)"]
    macros = [m["macro"] for m in macro_table] + [TLS_REFERENCE_AUC]
    err_lo = [m["macro"] - m["ci_lo"] for m in macro_table] + [0.064]  # exp14 CI [0.64, 0.77]
    err_hi = [m["ci_hi"] - m["macro"] for m in macro_table] + [0.066]
    colors = ["#7d8aa6"] * len(macro_table) + ["#cc4444"]
    x = np.arange(len(names))
    ax.bar(x, macros, color=colors, alpha=0.6, yerr=[err_lo, err_hi],
           capsize=5, edgecolor="black")
    for fi, m in enumerate(macro_table):
        for r_idx, region in enumerate(QUAL_REGIONS):
            if region not in m["per_region_aucs"]:
                continue
            ax.scatter(fi + 0.18 * (r_idx - 1.5) / 4, m["per_region_aucs"][region],
                       s=22, alpha=0.85, color=plt.cm.tab10(r_idx),
                       edgecolor="white", linewidth=0.5,
                       label=region if fi == 0 else None)
    ax.axhline(0.5, color="black", ls="--", lw=0.8, label="chance")
    ax.axhline(TLS_REFERENCE_AUC, color="#cc4444", ls=":", lw=1.2,
               label=f"TLS exp14 = {TLS_REFERENCE_AUC}")
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=8, rotation=20, ha="right")
    ax.set_ylim(0.30, 0.85)
    ax.set_ylabel("Macro AUC (LORO across 4 regions)")
    ax.set_title("exp18: trivial-baseline ablation — does TLS template add anything?")
    ax.legend(fontsize=8, loc="lower right", ncol=2)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(EXP_DIR / "baseline_comparison_plot.png")
    plt.close(fig)
    print(f"[{_ts()}] [plot] baseline_comparison_plot.png", flush=True)

    # === Verdict ===
    best_baseline = max(macro_table, key=lambda m: m["macro"]
                        if np.isfinite(m["macro"]) else -1)
    delta_vs_tls = TLS_REFERENCE_AUC - best_baseline["macro"]
    if delta_vs_tls < 0.02:
        verdict = ("TLS_ORNAMENTAL — best trivial baseline matches TLS within 0.02; "
                   "TLS template adds no information")
    elif delta_vs_tls < 0.05:
        verdict = ("TLS_MARGINAL — TLS exceeds best baseline by 0.02-0.05; "
                   "weak evidence the template adds anything")
    else:
        verdict = ("TLS_SUBSTANTIVE — TLS exceeds best baseline by ≥0.05; "
                   "the template demonstrably adds information")

    summary = {
        "experiment": "exp18_trivial_baseline_ablation",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "tls_reference_auc_exp14": TLS_REFERENCE_AUC,
        "baselines_evaluated": list(baselines.keys()),
        "best_baseline": best_baseline["baseline"],
        "best_baseline_macro": best_baseline["macro"],
        "best_baseline_per_region": best_baseline["per_region_aucs"],
        "delta_tls_minus_best_baseline": float(delta_vs_tls),
        "verdict": verdict,
        "macro_table": [{k: v for k, v in m.items()} for m in macro_table],
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[{_ts()}] [persist] summary.json", flush=True)

    print(f"\n[{_ts()}] [exp18] HEADLINE", flush=True)
    print(f"[{_ts()}]   TLS exp14 reference macro AUC: {TLS_REFERENCE_AUC}", flush=True)
    print(f"[{_ts()}]   Best trivial baseline:        {best_baseline['baseline']} = "
          f"{best_baseline['macro']:.3f}", flush=True)
    print(f"[{_ts()}]   Δ (TLS − best baseline):      {delta_vs_tls:+.3f}", flush=True)
    print(f"[{_ts()}]   VERDICT: {verdict}", flush=True)
    print(f"[{_ts()}] [exp18] done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
