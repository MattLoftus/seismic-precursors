"""exp16 — Mask-final-5-days control: does any 25-day-ahead signal exist?

exp14 found Benioff total trajectory passes 3σ + Bonferroni (macro AUC =
0.704), and exp15 confirmed by window-shift that the signal is foreshock-
specific. The natural follow-up: if we EXCLUDE the foreshock sub-window
from the trajectory comparison, does any residual signal survive?

This is the strongest reviewer-defense experiment: it directly tests
whether the result is purely foreshock-driven or has a 25-day-ahead
component.

Approach: same TLS LORO scan as exp14, but the trajectory is only the
first 5 sub-windows (days 0–25 of the [t-30, t) window). Sub-window 5
(days 25–30) is masked out. Templates and scores are computed on the
5-point trajectory.

If the masked AUC ≈ 0.5: the entire exp14 signal lives in days 25–30
(foreshock-exclusive). If the masked AUC > 0.55 with significance: a
25-day-ahead signal exists and the result extends beyond foreshocks.

Outputs:
    masked_macro_auc.json
    masked_template_plot.png
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
EXP06_DIR = ROOT / "experiments" / "exp06_cross_regional_macro"

MC_PER_REGION = {"California": 3.50, "Cascadia": 3.50, "Turkey": 3.10, "Italy": 3.50}
QUAL_REGIONS = ["California", "Cascadia", "Turkey", "Italy"]
NULL_KIND = "null_A"
N_SUBWINDOWS_FULL = 6
KEEP_SUBWINDOWS = 5  # drop the last (foreshock) sub-window
N_BOOT = 1000
N_PERM = 1000
RANDOM_SEED = 42


def _ts() -> str:
    return dt.datetime.now().strftime("%H:%M:%S")


def load_catalog(region: str) -> pd.DataFrame:
    for d in (EXP06_DIR, ROOT / "experiments" / "exp07_macro_pra2"):
        p = d / f"catalog_{region}.csv"
        if p.is_file():
            df = pd.read_csv(p)
            df["time"] = pd.to_datetime(df["time"], utc=True, format="ISO8601")
            return df
    raise FileNotFoundError(f"no catalog for {region}")


def benioff_traj(catalog, t_start, t_end, mc, n_subwindows):
    duration = (t_end - t_start).total_seconds()
    edges = [t_start + pd.Timedelta(seconds=duration * k / n_subwindows)
             for k in range(n_subwindows + 1)]
    out = np.zeros(n_subwindows)
    for k in range(n_subwindows):
        mask = ((catalog["time"] >= edges[k]) & (catalog["time"] < edges[k + 1])
                & (catalog["magnitude"] >= mc - 1e-9))
        sub = catalog.loc[mask]
        if len(sub):
            energies = 10 ** (1.5 * sub["magnitude"].to_numpy() + 4.8)
            out[k] = float(np.sum(np.sqrt(energies)))
    return np.log10(np.where(out > 0, out, 1.0))


def main() -> int:
    print(f"[{_ts()}] [exp16] start — mask final sub-window (days 25-30) from "
          f"Benioff TLS scan", flush=True)

    df = pd.read_csv(EXP07_CSV)
    df["t_start"] = pd.to_datetime(df["t_start"], utc=True, format="ISO8601")
    df["t_end"] = pd.to_datetime(df["t_end"], utc=True, format="ISO8601")
    df = df[df["region"].isin(QUAL_REGIONS) &
            df["window_kind"].isin(["precursor", NULL_KIND])].copy().reset_index(drop=True)

    catalogs = {region: load_catalog(region) for region in QUAL_REGIONS}

    # Compute full 6-point trajectories, then keep only first 5
    print(f"[{_ts()}] [exp16] computing 6-point trajectories then masking sub-window 5...",
          flush=True)
    trajs_full = []
    for i, row in df.iterrows():
        full = benioff_traj(catalogs[row["region"]], row["t_start"], row["t_end"],
                            MC_PER_REGION[row["region"]], N_SUBWINDOWS_FULL)
        trajs_full.append(full)
    df["benioff_full"] = trajs_full
    df["benioff_masked"] = [t[:KEEP_SUBWINDOWS] for t in trajs_full]

    # === LORO TLS on masked trajectory ===
    rng = np.random.default_rng(RANDOM_SEED)
    print(f"\n[{_ts()}] [exp16] === LORO TLS on masked trajectory (5 sub-windows) ===",
          flush=True)
    results = []
    templates = {}
    for held_out in QUAL_REGIONS:
        train = df[(df["region"] != held_out) & (df["window_kind"] == "precursor")]
        if len(train) < 5:
            continue
        template = np.mean(np.stack(train["benioff_masked"].values), axis=0)
        templates[held_out] = template

        test = df[df["region"] == held_out]
        scores = np.full(len(test), np.nan)
        for ti, t_row in enumerate(test.itertuples(index=False)):
            traj = t_row.benioff_masked
            if np.std(traj) == 0 or np.std(template) == 0:
                continue
            scores[ti] = float(np.corrcoef(traj, template)[0, 1])

        y_test = (test["window_kind"] == "precursor").astype(int).to_numpy()
        mask = np.isfinite(scores)
        auc = roc_auc_score(y_test[mask], scores[mask]) if mask.sum() > 4 else float("nan")

        valid = np.where(mask)[0]
        boot = np.empty(N_BOOT)
        for bi in range(N_BOOT):
            samp = rng.choice(valid, size=len(valid), replace=True)
            if (y_test[samp] == 1).sum() < 2 or (y_test[samp] == 0).sum() < 2:
                boot[bi] = np.nan; continue
            boot[bi] = roc_auc_score(y_test[samp], scores[samp])
        ci_lo, ci_hi = np.nanpercentile(boot, [2.5, 97.5])

        results.append({
            "held_out_region": held_out, "auc": float(auc),
            "ci_lo": float(ci_lo), "ci_hi": float(ci_hi),
            "y_test": y_test.tolist(), "scores": scores.tolist(),
        })
        print(f"[{_ts()}]   test={held_out:<11s}  AUC={auc:.3f}  "
              f"CI=[{ci_lo:.3f}, {ci_hi:.3f}]", flush=True)

    aucs = np.array([r["auc"] for r in results])
    macro = float(np.nanmean(aucs))

    boot = np.empty(N_BOOT)
    for bi in range(N_BOOT):
        idx = rng.integers(0, len(aucs), size=len(aucs))
        boot[bi] = np.nanmean(aucs[idx])
    ci_lo, ci_hi = np.nanpercentile(boot, [2.5, 97.5])

    perm = np.empty(N_PERM)
    for p_i in range(N_PERM):
        ps = []
        for r in results:
            y = np.array(r["y_test"]); s = np.array(r["scores"])
            m = np.isfinite(s)
            y2 = y[m]; s2 = s[m]
            rng.shuffle(y2)
            if (y2 == 1).sum() < 2 or (y2 == 0).sum() < 2:
                ps.append(np.nan); continue
            ps.append(roc_auc_score(y2, s2))
        perm[p_i] = float(np.nanmean(ps))
    perm_mean = float(np.nanmean(perm))
    perm_std = float(np.nanstd(perm))
    z = (macro - perm_mean) / perm_std if perm_std > 0 else float("nan")
    from scipy.stats import norm
    p_two = 2 * (1 - norm.cdf(abs(z)))

    ci_above = ci_lo > 0.5
    effect = abs(macro - 0.5)
    sigma3 = ci_above and effect > 0.05 and abs(z) > 3
    bonf = p_two < 0.05 / 4

    print(f"\n[{_ts()}] [headline] benioff_total_traj (masked: days 0-25 only)  "
          f"macro={macro:.3f}  CI=[{ci_lo:.3f}, {ci_hi:.3f}]  z={z:+.2f}  "
          f"p={p_two:.2e}  3σ={'Y' if sigma3 else 'n'}  "
          f"Bonf={'Y' if bonf else 'n'}", flush=True)

    # Comparison: full (exp14) vs masked (exp16)
    print(f"\n[{_ts()}] [comparison]", flush=True)
    print(f"[{_ts()}]   exp14 (days 0-30)   macro=0.704  z=+8.01", flush=True)
    print(f"[{_ts()}]   exp16 (days 0-25)   macro={macro:.3f}  z={z:+.2f}", flush=True)
    if abs(macro - 0.5) < 0.03:
        verdict = ("ENTIRELY FORESHOCK-DRIVEN: signal collapses to chance when "
                   "days 25-30 are masked")
    elif macro > 0.55 and sigma3:
        verdict = ("RESIDUAL 25-DAY-AHEAD SIGNAL EXISTS: macro AUC > 0.55 with "
                   "3σ even after foreshock masking")
    elif macro > 0.55:
        verdict = ("WEAK 25-DAY-AHEAD SIGNAL: macro AUC > 0.55 but doesn't pass 3σ")
    else:
        verdict = ("MOSTLY FORESHOCK-DRIVEN: small residual signal but inconclusive")
    print(f"[{_ts()}] [verdict] {verdict}", flush=True)

    # Plot template (5-point)
    fig, ax = plt.subplots(figsize=(7, 4), dpi=120)
    for r_idx, region in enumerate(QUAL_REGIONS):
        if region not in templates:
            continue
        ax.plot(np.arange(KEEP_SUBWINDOWS), templates[region],
                marker="o", lw=1.5, alpha=0.85,
                color=plt.cm.tab10(r_idx), label=f"hold-out {region}")
    ax.set_title("Benioff template, masked (days 0-25, foreshock excluded)", fontsize=10)
    ax.set_xlabel("sub-window (each 5d)"); ax.set_ylabel(r"log$_{10}$ $\Sigma\sqrt{E_J}$")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(EXP_DIR / "masked_template_plot.png")
    plt.close(fig)
    print(f"[{_ts()}] [plot] masked_template_plot.png", flush=True)

    summary = {
        "experiment": "exp16_mask_final5_control",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "purpose": ("Re-run exp14 Benioff TLS with sub-window 5 (days 25-30) "
                    "masked out. Tests whether any signal exists in days 0-25."),
        "exp14_unmasked_macro_auc": 0.704,
        "exp14_unmasked_z": 8.01,
        "exp16_masked_macro_auc": macro,
        "exp16_masked_ci": [float(ci_lo), float(ci_hi)],
        "exp16_masked_z": float(z),
        "exp16_masked_p_two": float(p_two),
        "exp16_masked_passes_3sigma": bool(sigma3),
        "exp16_masked_passes_bonferroni": bool(bonf),
        "per_region_aucs": {r["held_out_region"]: r["auc"] for r in results},
        "verdict": verdict,
    }
    with open(EXP_DIR / "masked_macro_auc.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[{_ts()}] [persist] masked_macro_auc.json", flush=True)
    print(f"[{_ts()}] [exp16] done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
