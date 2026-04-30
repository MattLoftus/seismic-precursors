"""exp15 — Window-shift control for the exp14 Benioff TLS finding.

exp14 found the Benioff total trajectory passes 3σ + Bonferroni cross-
regionally (macro AUC 0.70, z = +8.01) but the template visualization shows
the discriminative signal is concentrated in sub-window 5 (days 25–30
before target) — i.e., the foreshock period.

This control re-runs the same TLS scan but with the "precursor" window
redefined as [t − 60, t − 30) — shifted 30 days earlier so the window
ends 30 days before the target rather than at the target. If the signal
is a long-distance precursor, the shifted window should show the same
template-correlation effect. If the signal is foreshock-specific, the
shifted window will be null.

This experiment is NOT pre-registered. It is a deviation control,
explicitly to interpret exp14's finding. Reported as exploratory in the
paper.

Outputs:
    shifted_macro_auc_table.csv
    shifted_template_plot.png
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
EXP06_DIR = ROOT / "experiments" / "exp06_cross_regional_macro"

MC_PER_REGION = {
    "California": 3.50,
    "Cascadia": 3.50,
    "Turkey": 3.10,
    "Italy": 3.50,
}
QUAL_REGIONS = ["California", "Cascadia", "Turkey", "Italy"]
NULL_KIND = "null_A"

WINDOW_DAYS = 30
SHIFT_DAYS = 30                # shift the precursor window 30 days earlier
N_SUBWINDOWS = 6

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


def benioff_traj(catalog: pd.DataFrame, t_start, t_end, mc: float,
                 n_subwindows: int) -> np.ndarray:
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
    print(f"[{_ts()}] [exp15] start — window-shift control", flush=True)
    print(f"[{_ts()}] [exp15] precursor window: [t-{WINDOW_DAYS+SHIFT_DAYS}, "
          f"t-{SHIFT_DAYS}) (shifted from exp14's [t-{WINDOW_DAYS}, t))", flush=True)

    df = pd.read_csv(EXP07_CSV)
    df["t_start"] = pd.to_datetime(df["t_start"], utc=True, format="ISO8601")
    df["t_end"] = pd.to_datetime(df["t_end"], utc=True, format="ISO8601")
    df = df[df["region"].isin(QUAL_REGIONS) &
            df["window_kind"].isin(["precursor", NULL_KIND])].copy().reset_index(drop=True)

    # SHIFT only the precursor windows by 30 days earlier; keep nulls as-is
    pre_mask = df["window_kind"] == "precursor"
    df.loc[pre_mask, "t_start"] = df.loc[pre_mask, "t_start"] - pd.Timedelta(days=SHIFT_DAYS)
    df.loc[pre_mask, "t_end"] = df.loc[pre_mask, "t_end"] - pd.Timedelta(days=SHIFT_DAYS)
    print(f"[{_ts()}] [exp15] {pre_mask.sum()} precursor windows shifted; "
          f"{(~pre_mask).sum()} null A windows untouched", flush=True)

    catalogs = {region: load_catalog(region) for region in QUAL_REGIONS}

    # Compute Benioff trajectories
    print(f"[{_ts()}] [exp15] computing Benioff trajectories...", flush=True)
    trajs = []
    for i, row in df.iterrows():
        t = benioff_traj(catalogs[row["region"]], row["t_start"], row["t_end"],
                         MC_PER_REGION[row["region"]], N_SUBWINDOWS)
        trajs.append(t)
    df["benioff_total_traj"] = trajs

    # LORO TLS scan
    print(f"[{_ts()}] [exp15] === LORO TLS scan on shifted window ===", flush=True)
    rng = np.random.default_rng(RANDOM_SEED)
    results = []
    templates = {}
    for held_out in QUAL_REGIONS:
        train = df[(df["region"] != held_out) & (df["window_kind"] == "precursor")]
        if len(train) < 5:
            continue
        template = np.mean(np.stack(train["benioff_total_traj"].values), axis=0)
        templates[held_out] = template

        test = df[df["region"] == held_out]
        scores = np.full(len(test), np.nan)
        for ti, t_row in enumerate(test.itertuples(index=False)):
            traj = t_row.benioff_total_traj
            if np.std(traj) == 0:
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
                boot[bi] = np.nan
                continue
            boot[bi] = roc_auc_score(y_test[samp], scores[samp])
        ci_lo, ci_hi = np.nanpercentile(boot, [2.5, 97.5])

        results.append({
            "held_out_region": held_out,
            "auc": float(auc),
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

    print(f"\n[{_ts()}] [headline] benioff_total_traj (shifted) "
          f"macro={macro:.3f}  CI=[{ci_lo:.3f}, {ci_hi:.3f}]  "
          f"z={z:+.2f}  p={p_two:.2e}  3σ={'Y' if sigma3 else 'n'}  "
          f"Bonf={'Y' if bonf else 'n'}", flush=True)

    # === Plot template (shifted) ===
    fig, ax = plt.subplots(figsize=(7, 4), dpi=120)
    for r_idx, region in enumerate(QUAL_REGIONS):
        if region not in templates:
            continue
        ax.plot(np.arange(N_SUBWINDOWS), templates[region],
                marker="o", lw=1.5, alpha=0.85,
                color=plt.cm.tab10(r_idx), label=f"hold-out {region}")
    ax.set_title("Benioff template (shifted window: [t-60d, t-30d))", fontsize=10)
    ax.set_xlabel("sub-window (each 5d, starting 60d before target)")
    ax.set_ylabel(r"log$_{10}$ $\Sigma\sqrt{E_J}$")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(EXP_DIR / "shifted_template_plot.png")
    plt.close(fig)
    print(f"[{_ts()}] [plot] shifted_template_plot.png", flush=True)

    summary = {
        "experiment": "exp15_window_shift_control",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "purpose": ("Re-run exp14 Benioff TLS with the precursor window shifted "
                    "30 days earlier to t-60..t-30, leaving null windows untouched. "
                    "If signal disappears, exp14's pass was foreshock-specific."),
        "exp14_unshifted_macro_auc": 0.704,
        "exp14_unshifted_z": 8.01,
        "shifted_macro_auc": macro,
        "shifted_ci_lo": float(ci_lo),
        "shifted_ci_hi": float(ci_hi),
        "shifted_z": float(z),
        "shifted_p_two_sided": float(p_two),
        "shifted_passes_3sigma": bool(sigma3),
        "shifted_passes_bonferroni": bool(bonf),
        "per_region_aucs": {r["held_out_region"]: r["auc"] for r in results},
        "interpretation": (
            "If shifted macro AUC ≈ 0.5: the exp14 finding is foreshock-specific "
            "(signal in last sub-window of original window). "
            "If shifted macro AUC ≈ 0.70: the trajectory shape is intrinsic to "
            "precursor periods, not just foreshock leakage."
        ),
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[{_ts()}] [persist] summary.json", flush=True)

    print(f"\n[{_ts()}] [exp15] VERDICT", flush=True)
    if abs(macro - 0.5) < 0.05:
        print(f"[{_ts()}] [exp15]   shifted macro AUC = {macro:.3f} ≈ 0.5 — "
              f"the exp14 finding was FORESHOCK LEAKAGE in the last sub-window, "
              f"not a long-distance precursor signal", flush=True)
    elif macro > 0.6:
        print(f"[{_ts()}] [exp15]   shifted macro AUC = {macro:.3f} > 0.6 — "
              f"trajectory shape persists in shifted window; potential REAL "
              f"long-distance precursor signal", flush=True)
    else:
        print(f"[{_ts()}] [exp15]   shifted macro AUC = {macro:.3f} — partial; "
              f"some signal survives the shift but weakened", flush=True)
    print(f"[{_ts()}] [exp15] done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
