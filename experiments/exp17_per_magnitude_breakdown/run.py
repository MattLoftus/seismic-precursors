"""exp17 — Per-target-magnitude breakdown of exp14 Benioff TLS result.

Stratifies the precursor windows by target_M into two bins (M=4.5–5.0 vs
M=5.0+) and re-computes the LORO TLS scan to test whether the Benioff
trajectory signal scales with target magnitude. Physically motivated:
larger mainshocks should have stronger foreshock activity per Helmstetter
et al. 2003 ETAS-derived expectations.

Outputs:
    per_mag_table.csv
    per_mag_plot.png
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

MC_PER_REGION = {"California": 3.50, "Cascadia": 3.50, "Turkey": 3.10, "Italy": 3.50}
QUAL_REGIONS = ["California", "Cascadia", "Turkey", "Italy"]
NULL_KIND = "null_A"
N_SUBWINDOWS = 6
MAG_BINS = [(4.5, 5.0), (5.0, 10.0)]


def _ts() -> str:
    return dt.datetime.now().strftime("%H:%M:%S")


def load_catalog(region):
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
        m = ((catalog["time"] >= edges[k]) & (catalog["time"] < edges[k + 1])
             & (catalog["magnitude"] >= mc - 1e-9))
        sub = catalog.loc[m]
        if len(sub):
            energies = 10 ** (1.5 * sub["magnitude"].to_numpy() + 4.8)
            out[k] = float(np.sum(np.sqrt(energies)))
    return np.log10(np.where(out > 0, out, 1.0))


def main() -> int:
    print(f"[{_ts()}] [exp17] start — per-target-magnitude breakdown of "
          f"exp14 Benioff TLS", flush=True)

    df = pd.read_csv(EXP07_CSV)
    df["t_start"] = pd.to_datetime(df["t_start"], utc=True, format="ISO8601")
    df["t_end"] = pd.to_datetime(df["t_end"], utc=True, format="ISO8601")
    df = df[df["region"].isin(QUAL_REGIONS) &
            df["window_kind"].isin(["precursor", NULL_KIND])].copy().reset_index(drop=True)

    catalogs = {region: load_catalog(region) for region in QUAL_REGIONS}

    print(f"[{_ts()}] [exp17] computing Benioff trajectories...", flush=True)
    trajs = []
    for i, row in df.iterrows():
        t = benioff_traj(catalogs[row["region"]], row["t_start"], row["t_end"],
                         MC_PER_REGION[row["region"]], N_SUBWINDOWS)
        trajs.append(t)
    df["benioff_traj"] = trajs

    # Stratify precursor by target_M
    pre = df[df["window_kind"] == "precursor"]
    print(f"[{_ts()}] [exp17] precursor target_M distribution:", flush=True)
    print(f"[{_ts()}]   mean = {pre['target_M'].mean():.2f}, "
          f"median = {pre['target_M'].median():.2f}, "
          f"max = {pre['target_M'].max():.2f}", flush=True)
    for lo, hi in MAG_BINS:
        n = ((pre['target_M'] >= lo) & (pre['target_M'] < hi)).sum()
        print(f"[{_ts()}]   M=[{lo:.1f}, {hi:.1f}): n={n}", flush=True)

    # === LORO TLS scan per magnitude bin ===
    rows = []
    for lo, hi in MAG_BINS:
        bin_label = f"M[{lo:.1f},{hi:.1f})"
        # Restrict precursor to this bin; null A unchanged
        bin_pre = df[(df["window_kind"] == "precursor")
                      & (df["target_M"] >= lo) & (df["target_M"] < hi)]
        bin_null = df[df["window_kind"] == NULL_KIND]
        bin_df = pd.concat([bin_pre, bin_null], ignore_index=True)

        print(f"\n[{_ts()}] [exp17] === {bin_label}: "
              f"{len(bin_pre)} precursor, {len(bin_null)} null A ===", flush=True)
        if len(bin_pre) < 8:
            print(f"[{_ts()}] [exp17]   too few precursors; skip", flush=True)
            continue

        per_region_aucs = {}
        for held_out in QUAL_REGIONS:
            train = bin_df[(bin_df["region"] != held_out)
                           & (bin_df["window_kind"] == "precursor")]
            if len(train) < 5:
                per_region_aucs[held_out] = float("nan")
                continue
            template = np.mean(np.stack(train["benioff_traj"].values), axis=0)
            test = bin_df[bin_df["region"] == held_out]
            scores = np.full(len(test), np.nan)
            for ti, t_row in enumerate(test.itertuples(index=False)):
                traj = t_row.benioff_traj
                if np.std(traj) == 0:
                    continue
                scores[ti] = float(np.corrcoef(traj, template)[0, 1])
            y_test = (test["window_kind"] == "precursor").astype(int).to_numpy()
            mask = np.isfinite(scores)
            if mask.sum() < 5 or (y_test[mask] == 1).sum() < 2:
                per_region_aucs[held_out] = float("nan")
                continue
            auc = roc_auc_score(y_test[mask], scores[mask])
            per_region_aucs[held_out] = float(auc)
            print(f"[{_ts()}]   {bin_label} test={held_out:<11s}  AUC={auc:.3f}", flush=True)

        macro = float(np.nanmean(list(per_region_aucs.values())))
        print(f"[{_ts()}] [exp17] {bin_label} macro AUC = {macro:.3f}", flush=True)
        rows.append({
            "magnitude_bin": bin_label,
            "n_precursor": int(len(bin_pre)),
            "macro_auc": macro,
            "per_region_aucs": per_region_aucs,
        })

    # Plot
    if rows:
        fig, ax = plt.subplots(figsize=(8, 4.5), dpi=120)
        bin_labels = [r["magnitude_bin"] for r in rows]
        macros = [r["macro_auc"] for r in rows]
        ns = [r["n_precursor"] for r in rows]
        x = np.arange(len(bin_labels))
        ax.bar(x, macros, color="#7d8aa6", alpha=0.6, edgecolor="black")
        for i, (m, n) in enumerate(zip(macros, ns)):
            ax.text(i, m + 0.01, f"AUC={m:.3f}\nn={n}",
                    ha="center", fontsize=9)
        # Per-region scatter
        for ri, r in enumerate(rows):
            for rg_idx, region in enumerate(QUAL_REGIONS):
                auc_rg = r["per_region_aucs"].get(region, float("nan"))
                if np.isfinite(auc_rg):
                    ax.scatter(ri + 0.18 * (rg_idx - 1.5) / 4, auc_rg, s=22,
                               alpha=0.85, color=plt.cm.tab10(rg_idx),
                               edgecolor="white", linewidth=0.5,
                               label=region if ri == 0 else None)
        ax.axhline(0.5, color="black", ls="--", lw=0.8)
        ax.axhline(0.704, color="#c0392b", ls=":", lw=1.0,
                   label="exp14 all-M macro = 0.704")
        ax.set_xticks(x); ax.set_xticklabels(bin_labels, fontsize=9)
        ax.set_ylim(0.40, 0.85)
        ax.set_ylabel("Benioff TLS macro AUC")
        ax.set_title("exp17: Benioff TLS macro AUC by target magnitude")
        ax.legend(fontsize=8, loc="lower right")
        ax.grid(alpha=0.3, axis="y")
        fig.tight_layout()
        fig.savefig(EXP_DIR / "per_mag_plot.png")
        plt.close(fig)
        print(f"[{_ts()}] [plot] per_mag_plot.png", flush=True)

    pd.DataFrame(rows).to_csv(EXP_DIR / "per_mag_table.csv", index=False)
    summary = {
        "experiment": "exp17_per_magnitude_breakdown",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "magnitude_bins": MAG_BINS,
        "rows": rows,
        "interpretation": (
            "If macro AUC scales with target M, larger events have stronger "
            "foreshock signals (consistent with Helmstetter+ 2003 ETAS theory). "
            "If invariant, the signal is M-independent over our range."
        ),
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[{_ts()}] [persist] summary.json", flush=True)
    print(f"[{_ts()}] [exp17] done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
