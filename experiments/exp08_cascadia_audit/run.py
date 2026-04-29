"""exp08_cascadia_audit: forensic audit of the b-vs-null_A z=+3.41 from exp07.

Reads `experiments/exp07_macro_pra2/feature_summary.csv` and per-region
per-window-kind: counts windows, NaN rates on b, n_above_mc distributions.
Specifically tests the hypothesis that Cascadia's AUC=1.00 reflects null A
b-values being mostly NaN (so the AUC is computed on a tiny degenerate set).

If the audit confirms artifact: the b-vs-null_A z=+3.41 is null-of-the-null
contamination, not signal. PRA-2 result is fully null. Methods paper is the
next move.

If the audit shows real high-N b distributions in null A AND consistent
separation from precursor: signal is real, drill into mechanism.

Outputs:
    audit.png        — per-region b-distribution by window kind
    audit_table.csv  — counts and rates
    summary.json     — verdict
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

EXP_DIR = Path(__file__).resolve().parent
SOURCE_CSV = ROOT / "experiments" / "exp07_macro_pra2" / "feature_summary.csv"


def _ts() -> str:
    return dt.datetime.now().strftime("%H:%M:%S")


def per_feature_auc(p: np.ndarray, q: np.ndarray) -> float:
    p = p[np.isfinite(p)]
    q = q[np.isfinite(q)]
    if len(p) == 0 or len(q) == 0:
        return float("nan")
    n_p, n_q = len(p), len(q)
    arr = np.concatenate([p, q])
    order = np.argsort(arr)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, n_p + n_q + 1)
    return float((ranks[:n_p].sum() - n_p * (n_p + 1) / 2) / (n_p * n_q))


def main() -> int:
    print(f"[{_ts()}] [exp08] forensic audit of exp07 b-vs-null_A z=+3.41", flush=True)
    if not SOURCE_CSV.is_file():
        print(f"[{_ts()}] missing {SOURCE_CSV}", flush=True)
        return 2

    df = pd.read_csv(SOURCE_CSV)
    print(f"[{_ts()}] loaded {len(df)} rows from {SOURCE_CSV.name}", flush=True)
    print(f"[{_ts()}] columns: {list(df.columns)}", flush=True)
    print(f"[{_ts()}] regions present: {sorted(df['region'].unique())}", flush=True)
    print(f"[{_ts()}] window_kinds: {sorted(df['window_kind'].unique())}", flush=True)
    print()

    # === Per-region per-window-kind diagnostics ===
    audit_rows = []
    for region in sorted(df["region"].unique()):
        for wk in ("precursor", "null_A", "null_B"):
            sub = df[(df["region"] == region) & (df["window_kind"] == wk)]
            n = len(sub)
            n_b_finite = int(np.isfinite(sub["b"]).sum())
            n_above_mc_med = float(sub["n_above_mc"].median()) if n else float("nan")
            n_above_mc_p25 = float(sub["n_above_mc"].quantile(0.25)) if n else float("nan")
            n_above_mc_p75 = float(sub["n_above_mc"].quantile(0.75)) if n else float("nan")
            b_med = float(sub["b"].median()) if n_b_finite else float("nan")
            b_iqr_lo = float(sub["b"].quantile(0.25)) if n_b_finite else float("nan")
            b_iqr_hi = float(sub["b"].quantile(0.75)) if n_b_finite else float("nan")
            audit_rows.append({
                "region": region, "window_kind": wk,
                "n": n, "n_b_finite": n_b_finite,
                "frac_b_finite": n_b_finite / n if n else 0.0,
                "n_above_mc_p25": n_above_mc_p25,
                "n_above_mc_med": n_above_mc_med,
                "n_above_mc_p75": n_above_mc_p75,
                "b_med": b_med, "b_iqr_lo": b_iqr_lo, "b_iqr_hi": b_iqr_hi,
            })
    audit_df = pd.DataFrame(audit_rows)
    audit_df.to_csv(EXP_DIR / "audit_table.csv", index=False)

    # Print the table grouped by region
    print(f"{'region':<11s} {'kind':<10s} {'N':>4s} {'b_OK':>6s} {'%OK':>5s}  "
          f"{'n_above_mc P25/med/P75':>22s}  {'b med [IQR]':>20s}", flush=True)
    print("-" * 100, flush=True)
    for _, r in audit_df.iterrows():
        print(f"{r['region']:<11s} {r['window_kind']:<10s} {int(r['n']):>4d} "
              f"{int(r['n_b_finite']):>6d} {100*r['frac_b_finite']:>4.0f}%  "
              f"{r['n_above_mc_p25']:>5.0f}/{r['n_above_mc_med']:>4.0f}/{r['n_above_mc_p75']:>4.0f}  "
              f"{r['b_med']:>5.2f} [{r['b_iqr_lo']:.2f}, {r['b_iqr_hi']:.2f}]", flush=True)

    # === Cascadia-specific deep dive ===
    print(f"\n[{_ts()}] === CASCADIA DEEP DIVE ===", flush=True)
    cs = df[df["region"] == "Cascadia"].copy()
    cs_pre = cs[cs["window_kind"] == "precursor"]
    cs_a = cs[cs["window_kind"] == "null_A"]
    cs_b = cs[cs["window_kind"] == "null_B"]
    print(f"[{_ts()}] Cascadia: precursor={len(cs_pre)}, null_A={len(cs_a)}, "
          f"null_B={len(cs_b)}", flush=True)
    print(f"[{_ts()}] Cascadia: b finite — precursor={int(np.isfinite(cs_pre['b']).sum())}/"
          f"{len(cs_pre)}, null_A={int(np.isfinite(cs_a['b']).sum())}/{len(cs_a)}, "
          f"null_B={int(np.isfinite(cs_b['b']).sum())}/{len(cs_b)}", flush=True)

    # Recompute the AUC explicitly to verify
    pre_b = cs_pre["b"].to_numpy(dtype=float)
    a_b = cs_a["b"].to_numpy(dtype=float)
    auc = per_feature_auc(pre_b, a_b)
    pre_b_finite = pre_b[np.isfinite(pre_b)]
    a_b_finite = a_b[np.isfinite(a_b)]
    print(f"[{_ts()}] Cascadia b vs null_A: AUC = {auc:.4f} computed on "
          f"{len(pre_b_finite)} precursor × {len(a_b_finite)} null_A finite values", flush=True)
    if len(pre_b_finite) and len(a_b_finite):
        print(f"[{_ts()}] Cascadia precursor b range: [{pre_b_finite.min():.3f}, "
              f"{pre_b_finite.max():.3f}], median {np.median(pre_b_finite):.3f}", flush=True)
        print(f"[{_ts()}] Cascadia null_A    b range: [{a_b_finite.min():.3f}, "
              f"{a_b_finite.max():.3f}], median {np.median(a_b_finite):.3f}", flush=True)
        sep = pre_b_finite.min() - a_b_finite.max()
        if sep > 0:
            print(f"[{_ts()}] *** Cascadia: precursor MIN ({pre_b_finite.min():.3f}) > "
                  f"null_A MAX ({a_b_finite.max():.3f}); separation = {sep:.3f}", flush=True)
            print(f"[{_ts()}] AUC=1.00 is mathematically correct given the data, but check "
                  f"whether the small effective N is the cause", flush=True)
        else:
            overlap_n = int(((pre_b_finite[:, None] <= a_b_finite[None, :]).sum()))
            print(f"[{_ts()}] precursor and null_A b distributions overlap; "
                  f"AUC < 1.00 expected", flush=True)

    # === Plot: per-region b distributions by window kind ===
    regions_with_data = sorted(df["region"].unique())
    fig, axes = plt.subplots(1, len(regions_with_data),
                             figsize=(3.0 * len(regions_with_data), 4.2),
                             dpi=120, sharey=True)
    if len(regions_with_data) == 1:
        axes = [axes]
    for ax, region in zip(axes, regions_with_data):
        sub = df[df["region"] == region]
        for wk, color in [("precursor", "#c0392b"), ("null_A", "#7d8aa6"),
                          ("null_B", "#999000")]:
            vals = sub.loc[sub["window_kind"] == wk, "b"].to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]
            if len(vals) == 0:
                ax.axhline(0.0, color=color, lw=0)  # placeholder
                continue
            ax.hist(vals, bins=15, alpha=0.55, color=color,
                    label=f"{wk} (n={len(vals)})", density=True,
                    edgecolor="white", linewidth=0.4)
        ax.set_title(region, fontsize=10)
        ax.set_xlabel("b-value")
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("density")
    fig.suptitle("exp08 audit: per-region b-value distributions by window kind",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(EXP_DIR / "audit.png")
    plt.close(fig)
    print(f"\n[{_ts()}] [plot] audit.png", flush=True)

    # === Verdict ===
    cs_a_finite = int(np.isfinite(cs_a["b"]).sum())
    cs_pre_finite = int(np.isfinite(cs_pre["b"]).sum())
    cs_a_nan_rate = 1 - (cs_a_finite / len(cs_a)) if len(cs_a) else 0

    verdict = {}
    if cs_a_finite < 30:
        verdict["headline"] = "ARTIFACT — null A had < 30 finite b values for Cascadia"
        verdict["pra2_b_signal"] = "spurious"
    elif cs_a_finite < 100 and auc > 0.95:
        verdict["headline"] = ("LIKELY ARTIFACT — Cascadia null A had only "
                               f"{cs_a_finite} finite b values (out of {len(cs_a)}); "
                               f"degenerate AUC={auc:.2f}")
        verdict["pra2_b_signal"] = "likely-spurious"
    elif auc > 0.95 and cs_a_finite > 100:
        verdict["headline"] = ("SIGNAL — Cascadia AUC=1.00 holds with "
                               f"{cs_a_finite} finite null A b-values")
        verdict["pra2_b_signal"] = "real-needs-followup"
    else:
        verdict["headline"] = "INTERMEDIATE — needs visual inspection"
        verdict["pra2_b_signal"] = "ambiguous"

    print(f"\n[{_ts()}] === VERDICT ===", flush=True)
    print(f"[{_ts()}] {verdict['headline']}", flush=True)
    print(f"[{_ts()}] cascadia null_A b finite: {cs_a_finite}/{len(cs_a)} "
          f"({100*(1-cs_a_nan_rate):.1f}%)", flush=True)
    print(f"[{_ts()}] cascadia precursor b finite: {cs_pre_finite}/{len(cs_pre)}", flush=True)

    summary = {
        "experiment": "exp08_cascadia_audit",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "source": str(SOURCE_CSV.relative_to(ROOT)),
        "verdict": verdict,
        "cascadia": {
            "n_precursor": int(len(cs_pre)),
            "n_precursor_b_finite": cs_pre_finite,
            "n_null_a": int(len(cs_a)),
            "n_null_a_b_finite": cs_a_finite,
            "null_a_b_nan_rate": float(cs_a_nan_rate),
            "auc_b_vs_null_A_recomputed": float(auc),
        },
        "per_region_table": audit_df.to_dict(orient="records"),
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[{_ts()}] [persist] summary.json", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
