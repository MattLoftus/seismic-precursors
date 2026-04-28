"""exp02_parkfield_declustered: Reasenberg-style declustering via ZBZ 2013, then
b-value sensitivity to Mc and to the declustering decision.

Hypothesis (from exp01 diagnosis): the 9.67% rel-error gap between our
b=0.813 and Bakun's 0.9 is driven by the 2004 M6.0 mainshock-aftershock
sequence in our 2000-2024 window. If we (a) decluster, then (b) restrict to
the pre-2004 background period, our b should converge toward Bakun's 0.9.

Outputs:
- catalog.csv                — cached full event records (15k+ rows)
- bimodality_log_eta.png     — ZBZ proximity-statistic histogram, w/ threshold
- fmd_full_vs_decl.png       — overlaid FMDs for full / declustered / pre-2004
- bvalue_vs_mc.png           — b vs Mc curve for each catalog flavor
- summary.json               — machine-readable result
- run.log                    — console log
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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.data import fetch_comcat_catalog  # noqa: E402
from src.features.bvalue import (  # noqa: E402
    aki_bvalue,
    bvalue_with_bootstrap,
    magnitude_of_completeness_max_curvature,
)
from src.features.declustering import zbz_decluster  # noqa: E402

PARKFIELD_LAT = 35.85
PARKFIELD_LON = -120.40
PARKFIELD_RADIUS_KM = 50.0
M_MIN = 1.0
START = dt.datetime(2000, 1, 1)
END = dt.datetime(2025, 1, 1)
PRE_2004_END = dt.datetime(2004, 9, 28)  # 2004-09-28 = Parkfield M6.0 mainshock

BAKUN_2005_B = 0.9
B_TOLERANCE = 0.10
ZBZ_B_INPUT = 0.9   # use Bakun reference b for the proximity calc (decoupled from output b)

MC_GRID = [1.00, 1.15, 1.25, 1.35, 1.50]

EXP_DIR = Path(__file__).resolve().parent
DATA_CACHE = EXP_DIR / "catalog.csv"


def _ts() -> str:
    return dt.datetime.now().strftime("%H:%M:%S")


def plot_bimodality(decl_result, out_path: Path) -> None:
    log_eta = decl_result.log_eta_nn
    finite = log_eta[np.isfinite(log_eta)]

    fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=120)
    bins = np.linspace(np.percentile(finite, 0.5), np.percentile(finite, 99.5), 80)
    ax.hist(finite, bins=bins, color="#7d8aa6", edgecolor="white", linewidth=0.4)
    ax.axvline(decl_result.eta_threshold, color="#c0392b", lw=1.6,
               label=f"threshold log10(eta_0) = {decl_result.eta_threshold:.2f}")
    ax.set_xlabel("log10(eta_NN)   [proximity statistic — small = clustered, large = background]")
    ax.set_ylabel("# events")
    ax.set_title(
        f"ZBZ 2013 nearest-neighbor proximity histogram\n"
        f"N={decl_result.n_total}, "
        f"{decl_result.n_background} background ({100 * decl_result.n_background / decl_result.n_total:.1f}%), "
        f"{decl_result.n_clustered} clustered "
        f"(b_in={decl_result.b_used}, d_f={decl_result.d_f})"
    )
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_fmd_overlay(catalogs: dict[str, np.ndarray], mc: float, out_path: Path) -> None:
    dm = 0.1
    fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=120)
    colors = {"full": "#7d8aa6", "declustered": "#1f5fa6", "decl_pre2004": "#c0392b"}
    markers = {"full": "o", "declustered": "s", "decl_pre2004": "^"}
    for label, mags in catalogs.items():
        if len(mags) == 0:
            continue
        edges = np.arange(np.floor(mags.min() * 10) / 10,
                          np.ceil(mags.max() * 10) / 10 + dm, dm)
        counts, _ = np.histogram(mags, bins=edges)
        centers = edges[:-1] + dm / 2
        cum = np.cumsum(counts[::-1])[::-1]
        ax.semilogy(centers, np.where(cum > 0, cum, np.nan),
                    marker=markers.get(label, "o"), ms=3.8, ls="none",
                    color=colors.get(label, None),
                    label=f"{label} (N={len(mags)})")
    ax.axvline(mc, color="k", ls="--", lw=0.8, alpha=0.6, label=f"Mc = {mc:.2f}")
    ax.set_xlabel("Magnitude")
    ax.set_ylabel("cumulative N(M >= m)")
    ax.set_title(f"Parkfield FMD: full vs declustered vs declustered+pre-2004 (Mc={mc:.2f})")
    ax.legend(loc="lower left", fontsize=9, framealpha=0.95)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_b_vs_mc(table: dict[str, dict[float, dict]], out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=120)
    colors = {"full": "#7d8aa6", "declustered": "#1f5fa6", "decl_pre2004": "#c0392b"}
    for label, by_mc in table.items():
        mcs = sorted(by_mc.keys())
        bs = [by_mc[m]["b"] for m in mcs]
        b_lo = [by_mc[m]["ci_lo"] for m in mcs]
        b_hi = [by_mc[m]["ci_hi"] for m in mcs]
        ax.errorbar(mcs, bs,
                    yerr=[np.array(bs) - np.array(b_lo), np.array(b_hi) - np.array(bs)],
                    marker="o", ms=5, ls="-", lw=1.2, capsize=4,
                    color=colors.get(label, None), label=label)
    ax.axhline(BAKUN_2005_B, color="black", ls="--", lw=1.0,
               label=f"Bakun+ 2005 reference b = {BAKUN_2005_B}")
    ax.fill_between([min(MC_GRID) - 0.05, max(MC_GRID) + 0.05],
                    BAKUN_2005_B * (1 - B_TOLERANCE),
                    BAKUN_2005_B * (1 + B_TOLERANCE),
                    alpha=0.10, color="black", label=f"+/-{int(B_TOLERANCE * 100)}% gate")
    ax.set_xlim(min(MC_GRID) - 0.05, max(MC_GRID) + 0.05)
    ax.set_xlabel("Mc")
    ax.set_ylabel("b-value (Aki MLE, error bars = bootstrap 95% CI)")
    ax.set_title("b-value sensitivity to Mc and declustering / pre-2004 restriction")
    ax.legend(fontsize=8, framealpha=0.95)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def main() -> int:
    print(f"[{_ts()}] [exp02] start", flush=True)
    print(f"[{_ts()}] [exp02] cache: {DATA_CACHE}", flush=True)

    df = fetch_comcat_catalog(
        lat=PARKFIELD_LAT,
        lon=PARKFIELD_LON,
        radius_km=PARKFIELD_RADIUS_KM,
        start=START,
        end=END,
        m_min=M_MIN,
        chunk_years=2,
        cache_path=DATA_CACHE,
        log=lambda m: print(f"[{_ts()}] {m}", flush=True),
    )
    print(f"[{_ts()}] [exp02] catalog has {len(df)} events, "
          f"M [{df['magnitude'].min():.2f}, {df['magnitude'].max():.2f}]", flush=True)

    # === Decluster ===
    decl = zbz_decluster(
        df,
        b_value=ZBZ_B_INPUT,
        log=lambda m: print(f"[{_ts()}] {m}", flush=True),
    )
    plot_bimodality(decl, EXP_DIR / "bimodality_log_eta.png")
    print(f"[{_ts()}] [plot] wrote bimodality_log_eta.png", flush=True)

    bg_mask = decl.is_background
    decl_mags = df.loc[bg_mask, "magnitude"].to_numpy()
    decl_times = df.loc[bg_mask, "time"]

    pre2004_full = df["time"] < pd_utc(PRE_2004_END)
    pre2004_decl_mask = bg_mask & pre2004_full.to_numpy()
    pre2004_decl_mags = df.loc[pre2004_decl_mask, "magnitude"].to_numpy()

    full_mags = df["magnitude"].to_numpy()

    print(f"[{_ts()}] [exp02] catalog sizes: "
          f"full={len(full_mags)}, declustered={len(decl_mags)}, "
          f"declustered+pre-2004={len(pre2004_decl_mags)}", flush=True)

    # === b-value sweep over Mc x catalog flavor ===
    catalogs = {
        "full": full_mags,
        "declustered": decl_mags,
        "decl_pre2004": pre2004_decl_mags,
    }
    table: dict[str, dict[float, dict]] = {k: {} for k in catalogs}
    for label, mags in catalogs.items():
        for mc in MC_GRID:
            if (mags >= mc).sum() < 30:
                table[label][mc] = {"b": float("nan"), "ci_lo": float("nan"), "ci_hi": float("nan"),
                                    "n_above_mc": int((mags >= mc).sum())}
                continue
            res = bvalue_with_bootstrap(mags, mc=mc, n_boot=500, mc_correction=0.0)
            table[label][mc] = {
                "b": res.b,
                "b_se": res.b_se,
                "ci_lo": res.b_boot_lo,
                "ci_hi": res.b_boot_hi,
                "n_above_mc": res.n_above_mc,
            }
            print(f"[{_ts()}] [bvalue] {label:<14s} Mc={mc:.2f}  "
                  f"b={res.b:.4f} +/- {res.b_se:.4f}  "
                  f"CI=[{res.b_boot_lo:.4f}, {res.b_boot_hi:.4f}]  "
                  f"n>=Mc={res.n_above_mc}", flush=True)

    plot_b_vs_mc(table, EXP_DIR / "bvalue_vs_mc.png")
    print(f"[{_ts()}] [plot] wrote bvalue_vs_mc.png", flush=True)

    # FMD overlay at the canonical Mc=1.35 from exp01
    plot_fmd_overlay(
        {k: catalogs[k] for k in ["full", "declustered", "decl_pre2004"]},
        mc=1.35,
        out_path=EXP_DIR / "fmd_full_vs_decl.png",
    )
    print(f"[{_ts()}] [plot] wrote fmd_full_vs_decl.png", flush=True)

    # === Calibration gate(s) ===
    headline_mc = 1.35
    gates = {}
    for label in catalogs:
        b = table[label][headline_mc]["b"]
        ci_lo = table[label][headline_mc]["ci_lo"]
        ci_hi = table[label][headline_mc]["ci_hi"]
        if not np.isfinite(b):
            gates[label] = {"b": None, "rel_err": None, "in_ci": None, "passed": None}
            continue
        rel_err = abs(b - BAKUN_2005_B) / BAKUN_2005_B
        in_ci = ci_lo <= BAKUN_2005_B <= ci_hi
        passed = (rel_err <= B_TOLERANCE) or in_ci
        gates[label] = {
            "b": b, "rel_err": rel_err, "in_ci": bool(in_ci), "passed": bool(passed),
        }
        print(f"[{_ts()}] [gate@Mc=1.35] {label:<14s}  "
              f"b={b:.4f}  rel_err={rel_err:.2%}  Bakun_in_CI={in_ci}  "
              f"verdict={'PASS' if passed else 'FAIL'}", flush=True)

    summary = {
        "experiment": "exp02_parkfield_declustered",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "catalog": {
            "n_total": len(full_mags),
            "n_declustered": len(decl_mags),
            "n_decl_pre2004": len(pre2004_decl_mags),
        },
        "declustering": {
            "method": "Zaliapin & Ben-Zion 2013",
            "b_input_to_zbz": ZBZ_B_INPUT,
            "d_f": decl.d_f,
            "p": decl.p,
            "q": decl.q,
            "log_eta_threshold": decl.eta_threshold,
            "n_background": decl.n_background,
            "n_clustered": decl.n_clustered,
        },
        "bvalue_table": {
            label: {f"Mc={mc:.2f}": v for mc, v in by_mc.items()}
            for label, by_mc in table.items()
        },
        "calibration_gates_at_Mc_1p35": gates,
        "headline_decl_pre2004_b_at_Mc_1p35": gates["decl_pre2004"]["b"],
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=lambda o: float(o) if isinstance(o, np.floating) else o)
    print(f"[{_ts()}] [persist] wrote summary.json", flush=True)
    print(f"[{_ts()}] [exp02] done", flush=True)
    return 0


def pd_utc(ts: dt.datetime):
    """Convert a naive datetime to a tz-aware pandas timestamp (UTC)."""
    import pandas as pd
    return pd.Timestamp(ts, tz="UTC")


if __name__ == "__main__":
    sys.exit(main())
