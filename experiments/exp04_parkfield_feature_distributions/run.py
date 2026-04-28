"""exp04_parkfield_feature_distributions: catalog-feature sanity check.

For all Parkfield M>=4.5 target events 2000-2024, compute catalog-derived
features in their 30-day pre-event windows. Compare to features in
aftershock-free random null windows of the same length.

Three features per window (catalog-only; waveform features deferred):
    1. b-value (Aki MLE at fixed Mc=1.0; Bakun-comparable Mc=1.5 is too sparse for 30 days)
    2. b-value drift (slope of b across 3 equal-time sub-windows)
    3. Benioff strain (total, rate, curvature)
    4. n_events above Mc (background diagnostic)

Sanity checks:
    - Per-feature distributions plot for precursor vs null
    - Quick per-feature ROC-AUC (Round C preview only — these AUC numbers are
      NOT pre-registered and are not the headline result of the project)
    - Within-region n=8 precursor windows is far too small for inference;
      this is a mechanics check, not a science claim

Outputs:
    feature_distributions.png       — 2x2 grid of histograms
    feature_summary.csv             — per-window features
    summary.json                    — preview AUCs and gates
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

from src.features.bvalue import bvalue_drift, bvalue_with_bootstrap  # noqa: E402
from src.features.benioff import benioff_features  # noqa: E402

# Pre-cached catalog from exp02
CATALOG_CSV = ROOT / "experiments" / "exp02_parkfield_declustered" / "catalog.csv"

EXP_DIR = Path(__file__).resolve().parent

WINDOW_DAYS = 30
TARGET_M_MIN = 4.5
MC_FOR_FEATURES = 1.0       # lower than calibration Mc — needed for n>=10 events per sub-window
N_NULL_WINDOWS = 200
RANDOM_SEED = 42
NULL_BUFFER_DAYS_BEFORE = 30      # don't sample within 30 days before any M>=4.5
NULL_BUFFER_DAYS_AFTER = 60       # don't sample within 60 days after any M>=4.5 (aftershock-free)
N_DRIFT_SUBWINDOWS = 3
MIN_EVENTS_PER_SUBWINDOW = 10


def _ts() -> str:
    return dt.datetime.now().strftime("%H:%M:%S")


def load_catalog() -> pd.DataFrame:
    if not CATALOG_CSV.is_file():
        raise FileNotFoundError(f"missing catalog cache at {CATALOG_CSV} — run exp02 first")
    # ComCat times are mixed-format (some with microseconds, some without).
    # Use format="ISO8601" to handle both cleanly.
    df = pd.read_csv(CATALOG_CSV)
    df["time"] = pd.to_datetime(df["time"], utc=True, format="ISO8601")
    df = df.sort_values("time").reset_index(drop=True)
    return df


def find_targets(df: pd.DataFrame) -> pd.DataFrame:
    targets = df[df["magnitude"] >= TARGET_M_MIN].copy().sort_values("time").reset_index(drop=True)
    return targets


def sample_null_windows(df: pd.DataFrame, targets: pd.DataFrame, n: int,
                        seed: int = RANDOM_SEED) -> list[tuple]:
    """Sample n random aftershock-free 30-day windows.

    Excludes windows that overlap any (target - NULL_BUFFER_BEFORE_DAYS, target +
    NULL_BUFFER_AFTER_DAYS) interval.
    """
    rng = np.random.default_rng(seed)
    t_first = df["time"].iloc[0]
    t_last = df["time"].iloc[-1] - pd.Timedelta(days=WINDOW_DAYS)

    # Build forbidden intervals (Pandas Timedelta arithmetic)
    forbidden_starts = targets["time"] - pd.Timedelta(days=WINDOW_DAYS + NULL_BUFFER_DAYS_BEFORE)
    forbidden_ends = targets["time"] + pd.Timedelta(days=NULL_BUFFER_DAYS_AFTER)

    windows: list[tuple] = []
    attempts = 0
    span_seconds = (t_last - t_first).total_seconds()
    while len(windows) < n and attempts < n * 30:
        attempts += 1
        offset = rng.uniform(0, span_seconds)
        t_start = t_first + pd.Timedelta(seconds=offset)
        t_end = t_start + pd.Timedelta(days=WINDOW_DAYS)
        # Reject if (t_start, t_end) overlaps any forbidden interval
        overlap = ((t_start <= forbidden_ends) & (t_end >= forbidden_starts)).any()
        if not overlap:
            windows.append((t_start, t_end))
    return windows


def features_for_window(df: pd.DataFrame, t_start, t_end, mc: float) -> dict:
    """Compute all catalog-derived features for events in [t_start, t_end)."""
    mask = (df["time"] >= t_start) & (df["time"] < t_end)
    sub = df.loc[mask]
    n = int(len(sub))

    # Time-since-window-start in seconds for drift computation
    if n > 0:
        ts_seconds = (sub["time"] - t_start).dt.total_seconds().to_numpy()
        mags = sub["magnitude"].to_numpy()
    else:
        ts_seconds = np.array([])
        mags = np.array([])

    out: dict = {
        "t_start": t_start,
        "t_end": t_end,
        "n_events": n,
        "n_above_mc": int((mags >= mc).sum()) if n else 0,
    }

    # b-value (no bootstrap for speed; CIs would just confuse the sanity plot)
    if (mags >= mc).sum() >= 30:
        res = bvalue_with_bootstrap(mags, mc=mc, n_boot=0, mc_correction=0.0)
        out["b"] = float(res.b)
    else:
        out["b"] = float("nan")

    # b-value drift
    drift = bvalue_drift(
        ts_seconds, mags, mc=mc,
        n_subwindows=N_DRIFT_SUBWINDOWS,
        min_events_per_subwindow=MIN_EVENTS_PER_SUBWINDOW,
    )
    out["b_drift_per_day"] = float(drift["drift_slope_per_day"])
    out["b_drift_ok"] = bool(drift["ok"])

    # Benioff strain
    bf = benioff_features(ts_seconds, mags, window_seconds=WINDOW_DAYS * 86400.0)
    out["benioff_total_log10"] = float(np.log10(bf["benioff_total"])) if bf["benioff_total"] > 0 else float("nan")
    out["benioff_curv"] = float(bf["benioff_curv"])

    return out


def per_feature_auc(precursor_vals: np.ndarray, null_vals: np.ndarray) -> float:
    """Mann-Whitney U based AUC. Robust to NaN."""
    p = precursor_vals[np.isfinite(precursor_vals)]
    q = null_vals[np.isfinite(null_vals)]
    if len(p) == 0 or len(q) == 0:
        return float("nan")
    # AUC = P(X_precursor > X_null) + 0.5 * P(X_precursor == X_null)
    n_p, n_q = len(p), len(q)
    rank_arr = np.concatenate([p, q])
    order = np.argsort(rank_arr)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, n_p + n_q + 1)
    R_p = ranks[:n_p].sum()
    U = R_p - n_p * (n_p + 1) / 2
    return float(U / (n_p * n_q))


def main() -> int:
    print(f"[{_ts()}] [exp04] loading catalog from {CATALOG_CSV.relative_to(ROOT)}", flush=True)
    df = load_catalog()
    print(f"[{_ts()}] [exp04] {len(df)} events, {df['time'].iloc[0]} .. {df['time'].iloc[-1]}",
          flush=True)

    targets = find_targets(df)
    print(f"[{_ts()}] [exp04] M>={TARGET_M_MIN} targets: {len(targets)}", flush=True)
    for _, row in targets.iterrows():
        print(f"[{_ts()}]   target  {row['time']}  M={row['magnitude']:.2f}  "
              f"({row['latitude']:.3f},{row['longitude']:.3f})  depth={row['depth_km']:.1f}km",
              flush=True)

    if len(targets) < 2:
        print(f"[{_ts()}] [exp04] too few targets for distribution comparison; abort", flush=True)
        return 2

    # Precursor windows: 30 days BEFORE each target
    print(f"[{_ts()}] [exp04] computing features for {len(targets)} precursor windows...", flush=True)
    pre_rows = []
    for _, row in targets.iterrows():
        t_end = row["time"]
        t_start = t_end - pd.Timedelta(days=WINDOW_DAYS)
        feat = features_for_window(df, t_start, t_end, mc=MC_FOR_FEATURES)
        feat["target_M"] = float(row["magnitude"])
        feat["target_time"] = t_end
        feat["window_kind"] = "precursor"
        pre_rows.append(feat)

    # Null windows
    print(f"[{_ts()}] [exp04] sampling {N_NULL_WINDOWS} aftershock-free null windows...", flush=True)
    null_intervals = sample_null_windows(df, targets, N_NULL_WINDOWS)
    print(f"[{_ts()}] [exp04] sampled {len(null_intervals)} null intervals", flush=True)
    null_rows = []
    for t_start, t_end in null_intervals:
        feat = features_for_window(df, t_start, t_end, mc=MC_FOR_FEATURES)
        feat["target_M"] = float("nan")
        feat["target_time"] = pd.NaT
        feat["window_kind"] = "null"
        null_rows.append(feat)

    rows = pre_rows + null_rows
    feat_df = pd.DataFrame(rows)
    feat_df.to_csv(EXP_DIR / "feature_summary.csv", index=False)
    print(f"[{_ts()}] [exp04] wrote feature_summary.csv ({len(feat_df)} rows)", flush=True)

    # Per-feature AUCs (preview only — NOT pre-registered)
    feat_cols = ["b", "b_drift_per_day", "benioff_total_log10", "benioff_curv", "n_above_mc"]
    aucs = {}
    for col in feat_cols:
        p = feat_df.loc[feat_df["window_kind"] == "precursor", col].to_numpy(dtype=float)
        q = feat_df.loc[feat_df["window_kind"] == "null", col].to_numpy(dtype=float)
        auc = per_feature_auc(p, q)
        aucs[col] = auc
        n_p_fin = int(np.isfinite(p).sum())
        n_q_fin = int(np.isfinite(q).sum())
        print(f"[{_ts()}] [auc] {col:<28s}  AUC={auc:.3f}  "
              f"(n_precursor={n_p_fin}/{len(p)}, n_null={n_q_fin}/{len(q)})", flush=True)

    # Plot
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), dpi=120)
    plot_specs = [
        ("b", "b-value (Mc=1.0)", axes[0, 0]),
        ("b_drift_per_day", "b-value drift (per day)", axes[0, 1]),
        ("benioff_total_log10", "log10 Benioff strain (J^0.5)", axes[0, 2]),
        ("benioff_curv", "Benioff cumulative curvature", axes[1, 0]),
        ("n_above_mc", "n events >= Mc", axes[1, 1]),
    ]
    axes[1, 2].axis("off")  # leave one panel for the AUC table
    for col, label, ax in plot_specs:
        p = feat_df.loc[feat_df["window_kind"] == "precursor", col].to_numpy(dtype=float)
        q = feat_df.loc[feat_df["window_kind"] == "null", col].to_numpy(dtype=float)
        p_fin = p[np.isfinite(p)]
        q_fin = q[np.isfinite(q)]
        all_vals = np.concatenate([p_fin, q_fin])
        if len(all_vals) > 0:
            bins = np.linspace(np.min(all_vals), np.max(all_vals), 20)
            ax.hist(q_fin, bins=bins, alpha=0.55, color="#7d8aa6", density=True,
                    label=f"null (n={len(q_fin)})", edgecolor="white", linewidth=0.4)
            for v in p_fin:
                ax.axvline(v, color="#c0392b", lw=1.5, alpha=0.8)
        ax.set_title(f"{label}\nAUC = {aucs[col]:.3f}")
        ax.set_xlabel(label)
        ax.set_ylabel("density")
        ax.legend(fontsize=8, loc="upper right",
                  handles=[plt.Rectangle((0, 0), 1, 1, color="#7d8aa6", alpha=0.55, label=f"null (n={len(q_fin)})"),
                           plt.Line2D([0], [0], color="#c0392b", lw=1.5,
                                      label=f"precursor (n={len(p_fin)})")])
        ax.grid(alpha=0.3)

    # AUC table panel
    ax_t = axes[1, 2]
    ax_t.axis("off")
    table_text = "preview AUC (NOT pre-registered)\n\n"
    for col in feat_cols:
        table_text += f"  {col:<24s}  {aucs[col]:.3f}\n"
    table_text += (f"\nn_precursor windows: {(feat_df['window_kind']=='precursor').sum()}\n"
                   f"n_null windows:      {(feat_df['window_kind']=='null').sum()}\n"
                   f"Mc used:             {MC_FOR_FEATURES}\n"
                   f"window:              {WINDOW_DAYS} days")
    ax_t.text(0.0, 0.95, table_text, family="monospace", fontsize=9, va="top")

    fig.suptitle(f"exp04: Parkfield 30-day catalog feature distributions  "
                 f"(precursor red, null gray)", y=1.0, fontsize=11)
    fig.tight_layout()
    fig.savefig(EXP_DIR / "feature_distributions.png")
    plt.close(fig)
    print(f"[{_ts()}] [plot] wrote feature_distributions.png", flush=True)

    # === Sanity gates ===
    feat_computable = {
        col: int(np.isfinite(feat_df[col].to_numpy(dtype=float)).sum())
        for col in feat_cols
    }
    n_precursor = int((feat_df["window_kind"] == "precursor").sum())
    n_null = int((feat_df["window_kind"] == "null").sum())
    n_drift_ok_precursor = int((feat_df.loc[feat_df["window_kind"] == "precursor", "b_drift_ok"]).sum())

    sanity_pass = (
        feat_computable["b"] >= n_precursor                       # b is computable for ALL precursor windows
        and feat_computable["benioff_total_log10"] >= n_precursor
        and n_null > 0
    )
    print(f"\n[{_ts()}] [gate] feature computability:", flush=True)
    for col, n_ok in feat_computable.items():
        print(f"[{_ts()}]   {col:<24s}  {n_ok:>4d} / {n_precursor + n_null} computable", flush=True)
    print(f"[{_ts()}] [gate] precursor b_drift OK:  {n_drift_ok_precursor} / {n_precursor}", flush=True)
    print(f"[{_ts()}] [gate] >>> exp04 SANITY {'PASS' if sanity_pass else 'FAIL'} <<<\n", flush=True)

    summary = {
        "experiment": "exp04_parkfield_feature_distributions",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "params": {
            "window_days": WINDOW_DAYS,
            "target_M_min": TARGET_M_MIN,
            "Mc_for_features": MC_FOR_FEATURES,
            "n_null_target": N_NULL_WINDOWS,
            "n_drift_subwindows": N_DRIFT_SUBWINDOWS,
            "min_events_per_subwindow": MIN_EVENTS_PER_SUBWINDOW,
            "null_buffer_days_before": NULL_BUFFER_DAYS_BEFORE,
            "null_buffer_days_after": NULL_BUFFER_DAYS_AFTER,
        },
        "counts": {
            "n_precursor": n_precursor,
            "n_null": n_null,
            "n_drift_ok_precursor": n_drift_ok_precursor,
            "feature_computability": feat_computable,
        },
        "preview_aucs_NOT_PREREGISTERED": aucs,
        "gates": {
            "sanity_pass": bool(sanity_pass),
        },
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[{_ts()}] [persist] wrote summary.json", flush=True)
    return 0 if sanity_pass else 2


if __name__ == "__main__":
    sys.exit(main())
