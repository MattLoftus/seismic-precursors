"""exp05_california_features: full pre-reg protocol applied to California (one region preview).

Implements the locked protocol from `papers/pre_registration.md` end-to-end on
the California polygon (32-42 N, 125-114 W). This is a single-region PREVIEW of
Round D, not the headline result.

DEVIATIONS FROM PRE-REG (preview-only, NOT carried into Round D):
1. Catalog M_min = 2.5 instead of 1.5. Reason: California M≥1.5 over 25 years
   is ~1.25M events, exceeding ComCat per-query limits and requiring monthly
   chunking. M≥2.5 gives ~12k events in a single query, ~10-15 min runtime.
   For Round D headline, the protocol uses M≥1.5 with monthly chunking.
2. Mc estimated from M≥2.5 data; we cannot see the rolloff below 2.5, so the
   max-curvature method may be biased. Round D will not have this issue.
3. Null C (colored noise) skipped; this is catalog-only, no waveforms.

Outputs:
    catalog.csv                    — full California catalog (gitignored)
    bvalue_vs_mc.png               — Mc plateau diagnostic
    zbz_bimodality.png             — ZBZ proximity histogram, log10(η_0)=-5.0
    feature_distributions.png      — 6-panel histogram (b, b-drift, Benioff, n)
    auc_summary.png                — AUC + bootstrap + permutation per feature
    summary.json                   — machine-readable
    feature_summary.csv            — per-window features
"""
from __future__ import annotations

import datetime as dt
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.data import fetch_comcat_catalog_bbox  # noqa: E402
from src.features.benioff import benioff_features  # noqa: E402
from src.features.bvalue import (  # noqa: E402
    bvalue_drift,
    bvalue_with_bootstrap,
    magnitude_of_completeness_max_curvature,
    aki_bvalue,
)
from src.features.declustering import zbz_decluster  # noqa: E402
from src.regions import CALIFORNIA  # noqa: E402

# === Pre-reg locked parameters ===
WINDOW_DAYS = 30
TARGET_M_MIN = 4.5
ZBZ_LOG_ETA_THRESHOLD = -5.0      # ZBZ 2013 published default
ZBZ_B_INPUT = 1.0                 # ZBZ 2013 published default
N_DRIFT_SUBWINDOWS = 3
MIN_EVENTS_PER_SUBWINDOW = 10
NULL_BUFFER_DAYS_BEFORE = 30
NULL_BUFFER_DAYS_AFTER = 60
N_NULL_A = 200                    # per pre-reg §5.2
N_NULL_B = 200                    # per pre-reg §5.3
N_BOOT = 1000
N_PERM = 1000
RANDOM_SEED = 42

# === Preview deviations (not Round D) ===
CATALOG_M_MIN = 2.5               # pre-reg says 1.5; preview deviation
START = dt.datetime(2000, 1, 1)
END = dt.datetime(2025, 1, 1)

# Mc grid for the plateau check
MC_GRID = [2.5, 2.7, 2.9, 3.1, 3.3, 3.5]

EXP_DIR = Path(__file__).resolve().parent
CATALOG_CACHE = EXP_DIR / "catalog.csv"


def _ts() -> str:
    return dt.datetime.now().strftime("%H:%M:%S")


def estimate_mc_with_plateau_check(magnitudes: np.ndarray, mc_grid: list[float]) -> tuple[float, dict]:
    """Per pre-reg §1.2: max-curvature + Woessner +0.20, then iterate Mc up by 0.15
    until the b-vs-Mc curve plateaus.

    Returns (Mc_chosen, diagnostic dict).
    """
    mc_raw = magnitude_of_completeness_max_curvature(magnitudes)
    mc_woessner = mc_raw + 0.20
    bs = []
    for mc in mc_grid:
        if (magnitudes >= mc).sum() < 100:
            bs.append(float("nan"))
            continue
        b, _, _ = aki_bvalue(magnitudes, mc=mc)
        bs.append(float(b))
    # Look for the plateau: b-stops-rising. Stop when consecutive b values change
    # by < 0.02 (b is "flat enough"); pick that Mc as the chosen value. If never
    # flat, use the highest Mc in the grid.
    chosen = None
    for i in range(1, len(bs)):
        if not (math.isfinite(bs[i]) and math.isfinite(bs[i - 1])):
            continue
        if abs(bs[i] - bs[i - 1]) < 0.02:
            chosen = mc_grid[i]
            break
    if chosen is None:
        # Fall back to the highest Mc with finite b
        for i in reversed(range(len(bs))):
            if math.isfinite(bs[i]):
                chosen = mc_grid[i]
                break
    if chosen is None:
        chosen = mc_woessner
    # Final clamp: never go below the Woessner value
    chosen = max(chosen, mc_woessner)
    return float(chosen), {
        "mc_raw": float(mc_raw),
        "mc_woessner": float(mc_woessner),
        "mc_grid": list(mc_grid),
        "b_at_mc_grid": bs,
        "chosen_via": "plateau-stop" if any(
            i > 0 and math.isfinite(bs[i]) and math.isfinite(bs[i - 1]) and abs(bs[i] - bs[i - 1]) < 0.02
            for i in range(1, len(bs))
        ) else "highest-finite-Mc",
    }


def features_for_window(df: pd.DataFrame, t_start, t_end, mc: float) -> dict:
    mask = (df["time"] >= t_start) & (df["time"] < t_end)
    sub = df.loc[mask]
    n = int(len(sub))
    if n > 0:
        ts_seconds = (sub["time"] - t_start).dt.total_seconds().to_numpy()
        mags = sub["magnitude"].to_numpy()
    else:
        ts_seconds = np.array([])
        mags = np.array([])

    out = {
        "t_start": t_start,
        "t_end": t_end,
        "n_events": n,
        "n_above_mc": int((mags >= mc).sum()) if n else 0,
    }
    if (mags >= mc).sum() >= 30:
        res = bvalue_with_bootstrap(mags, mc=mc, n_boot=0, mc_correction=0.0)
        out["b"] = float(res.b)
    else:
        out["b"] = float("nan")
    drift = bvalue_drift(
        ts_seconds, mags, mc=mc,
        n_subwindows=N_DRIFT_SUBWINDOWS,
        min_events_per_subwindow=MIN_EVENTS_PER_SUBWINDOW,
    )
    out["b_drift_per_day"] = float(drift["drift_slope_per_day"])
    out["b_drift_ok"] = bool(drift["ok"])
    bf = benioff_features(ts_seconds, mags, window_seconds=WINDOW_DAYS * 86400.0)
    out["benioff_total_log10"] = float(np.log10(bf["benioff_total"])) if bf["benioff_total"] > 0 else float("nan")
    out["benioff_curv"] = float(bf["benioff_curv"])
    return out


def sample_null_windows(df: pd.DataFrame, target_times, n: int, seed: int,
                        kind: str) -> list[tuple]:
    """Sample n random null windows of WINDOW_DAYS length.

    kind="A": aftershock-free; window must be ≥30d before AND ≥60d after every target.
    kind="B": minimal exclusion; window must not overlap a 7-day post-event of any target.
    """
    rng = np.random.default_rng(seed)
    t_first = df["time"].iloc[0]
    t_last = df["time"].iloc[-1] - pd.Timedelta(days=WINDOW_DAYS)
    target_times = pd.to_datetime(target_times, utc=True)
    if kind == "A":
        forbid_starts = target_times - pd.Timedelta(days=WINDOW_DAYS + NULL_BUFFER_DAYS_BEFORE)
        forbid_ends = target_times + pd.Timedelta(days=NULL_BUFFER_DAYS_AFTER)
    elif kind == "B":
        forbid_starts = target_times - pd.Timedelta(days=WINDOW_DAYS)
        forbid_ends = target_times + pd.Timedelta(days=7)
    else:
        raise ValueError(f"unknown kind {kind}")

    span_seconds = (t_last - t_first).total_seconds()
    windows: list[tuple] = []
    attempts = 0
    while len(windows) < n and attempts < n * 30:
        attempts += 1
        offset = rng.uniform(0, span_seconds)
        t_start = t_first + pd.Timedelta(seconds=offset)
        t_end = t_start + pd.Timedelta(days=WINDOW_DAYS)
        overlap = ((t_start <= forbid_ends) & (t_end >= forbid_starts)).any()
        if not overlap:
            windows.append((t_start, t_end))
    return windows


def per_feature_auc(precursor_vals: np.ndarray, null_vals: np.ndarray) -> float:
    p = precursor_vals[np.isfinite(precursor_vals)]
    q = null_vals[np.isfinite(null_vals)]
    if len(p) == 0 or len(q) == 0:
        return float("nan")
    n_p, n_q = len(p), len(q)
    arr = np.concatenate([p, q])
    order = np.argsort(arr)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, n_p + n_q + 1)
    R_p = ranks[:n_p].sum()
    U = R_p - n_p * (n_p + 1) / 2
    return float(U / (n_p * n_q))


def auc_with_significance(precursor_vals: np.ndarray, null_vals: np.ndarray,
                          n_boot: int = N_BOOT, n_perm: int = N_PERM,
                          seed: int = RANDOM_SEED) -> dict:
    """Pre-reg §6.3: bootstrap CI95 + permutation z + Bonferroni-corrected p."""
    rng = np.random.default_rng(seed)
    auc = per_feature_auc(precursor_vals, null_vals)
    p = precursor_vals[np.isfinite(precursor_vals)]
    q = null_vals[np.isfinite(null_vals)]
    if len(p) < 2 or len(q) < 2:
        return {"auc": auc, "ci_lo": float("nan"), "ci_hi": float("nan"),
                "perm_z": float("nan"), "perm_p": float("nan"), "n_p": len(p), "n_q": len(q)}

    # Bootstrap
    boot_aucs = np.empty(n_boot)
    for i in range(n_boot):
        ps = rng.choice(p, size=len(p), replace=True)
        qs = rng.choice(q, size=len(q), replace=True)
        boot_aucs[i] = per_feature_auc(ps, qs)
    ci_lo, ci_hi = np.nanpercentile(boot_aucs, [2.5, 97.5])

    # Permutation
    pooled = np.concatenate([p, q])
    perm_aucs = np.empty(n_perm)
    for i in range(n_perm):
        rng.shuffle(pooled)
        perm_aucs[i] = per_feature_auc(pooled[:len(p)], pooled[len(p):])
    perm_mean = float(np.nanmean(perm_aucs))
    perm_std = float(np.nanstd(perm_aucs))
    z = (auc - perm_mean) / perm_std if perm_std > 0 else float("nan")
    # Two-sided p from |z|
    from scipy.stats import norm
    p_val = 2 * (1 - norm.cdf(abs(z)))
    return {
        "auc": float(auc), "ci_lo": float(ci_lo), "ci_hi": float(ci_hi),
        "perm_mean": perm_mean, "perm_std": perm_std,
        "perm_z": float(z), "perm_p": float(p_val),
        "n_p": int(len(p)), "n_q": int(len(q)),
    }


def main() -> int:
    print(f"[{_ts()}] [exp05] start, region={CALIFORNIA.name}", flush=True)
    print(f"[{_ts()}] [exp05] bbox lat=[{CALIFORNIA.lat_min},{CALIFORNIA.lat_max}] "
          f"lon=[{CALIFORNIA.lon_min},{CALIFORNIA.lon_max}]", flush=True)
    print(f"[{_ts()}] [exp05] catalog M_min={CATALOG_M_MIN} (PREVIEW deviation; pre-reg says 1.5)",
          flush=True)

    # === Catalog ===
    df = fetch_comcat_catalog_bbox(
        lat_min=CALIFORNIA.lat_min, lat_max=CALIFORNIA.lat_max,
        lon_min=CALIFORNIA.lon_min, lon_max=CALIFORNIA.lon_max,
        start=START, end=END, m_min=CATALOG_M_MIN,
        chunk_years=2, cache_path=CATALOG_CACHE,
        log=lambda m: print(f"[{_ts()}] {m}", flush=True),
    )
    print(f"[{_ts()}] [exp05] catalog has {len(df)} events", flush=True)

    # === Mc plateau check ===
    mc, mc_diag = estimate_mc_with_plateau_check(df["magnitude"].to_numpy(), MC_GRID)
    print(f"[{_ts()}] [Mc] raw={mc_diag['mc_raw']:.2f}, woessner={mc_diag['mc_woessner']:.2f}, "
          f"chosen={mc:.2f} via {mc_diag['chosen_via']}", flush=True)
    print(f"[{_ts()}] [Mc] b-vs-Mc grid: {list(zip(mc_diag['mc_grid'], mc_diag['b_at_mc_grid']))}",
          flush=True)

    # === Decluster catalog ===
    print(f"[{_ts()}] [zbz] declustering catalog (N={len(df)}, log10(η_0)={ZBZ_LOG_ETA_THRESHOLD})",
          flush=True)
    decl = zbz_decluster(df, b_value=ZBZ_B_INPUT, eta_threshold=ZBZ_LOG_ETA_THRESHOLD,
                         log=lambda m: print(f"[{_ts()}] {m}", flush=True))
    df_bg = df.loc[decl.is_background].reset_index(drop=True)
    print(f"[{_ts()}] [zbz] background events: {len(df_bg)} / {len(df)}", flush=True)

    # === Targets: M>=4.5 background-only ===
    targets_all = df.loc[df["magnitude"] >= TARGET_M_MIN].reset_index(drop=True)
    targets_bg = df_bg.loc[df_bg["magnitude"] >= TARGET_M_MIN].reset_index(drop=True)
    print(f"[{_ts()}] [targets] M>=4.5 total={len(targets_all)}, "
          f"background-only={len(targets_bg)}", flush=True)

    if len(targets_bg) < 8:
        print(f"[{_ts()}] [exp05] WARNING: only {len(targets_bg)} declustered targets "
              f"(pre-reg minimum: 8)", flush=True)

    # === Precursor windows (with overlap rejection per pre-reg §5.1) ===
    # The forbidden zone per OTHER target t' is [t'-30, t'+60] (90 days, the same
    # zone used for nulls). Reject our precursor window [t-30, t) if it overlaps
    # any other target's forbidden zone.
    # Standard interval overlap: [a,b] overlaps [c,d] iff a<=d AND b>=c.
    print(f"[{_ts()}] [windows] building precursor windows with overlap rejection "
          f"(forbidden zone per target = [-30d, +60d])", flush=True)
    forbid_starts = targets_bg["time"] - pd.Timedelta(days=NULL_BUFFER_DAYS_BEFORE)   # t' - 30
    forbid_ends = targets_bg["time"] + pd.Timedelta(days=NULL_BUFFER_DAYS_AFTER)      # t' + 60
    pre_rows = []
    rejected = 0
    for _, row in targets_bg.iterrows():
        t_end = row["time"]
        t_start = t_end - pd.Timedelta(days=WINDOW_DAYS)
        own_mask = forbid_starts == t_end - pd.Timedelta(days=NULL_BUFFER_DAYS_BEFORE)
        other_starts = forbid_starts[~own_mask]
        other_ends = forbid_ends[~own_mask]
        # Window [t_start, t_end) overlaps [other_starts, other_ends] iff
        # t_start <= other_ends AND t_end >= other_starts.
        overlap = ((t_start <= other_ends) & (t_end >= other_starts)).any()
        if overlap:
            rejected += 1
            continue
        feat = features_for_window(df, t_start, t_end, mc=mc)
        feat["target_M"] = float(row["magnitude"])
        feat["target_time"] = t_end
        feat["window_kind"] = "precursor"
        pre_rows.append(feat)
    print(f"[{_ts()}] [windows] precursor: kept {len(pre_rows)}, rejected {rejected}", flush=True)

    # === Null A and Null B windows ===
    print(f"[{_ts()}] [windows] sampling Null A (aftershock-free)...", flush=True)
    null_a = sample_null_windows(df, targets_bg["time"].tolist(), N_NULL_A, RANDOM_SEED, kind="A")
    null_b = sample_null_windows(df, targets_all["time"].tolist(), N_NULL_B, RANDOM_SEED + 1, kind="B")
    print(f"[{_ts()}] [windows] sampled Null A: {len(null_a)}, Null B: {len(null_b)}", flush=True)

    null_a_rows = []
    for t_start, t_end in null_a:
        feat = features_for_window(df, t_start, t_end, mc=mc)
        feat["target_M"] = float("nan")
        feat["target_time"] = pd.NaT
        feat["window_kind"] = "null_A"
        null_a_rows.append(feat)
    null_b_rows = []
    for t_start, t_end in null_b:
        feat = features_for_window(df, t_start, t_end, mc=mc)
        feat["target_M"] = float("nan")
        feat["target_time"] = pd.NaT
        feat["window_kind"] = "null_B"
        null_b_rows.append(feat)

    feat_df = pd.DataFrame(pre_rows + null_a_rows + null_b_rows)
    feat_df.to_csv(EXP_DIR / "feature_summary.csv", index=False)
    print(f"[{_ts()}] [persist] feature_summary.csv ({len(feat_df)} rows)", flush=True)

    # === AUC + significance per feature, per null type ===
    feat_cols = ["b", "b_drift_per_day", "benioff_total_log10", "benioff_curv", "n_above_mc"]
    auc_table = {}
    for null_kind in ("null_A", "null_B"):
        for col in feat_cols:
            p = feat_df.loc[feat_df["window_kind"] == "precursor", col].to_numpy(dtype=float)
            q = feat_df.loc[feat_df["window_kind"] == null_kind, col].to_numpy(dtype=float)
            res = auc_with_significance(p, q)
            auc_table[f"{col}__{null_kind}"] = res
            sig_3sigma = res["ci_lo"] > 0.5 and abs(res["auc"] - 0.5) > 0.05 and res["perm_z"] > 3
            sig_bonferroni = res["perm_p"] < 0.05 / 36
            print(f"[{_ts()}] [auc] {col:<22s} vs {null_kind:6s}  "
                  f"AUC={res['auc']:.3f} CI=[{res['ci_lo']:.3f},{res['ci_hi']:.3f}] "
                  f"z={res['perm_z']:+.2f} p={res['perm_p']:.2e}  "
                  f"3σ={'Y' if sig_3sigma else 'n'}  Bonf={'Y' if sig_bonferroni else 'n'}",
                  flush=True)

    # === Plot: Mc curve ===
    fig, ax = plt.subplots(figsize=(7, 4), dpi=120)
    ax.plot(mc_diag["mc_grid"], mc_diag["b_at_mc_grid"], "o-", color="#1f5fa6")
    ax.axvline(mc, color="#c0392b", ls="--", label=f"chosen Mc={mc:.2f}")
    ax.axvline(mc_diag["mc_woessner"], color="black", ls=":", label=f"Woessner={mc_diag['mc_woessner']:.2f}")
    ax.set_xlabel("Mc")
    ax.set_ylabel("b-value")
    ax.set_title(f"b-vs-Mc plateau check, California (N={len(df)})")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(EXP_DIR / "bvalue_vs_mc.png")
    plt.close(fig)

    # === Plot: ZBZ bimodality ===
    fig, ax = plt.subplots(figsize=(8, 4), dpi=120)
    finite = decl.log_eta_nn[np.isfinite(decl.log_eta_nn)]
    ax.hist(finite, bins=80, color="#7d8aa6", edgecolor="white", linewidth=0.4)
    ax.axvline(ZBZ_LOG_ETA_THRESHOLD, color="#c0392b", lw=1.5,
               label=f"pre-reg log10(η_0)={ZBZ_LOG_ETA_THRESHOLD}")
    ax.set_xlabel("log10(η_NN)")
    ax.set_ylabel("count")
    ax.set_title(f"ZBZ proximity histogram, California N={len(df)} "
                 f"(bg={int(decl.is_background.sum())}, cl={int((~decl.is_background).sum())})")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(EXP_DIR / "zbz_bimodality.png")
    plt.close(fig)

    # === Plot: feature distributions ===
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), dpi=120)
    plot_specs = [
        ("b", "b-value", axes[0, 0]),
        ("b_drift_per_day", "b-drift /day", axes[0, 1]),
        ("benioff_total_log10", "log10 Benioff", axes[0, 2]),
        ("benioff_curv", "Benioff curvature", axes[1, 0]),
        ("n_above_mc", "n events ≥ Mc", axes[1, 1]),
    ]
    axes[1, 2].axis("off")
    for col, label, ax in plot_specs:
        p = feat_df.loc[feat_df["window_kind"] == "precursor", col].to_numpy(dtype=float)
        q_a = feat_df.loc[feat_df["window_kind"] == "null_A", col].to_numpy(dtype=float)
        q_b = feat_df.loc[feat_df["window_kind"] == "null_B", col].to_numpy(dtype=float)
        all_vals = np.concatenate([p[np.isfinite(p)], q_a[np.isfinite(q_a)], q_b[np.isfinite(q_b)]])
        if len(all_vals) > 0:
            bins = np.linspace(np.min(all_vals), np.max(all_vals), 25)
            ax.hist(q_a[np.isfinite(q_a)], bins=bins, alpha=0.45, color="#7d8aa6",
                    density=True, label=f"null A (n={int(np.isfinite(q_a).sum())})",
                    edgecolor="white", linewidth=0.3)
            ax.hist(q_b[np.isfinite(q_b)], bins=bins, alpha=0.35, color="#999000",
                    density=True, label=f"null B (n={int(np.isfinite(q_b).sum())})",
                    edgecolor="white", linewidth=0.3)
            ax.hist(p[np.isfinite(p)], bins=bins, alpha=0.55, color="#c0392b",
                    density=True, label=f"precursor (n={int(np.isfinite(p).sum())})",
                    edgecolor="white", linewidth=0.4)
        auc_a = auc_table[f"{col}__null_A"]["auc"]
        auc_b = auc_table[f"{col}__null_B"]["auc"]
        ax.set_title(f"{label}\nAUC vs A: {auc_a:.3f}   vs B: {auc_b:.3f}")
        ax.set_xlabel(label)
        ax.set_ylabel("density")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    # AUC table panel
    ax_t = axes[1, 2]
    ax_t.axis("off")
    txt = ("preview AUCs (NOT pre-registered)\n\n"
           f"{'feature':<22s}  {'vs A':>6s}  {'z(A)':>6s}  {'vs B':>6s}  {'z(B)':>6s}\n")
    for col in feat_cols:
        a = auc_table[f"{col}__null_A"]
        b = auc_table[f"{col}__null_B"]
        txt += f"  {col:<20s}  {a['auc']:.3f}  {a['perm_z']:+.2f}  {b['auc']:.3f}  {b['perm_z']:+.2f}\n"
    txt += (f"\nN_precursor:  {int((feat_df['window_kind']=='precursor').sum())}\n"
            f"N_null_A:     {int((feat_df['window_kind']=='null_A').sum())}\n"
            f"N_null_B:     {int((feat_df['window_kind']=='null_B').sum())}\n"
            f"Mc used:      {mc:.2f}\n"
            f"window:       {WINDOW_DAYS} days\n"
            f"3σ pass: CI_lo>0.5 AND |effect|>0.05 AND z>3\n"
            f"Bonferroni α: 0.05/36 = {0.05/36:.2e}")
    ax_t.text(0.0, 0.95, txt, family="monospace", fontsize=8.5, va="top")
    fig.suptitle(f"exp05 California (preview): catalog feature distributions vs Null A and Null B",
                 fontsize=11, y=1.0)
    fig.tight_layout()
    fig.savefig(EXP_DIR / "feature_distributions.png")
    plt.close(fig)
    print(f"[{_ts()}] [plot] wrote bvalue_vs_mc.png, zbz_bimodality.png, feature_distributions.png",
          flush=True)

    # === Headline summary ===
    headline_features = []
    for col in feat_cols:
        for null_kind in ("null_A", "null_B"):
            res = auc_table[f"{col}__{null_kind}"]
            sig3 = res["ci_lo"] > 0.5 and abs(res["auc"] - 0.5) > 0.05 and res["perm_z"] > 3
            sigB = res["perm_p"] < 0.05 / 36
            if sig3 or sigB:
                headline_features.append({
                    "feature": col, "null": null_kind,
                    "auc": res["auc"], "z": res["perm_z"], "p": res["perm_p"],
                    "passes_3sigma": bool(sig3), "passes_bonferroni": bool(sigB),
                })
    print(f"\n[{_ts()}] [headline] features passing 3σ or Bonferroni: {len(headline_features)}", flush=True)
    for f in headline_features:
        print(f"[{_ts()}]   {f['feature']:<24s} vs {f['null']:6s}  AUC={f['auc']:.3f} "
              f"z={f['z']:+.2f}  3σ={'Y' if f['passes_3sigma'] else 'n'} "
              f"Bonf={'Y' if f['passes_bonferroni'] else 'n'}", flush=True)

    summary = {
        "experiment": "exp05_california_features",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "preview_deviations": [
            "Catalog M_min=2.5 (pre-reg says 1.5; preview-only)",
            "Mc plateau check works only on M>=2.5 distribution",
            "Null C (colored noise) skipped — catalog-only",
        ],
        "params": {
            "region": CALIFORNIA.name,
            "bbox": [CALIFORNIA.lat_min, CALIFORNIA.lat_max,
                     CALIFORNIA.lon_min, CALIFORNIA.lon_max],
            "window_days": WINDOW_DAYS,
            "target_M_min": TARGET_M_MIN,
            "Mc_chosen": mc, "Mc_diag": mc_diag,
            "zbz_log_eta_threshold": ZBZ_LOG_ETA_THRESHOLD,
        },
        "counts": {
            "n_total": int(len(df)),
            "n_background": int(decl.is_background.sum()),
            "n_clustered": int((~decl.is_background).sum()),
            "n_targets_all": int(len(targets_all)),
            "n_targets_background": int(len(targets_bg)),
            "n_precursor_kept": len(pre_rows),
            "n_precursor_rejected": rejected,
            "n_null_A": len(null_a_rows),
            "n_null_B": len(null_b_rows),
        },
        "auc_table": auc_table,
        "headline_features": headline_features,
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[{_ts()}] [persist] summary.json", flush=True)
    print(f"[{_ts()}] [exp05] done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
