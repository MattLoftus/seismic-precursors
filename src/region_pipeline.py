"""Per-region pipeline: catalog → Mc → ZBZ → precursor/null windows → features → AUC.

Factored out of `experiments/exp05_california_features/run.py` so exp06+ can run
the same protocol on multiple regions without duplication. The main entry point
is `run_region_pipeline(region, start, end, params, catalog_cache) -> RegionResult`.

Implements the locked protocol from `papers/pre_registration.md` at SHA
a4f1c6fa093ecea8bd00077fb8dbcab008c70eaa, with the same declared deviations
documented at the experiment level (catalog M_min and Null C handling).
"""
from __future__ import annotations

import datetime as dt
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from src.data import fetch_comcat_catalog_bbox, fetch_isc_catalog_bbox
from src.features.benioff import benioff_features
from src.features.bvalue import (
    aki_bvalue,
    bvalue_drift,
    bvalue_with_bootstrap,
    magnitude_of_completeness_max_curvature,
)
from src.features.declustering import zbz_decluster
from src.regions import Region


@dataclass
class PipelineParams:
    """All knobs in one place. Defaults match `papers/pre_registration.md` v1."""
    window_days: int = 30
    target_m_min: float = 4.5
    catalog_m_min: float = 2.5            # PREVIEW deviation; pre-reg = 1.5
    catalog_source: str = "ComCat"        # "ComCat" (v1) or "ISC" (PRA-2 amendment 1)
    zbz_log_eta_threshold: float = -5.0
    zbz_b_input: float = 1.0
    n_drift_subwindows: int = 3
    min_events_per_subwindow: int = 10
    null_buffer_days_before: int = 30     # set to 0 for PRA-2 amendment 2
    null_buffer_days_after: int = 60
    n_null_a: int = 200
    n_null_b: int = 200
    n_boot: int = 1000
    n_perm: int = 1000
    random_seed: int = 42
    chunk_years: int = 1                  # ComCat per-call limit ~20k events
    mc_grid: tuple[float, ...] = (2.5, 2.7, 2.9, 3.1, 3.3, 3.5)
    min_kept_precursor_per_region: int = 8  # PRA-2 amendment 3


FEATURE_COLS = ["b", "b_drift_per_day", "benioff_total_log10", "benioff_curv", "n_above_mc"]


@dataclass
class RegionResult:
    region_name: str
    mc: float
    mc_diag: dict
    n_total: int
    n_background: int
    n_targets_all: int
    n_targets_background: int
    n_precursor: int
    n_precursor_rejected: int
    n_null_a: int
    n_null_b: int
    feat_df: pd.DataFrame                 # all windows, with window_kind ∈ {precursor, null_A, null_B}
    auc_table: dict                       # {f"{feat}__{null_kind}": auc_with_significance dict}


def estimate_mc_with_plateau_check(magnitudes: np.ndarray,
                                   mc_grid: tuple[float, ...]) -> tuple[float, dict]:
    """Pre-reg §1.2: max-curvature + Woessner +0.20 + plateau iteration."""
    mc_raw = magnitude_of_completeness_max_curvature(magnitudes)
    mc_woessner = mc_raw + 0.20
    bs = []
    for mc in mc_grid:
        if (magnitudes >= mc).sum() < 100:
            bs.append(float("nan"))
            continue
        b, _, _ = aki_bvalue(magnitudes, mc=mc)
        bs.append(float(b))
    chosen = None
    for i in range(1, len(bs)):
        if not (math.isfinite(bs[i]) and math.isfinite(bs[i - 1])):
            continue
        if abs(bs[i] - bs[i - 1]) < 0.02:
            chosen = mc_grid[i]
            break
    if chosen is None:
        for i in reversed(range(len(bs))):
            if math.isfinite(bs[i]):
                chosen = mc_grid[i]
                break
    if chosen is None:
        chosen = mc_woessner
    chosen = max(chosen, mc_woessner)
    return float(chosen), {
        "mc_raw": float(mc_raw),
        "mc_woessner": float(mc_woessner),
        "mc_grid": list(mc_grid),
        "b_at_mc_grid": bs,
    }


def features_for_window(df: pd.DataFrame, t_start, t_end, mc: float,
                        params: PipelineParams) -> dict:
    mask = (df["time"] >= t_start) & (df["time"] < t_end)
    sub = df.loc[mask]
    n = int(len(sub))
    if n > 0:
        ts_seconds = (sub["time"] - t_start).dt.total_seconds().to_numpy()
        mags = sub["magnitude"].to_numpy()
    else:
        ts_seconds = np.array([])
        mags = np.array([])

    out = {"t_start": t_start, "t_end": t_end, "n_events": n,
           "n_above_mc": int((mags >= mc).sum()) if n else 0}
    if (mags >= mc).sum() >= 30:
        res = bvalue_with_bootstrap(mags, mc=mc, n_boot=0, mc_correction=0.0)
        out["b"] = float(res.b)
    else:
        out["b"] = float("nan")
    drift = bvalue_drift(
        ts_seconds, mags, mc=mc,
        n_subwindows=params.n_drift_subwindows,
        min_events_per_subwindow=params.min_events_per_subwindow,
    )
    out["b_drift_per_day"] = float(drift["drift_slope_per_day"])
    out["b_drift_ok"] = bool(drift["ok"])
    bf = benioff_features(ts_seconds, mags,
                          window_seconds=params.window_days * 86400.0)
    out["benioff_total_log10"] = (
        float(np.log10(bf["benioff_total"])) if bf["benioff_total"] > 0 else float("nan")
    )
    out["benioff_curv"] = float(bf["benioff_curv"])
    return out


def sample_null_windows(df: pd.DataFrame, target_times, n: int, seed: int,
                        kind: str, params: PipelineParams) -> list[tuple]:
    """kind='A': aftershock-free per pre-reg §5.2.  kind='B': minimal exclusion per §5.3."""
    rng = np.random.default_rng(seed)
    t_first = df["time"].iloc[0]
    t_last = df["time"].iloc[-1] - pd.Timedelta(days=params.window_days)
    target_times = pd.to_datetime(target_times, utc=True)
    if kind == "A":
        forbid_starts = target_times - pd.Timedelta(
            days=params.window_days + params.null_buffer_days_before
        )
        forbid_ends = target_times + pd.Timedelta(days=params.null_buffer_days_after)
    elif kind == "B":
        forbid_starts = target_times - pd.Timedelta(days=params.window_days)
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
        t_end = t_start + pd.Timedelta(days=params.window_days)
        overlap = ((t_start <= forbid_ends) & (t_end >= forbid_starts)).any()
        if not overlap:
            windows.append((t_start, t_end))
    return windows


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


def auc_with_significance(p: np.ndarray, q: np.ndarray,
                          n_boot: int, n_perm: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    auc = per_feature_auc(p, q)
    p_fin = p[np.isfinite(p)]
    q_fin = q[np.isfinite(q)]
    if len(p_fin) < 2 or len(q_fin) < 2:
        return {"auc": auc, "ci_lo": float("nan"), "ci_hi": float("nan"),
                "perm_z": float("nan"), "perm_p": float("nan"),
                "n_p": int(len(p_fin)), "n_q": int(len(q_fin))}
    boot = np.empty(n_boot)
    for i in range(n_boot):
        ps = rng.choice(p_fin, size=len(p_fin), replace=True)
        qs = rng.choice(q_fin, size=len(q_fin), replace=True)
        boot[i] = per_feature_auc(ps, qs)
    ci_lo, ci_hi = np.nanpercentile(boot, [2.5, 97.5])
    pooled = np.concatenate([p_fin, q_fin])
    perm = np.empty(n_perm)
    for i in range(n_perm):
        rng.shuffle(pooled)
        perm[i] = per_feature_auc(pooled[:len(p_fin)], pooled[len(p_fin):])
    perm_mean = float(np.nanmean(perm))
    perm_std = float(np.nanstd(perm))
    z = (auc - perm_mean) / perm_std if perm_std > 0 else float("nan")
    from scipy.stats import norm
    p_val = 2 * (1 - norm.cdf(abs(z)))
    return {"auc": float(auc), "ci_lo": float(ci_lo), "ci_hi": float(ci_hi),
            "perm_mean": perm_mean, "perm_std": perm_std,
            "perm_z": float(z), "perm_p": float(p_val),
            "n_p": int(len(p_fin)), "n_q": int(len(q_fin))}


def run_region_pipeline(
    region: Region,
    start: dt.datetime,
    end: dt.datetime,
    params: PipelineParams,
    catalog_cache: Path | None = None,
    log: callable | None = print,
) -> RegionResult:
    """Apply the full pre-reg protocol to one region; return RegionResult."""
    if log:
        log(f"[pipe] === region {region.name} (catalog={params.catalog_source}) ===")

    if params.catalog_source.upper() == "ISC":
        fetcher = fetch_isc_catalog_bbox
    elif params.catalog_source.upper() == "COMCAT":
        fetcher = fetch_comcat_catalog_bbox
    else:
        raise ValueError(f"unknown catalog_source {params.catalog_source!r}; "
                         f"expected 'ComCat' or 'ISC'")
    df = fetcher(
        lat_min=region.lat_min, lat_max=region.lat_max,
        lon_min=region.lon_min, lon_max=region.lon_max,
        start=start, end=end, m_min=params.catalog_m_min,
        chunk_years=params.chunk_years, cache_path=catalog_cache, log=log,
    )
    if log:
        log(f"[pipe] {region.name}: {len(df)} events")

    mc, mc_diag = estimate_mc_with_plateau_check(
        df["magnitude"].to_numpy(), params.mc_grid
    )
    if log:
        log(f"[pipe] {region.name}: Mc={mc:.2f} "
            f"(raw={mc_diag['mc_raw']:.2f}, woessner={mc_diag['mc_woessner']:.2f})")

    decl = zbz_decluster(
        df, b_value=params.zbz_b_input,
        eta_threshold=params.zbz_log_eta_threshold,
        log=None,                          # silence per-region progress to keep output readable
    )
    df_bg = df.loc[decl.is_background].reset_index(drop=True)
    if log:
        log(f"[pipe] {region.name}: ZBZ background {len(df_bg)}/{len(df)} "
            f"({100*len(df_bg)/max(1,len(df)):.1f}%)")

    targets_all = df.loc[df["magnitude"] >= params.target_m_min].reset_index(drop=True)
    targets_bg = df_bg.loc[df_bg["magnitude"] >= params.target_m_min].reset_index(drop=True)
    if log:
        log(f"[pipe] {region.name}: targets all={len(targets_all)} background={len(targets_bg)}")

    # Precursor windows with overlap rejection per pre-reg §5.1 (PRA-2 amendment 2:
    # forbidden zone is [t', t'+60] when null_buffer_days_before=0, dropping the pre-event side).
    forbid_starts = targets_bg["time"] - pd.Timedelta(days=params.null_buffer_days_before)
    forbid_ends = targets_bg["time"] + pd.Timedelta(days=params.null_buffer_days_after)
    pre_rows = []
    rejected = 0
    for _, row in targets_bg.iterrows():
        t_end = row["time"]
        t_start = t_end - pd.Timedelta(days=params.window_days)
        own = forbid_starts == t_end - pd.Timedelta(days=params.null_buffer_days_before)
        other_starts = forbid_starts[~own]
        other_ends = forbid_ends[~own]
        overlap = ((t_start <= other_ends) & (t_end >= other_starts)).any()
        if overlap:
            rejected += 1
            continue
        feat = features_for_window(df, t_start, t_end, mc, params)
        feat["target_M"] = float(row["magnitude"])
        feat["target_time"] = t_end
        feat["window_kind"] = "precursor"
        feat["region"] = region.name
        pre_rows.append(feat)
    if log:
        log(f"[pipe] {region.name}: precursor windows kept {len(pre_rows)} (rejected {rejected})  "
            f"[overlap_buffer_pre={params.null_buffer_days_before}d post={params.null_buffer_days_after}d]")

    # Null A and Null B
    null_a = sample_null_windows(df, targets_bg["time"].tolist(), params.n_null_a,
                                 params.random_seed, kind="A", params=params)
    null_b = sample_null_windows(df, targets_all["time"].tolist(), params.n_null_b,
                                 params.random_seed + 1, kind="B", params=params)
    if log:
        log(f"[pipe] {region.name}: null_A {len(null_a)}, null_B {len(null_b)}")

    null_a_rows = []
    for t_start, t_end in null_a:
        feat = features_for_window(df, t_start, t_end, mc, params)
        feat["target_M"] = float("nan"); feat["target_time"] = pd.NaT
        feat["window_kind"] = "null_A"; feat["region"] = region.name
        null_a_rows.append(feat)
    null_b_rows = []
    for t_start, t_end in null_b:
        feat = features_for_window(df, t_start, t_end, mc, params)
        feat["target_M"] = float("nan"); feat["target_time"] = pd.NaT
        feat["window_kind"] = "null_B"; feat["region"] = region.name
        null_b_rows.append(feat)

    all_rows = pre_rows + null_a_rows + null_b_rows
    if not all_rows:
        feat_df = pd.DataFrame(columns=["window_kind", "region"] + FEATURE_COLS)
    else:
        feat_df = pd.DataFrame(all_rows)
        if "window_kind" not in feat_df.columns:
            feat_df["window_kind"] = []  # defensive

    # Per-feature AUC for this region. If no precursor windows, all AUCs are NaN.
    auc_table = {}
    nan_res = {"auc": float("nan"), "ci_lo": float("nan"), "ci_hi": float("nan"),
               "perm_z": float("nan"), "perm_p": float("nan"), "n_p": 0, "n_q": 0}
    for null_kind in ("null_A", "null_B"):
        for col in FEATURE_COLS:
            if not pre_rows:
                auc_table[f"{col}__{null_kind}"] = dict(nan_res)
                continue
            p = feat_df.loc[feat_df["window_kind"] == "precursor", col].to_numpy(dtype=float)
            q = feat_df.loc[feat_df["window_kind"] == null_kind, col].to_numpy(dtype=float)
            res = auc_with_significance(p, q, params.n_boot, params.n_perm, params.random_seed)
            auc_table[f"{col}__{null_kind}"] = res
            if log:
                log(f"[pipe] {region.name:<11s} [auc] {col:<22s} vs {null_kind:6s}  "
                    f"AUC={res['auc']:.3f}  z={res['perm_z']:+.2f}  p={res['perm_p']:.2e}")
    if not pre_rows and log:
        log(f"[pipe] {region.name}: 0 precursor windows — all AUCs are NaN; "
            f"region effectively excluded from macro pool")

    return RegionResult(
        region_name=region.name,
        mc=mc, mc_diag=mc_diag,
        n_total=int(len(df)), n_background=int(decl.is_background.sum()),
        n_targets_all=int(len(targets_all)), n_targets_background=int(len(targets_bg)),
        n_precursor=len(pre_rows), n_precursor_rejected=rejected,
        n_null_a=len(null_a_rows), n_null_b=len(null_b_rows),
        feat_df=feat_df, auc_table=auc_table,
    )
