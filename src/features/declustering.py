"""Zaliapin & Ben-Zion 2013 nearest-neighbor declustering.

For each event j, find the parent event i (with t_i < t_j) that minimizes the
space-time-magnitude proximity:

    eta_ij = T_ij * R_ij
    T_ij   = (t_j - t_i) * 10^(-q * b * M_i)         (years)
    R_ij   = r_ij^d_f * 10^(-p * b * M_i)            (km)

with p + q = 1 (we use p = q = 0.5), b is the Gutenberg-Richter b-value, and
d_f is the fractal dimension of the epicenter distribution (~1.6).

The histogram of log10 eta_ij is bimodal: a "clustered" mode at small eta and
a "background" mode at large eta. We separate them by minimum-density
threshold eta_0 between the two modes (computed from a 1D KDE / Gaussian
mixture). Events with their nearest-neighbor proximity below eta_0 are tagged
"clustered" (offspring); the rest are "background."

References:
- Zaliapin, I. & Ben-Zion, Y. (2013). Earthquake clusters in southern California
  I: Identification and stability. JGR, 118(6), 2847-2864.
- Zaliapin, I., Gabrielov, A., Keilis-Borok, V., & Wong, H. (2008). Clustering
  analysis of seismicity and aftershock identification. PRL, 101(1), 018501.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd


# Constants for ZBZ 2013 (standard published values)
DEFAULT_DF = 1.6     # fractal dimension of the epicenter distribution
DEFAULT_Q = 0.5
DEFAULT_P = 0.5      # p + q = 1


@dataclass
class DeclusterResult:
    is_background: np.ndarray   # bool, length N
    parent_idx: np.ndarray      # int, length N (-1 for events with no parent / earliest)
    log_eta_nn: np.ndarray      # float, length N (log10 of nearest-neighbor eta; -inf for orphans)
    log_T_nn: np.ndarray
    log_R_nn: np.ndarray
    eta_threshold: float        # log10 eta_0 cutoff between clustered and background modes
    n_total: int
    n_background: int
    n_clustered: int
    b_used: float
    d_f: float
    p: float
    q: float


def _haversine_km(lat1: np.ndarray, lon1: np.ndarray, lat2: float, lon2: float) -> np.ndarray:
    """Vectorized great-circle distance in km."""
    R = 6371.0
    lat1r = np.radians(lat1)
    lon1r = np.radians(lon1)
    lat2r = math.radians(lat2)
    lon2r = math.radians(lon2)
    dlat = lat2r - lat1r
    dlon = lon2r - lon1r
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1r) * math.cos(lat2r) * np.sin(dlon / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def _find_eta_threshold(log_eta: np.ndarray, n_grid: int = 400) -> float:
    """Find the trough between the two modes of a bimodal log10(eta) histogram.

    Uses a KDE on a 1D grid; returns the log_eta location of the minimum-density
    point in the central region between the median and a high quantile.
    """
    valid = log_eta[np.isfinite(log_eta)]
    if len(valid) < 50:
        # Not enough data for stable KDE; fall back to a fixed -5 (ZBZ 2013 default).
        return -5.0
    from scipy.stats import gaussian_kde

    kde = gaussian_kde(valid, bw_method=0.20)
    grid = np.linspace(np.percentile(valid, 1), np.percentile(valid, 99), n_grid)
    density = kde(grid)

    # Search for the trough only in the central 30%-90% range. The very-low end
    # is dominated by the cluster mode (we don't want to call its low tail a
    # "trough"); the very-high end is just KDE rolloff.
    lo = int(0.30 * n_grid)
    hi = int(0.90 * n_grid)
    trough_idx = lo + int(np.argmin(density[lo:hi]))
    return float(grid[trough_idx])


def zbz_decluster(
    df: pd.DataFrame,
    b_value: float,
    d_f: float = DEFAULT_DF,
    p: float = DEFAULT_P,
    q: float = DEFAULT_Q,
    eta_threshold: float | None = None,
    log: callable | None = print,
) -> DeclusterResult:
    """Decluster a catalog DataFrame.

    Required DataFrame columns: time (datetime, UTC), latitude, longitude, magnitude.
    Returns a DeclusterResult with one entry per row of df, in the same order.
    """
    if not {"time", "latitude", "longitude", "magnitude"}.issubset(df.columns):
        raise ValueError("df must contain columns: time, latitude, longitude, magnitude")
    if not df["time"].is_monotonic_increasing:
        raise ValueError("df must be sorted by time ascending")

    n = len(df)
    times = df["time"].astype("int64").to_numpy() / 1e9 / (365.25 * 86400.0)  # years since epoch
    lats = df["latitude"].to_numpy()
    lons = df["longitude"].to_numpy()
    mags = df["magnitude"].to_numpy()

    parent = np.full(n, -1, dtype=np.int64)
    log_eta_nn = np.full(n, -np.inf)
    log_T_nn = np.full(n, np.nan)
    log_R_nn = np.full(n, np.nan)

    if log:
        log(f"[zbz] computing nearest-neighbor proximity for N={n} (vectorized, O(N^2))")

    # Vectorize j-by-j: for each j, evaluate eta over all candidate parents i with t_i < t_j.
    # Memory: each j needs an array of length j-1 of eta values, computed on the fly.
    progress_every = max(1, n // 20)
    for j in range(1, n):
        # candidate parents are events 0..j-1 (already sorted by time)
        dt_yr = times[j] - times[:j]
        # Guard: dt_yr should be > 0 since strictly increasing-ish; if zero (duplicates), set tiny epsilon
        dt_yr = np.where(dt_yr <= 0, 1e-9, dt_yr)
        r_km = _haversine_km(lats[:j], lons[:j], lats[j], lons[j])
        # Avoid log(0) for collocated pairs
        r_km = np.where(r_km < 1e-3, 1e-3, r_km)

        T_ij = dt_yr * 10 ** (-q * b_value * mags[:j])
        R_ij = (r_km ** d_f) * 10 ** (-p * b_value * mags[:j])
        eta_ij = T_ij * R_ij

        # Nearest neighbor in proximity (smallest eta)
        i_star = int(np.argmin(eta_ij))
        parent[j] = i_star
        log_eta_nn[j] = math.log10(eta_ij[i_star]) if eta_ij[i_star] > 0 else -np.inf
        log_T_nn[j] = math.log10(T_ij[i_star]) if T_ij[i_star] > 0 else np.nan
        log_R_nn[j] = math.log10(R_ij[i_star]) if R_ij[i_star] > 0 else np.nan

        if log and j % progress_every == 0:
            log(f"[zbz]   {j:>6d} / {n}  ({100 * j / n:5.1f}%)")

    if eta_threshold is None:
        eta_threshold = _find_eta_threshold(log_eta_nn)
        if log:
            log(f"[zbz] auto threshold log10(eta_0) = {eta_threshold:.3f}")
    else:
        if log:
            log(f"[zbz] user threshold log10(eta_0) = {eta_threshold:.3f}")

    is_background = log_eta_nn >= eta_threshold
    is_background[0] = True  # first event has no parent; treat as background by convention

    if log:
        log(
            f"[zbz] background = {int(is_background.sum())} / {n}  "
            f"({100 * is_background.mean():.1f}%); clustered = {int((~is_background).sum())}"
        )

    return DeclusterResult(
        is_background=is_background,
        parent_idx=parent,
        log_eta_nn=log_eta_nn,
        log_T_nn=log_T_nn,
        log_R_nn=log_R_nn,
        eta_threshold=eta_threshold,
        n_total=n,
        n_background=int(is_background.sum()),
        n_clustered=int((~is_background).sum()),
        b_used=b_value,
        d_f=d_f,
        p=p,
        q=q,
    )
