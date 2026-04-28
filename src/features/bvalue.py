"""Gutenberg-Richter b-value via Aki 1965 MLE, with Mc via Wiemer & Wyss 2000 max-curvature.

References:
- Aki, K. (1965). Maximum likelihood estimate of b in the formula log N = a - bM
  and its confidence limits. Bull. Earthq. Res. Inst., 43, 237-239.
- Wiemer, S. & Wyss, M. (2000). Minimum magnitude of completeness in earthquake
  catalogs: examples from Alaska, the western United States, and Japan.
  BSSA, 90(4), 859-869.
- Marzocchi, F. & Sandri, L. (2003). A review and new insights on the estimation
  of the b-value and its uncertainty. Annals of Geophysics, 46(6), 1271-1282.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


DM_DEFAULT = 0.1  # standard catalog magnitude binning


@dataclass
class BValueResult:
    b: float
    b_se: float            # 1-sigma standard error (Aki / Marzocchi-Sandri)
    b_boot_lo: float       # 95% bootstrap CI lower
    b_boot_hi: float       # 95% bootstrap CI upper
    a: float               # log10 N(M=0) intercept
    mc: float              # magnitude of completeness
    n_above_mc: int        # number of events above Mc
    n_total: int           # total catalog events used
    method_mc: str         # method tag for Mc
    dm: float              # magnitude binning


def magnitude_of_completeness_max_curvature(magnitudes: np.ndarray, dm: float = DM_DEFAULT) -> float:
    """Wiemer & Wyss 2000 max-curvature: Mc = mode of the binned frequency-magnitude distribution.

    A common +0.2 correction (Woessner & Wiemer 2005) gives a less biased Mc estimate.
    We return the raw mode here; the caller may add the correction.
    """
    if len(magnitudes) == 0:
        return float("nan")
    m_min = math.floor(np.min(magnitudes) / dm) * dm
    m_max = math.ceil(np.max(magnitudes) / dm) * dm + dm
    edges = np.arange(m_min, m_max + dm / 2, dm)
    counts, _ = np.histogram(magnitudes, bins=edges)
    if counts.sum() == 0:
        return float("nan")
    mode_idx = int(np.argmax(counts))
    # Bin center for the modal bin
    return float(edges[mode_idx] + dm / 2)


def aki_bvalue(magnitudes: np.ndarray, mc: float, dm: float = DM_DEFAULT) -> tuple[float, float, int]:
    """Aki 1965 MLE b-value with Marzocchi-Sandri 2003 binning correction.

    b = log10(e) / (<M> - (Mc - dm/2))
    sigma_b = b / sqrt(N)  (Aki uncertainty; Shi & Bolt 1982 gives a tighter form)

    Returns (b, sigma_b, n_used).
    """
    above = magnitudes[magnitudes >= mc - 1e-9]
    n = int(len(above))
    if n < 2:
        return float("nan"), float("nan"), n
    mean_m = float(np.mean(above))
    denom = mean_m - (mc - dm / 2)
    if denom <= 0:
        return float("nan"), float("nan"), n
    b = math.log10(math.e) / denom
    # Shi & Bolt 1982 standard error (sharper than Aki's sigma_b = b/sqrt(N))
    var_m = float(np.var(above, ddof=1)) if n > 1 else 0.0
    if var_m > 0:
        sigma_b = 2.30 * (b ** 2) * math.sqrt(var_m / n)
    else:
        sigma_b = b / math.sqrt(n)
    return b, sigma_b, n


def aki_avalue(magnitudes: np.ndarray, b: float, mc: float) -> float:
    """log10 N(M >= 0) intercept assuming N(M >= Mc) is observed."""
    n_above = int(np.sum(magnitudes >= mc - 1e-9))
    if n_above == 0 or not math.isfinite(b):
        return float("nan")
    return math.log10(n_above) + b * mc


def bvalue_with_bootstrap(
    magnitudes: np.ndarray,
    mc: float | None = None,
    dm: float = DM_DEFAULT,
    n_boot: int = 1000,
    rng_seed: int = 42,
    mc_correction: float = 0.2,
) -> BValueResult:
    """Compute b-value with bootstrap CI. If mc is None, estimate via max-curvature + correction.

    The mc_correction (default +0.2) follows Woessner & Wiemer 2005 to debias the
    raw max-curvature estimate.
    """
    magnitudes = np.asarray(magnitudes, dtype=float)
    if mc is None:
        mc_raw = magnitude_of_completeness_max_curvature(magnitudes, dm=dm)
        mc = mc_raw + mc_correction
        method_mc = f"max-curvature + {mc_correction:+.2f} (Woessner-Wiemer 2005)"
    else:
        method_mc = "user-specified"

    b, b_se, n_above = aki_bvalue(magnitudes, mc=mc, dm=dm)
    a = aki_avalue(magnitudes, b=b, mc=mc)

    above = magnitudes[magnitudes >= mc - 1e-9]
    rng = np.random.default_rng(rng_seed)
    if len(above) >= 2 and n_boot > 0:
        boot_bs = np.empty(n_boot)
        for i in range(n_boot):
            sample = rng.choice(above, size=len(above), replace=True)
            bb, _, _ = aki_bvalue(sample, mc=mc, dm=dm)
            boot_bs[i] = bb
        lo, hi = np.nanpercentile(boot_bs, [2.5, 97.5])
    else:
        lo = hi = float("nan")

    return BValueResult(
        b=b,
        b_se=b_se,
        b_boot_lo=float(lo),
        b_boot_hi=float(hi),
        a=a,
        mc=mc,
        n_above_mc=n_above,
        n_total=int(len(magnitudes)),
        method_mc=method_mc,
        dm=dm,
    )


def bvalue_drift(
    times_seconds: np.ndarray,
    magnitudes: np.ndarray,
    mc: float,
    n_subwindows: int = 5,
    min_events_per_subwindow: int = 30,
    dm: float = DM_DEFAULT,
) -> dict:
    """Time-resolved b within a window: split into n_subwindows equal-time chunks
    and fit b in each. Returns the linear slope (b vs sub-window center time)
    along with the per-sub-window b values.

    The slope is the precursor feature ("b-value drift", Gulia & Wiemer 2019
    style — the sign and magnitude of b-change leading up to a target event).

    times_seconds: float array of event times in seconds (any monotonic origin).
    magnitudes:    float array of magnitudes, same length.

    Returns dict with: drift_slope (b per second), drift_slope_per_day, b_series,
    t_centers, n_above_mc, ok (bool flag — False if too sparse).
    """
    times_seconds = np.asarray(times_seconds, dtype=float)
    magnitudes = np.asarray(magnitudes, dtype=float)
    if len(times_seconds) != len(magnitudes):
        raise ValueError("times_seconds and magnitudes must have same length")

    above = magnitudes >= (mc - 1e-9)
    n_above = int(above.sum())
    if n_above < n_subwindows * min_events_per_subwindow:
        return {
            "ok": False,
            "drift_slope_per_s": float("nan"),
            "drift_slope_per_day": float("nan"),
            "b_series": np.array([]),
            "t_centers": np.array([]),
            "n_above_mc": n_above,
            "reason": f"insufficient events above Mc ({n_above} < "
                      f"{n_subwindows * min_events_per_subwindow})",
        }

    t_min, t_max = float(times_seconds.min()), float(times_seconds.max())
    edges = np.linspace(t_min, t_max, n_subwindows + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    bs = np.full(n_subwindows, np.nan)
    for i in range(n_subwindows):
        mask = (times_seconds >= edges[i]) & (times_seconds < edges[i + 1] if i < n_subwindows - 1
                                               else times_seconds <= edges[i + 1])
        if int((mask & above).sum()) >= min_events_per_subwindow:
            bb, _, _ = aki_bvalue(magnitudes[mask], mc=mc, dm=dm)
            bs[i] = bb

    valid = np.isfinite(bs)
    if valid.sum() < 2:
        return {
            "ok": False,
            "drift_slope_per_s": float("nan"),
            "drift_slope_per_day": float("nan"),
            "b_series": bs,
            "t_centers": centers,
            "n_above_mc": n_above,
            "reason": f"only {int(valid.sum())} sub-windows had >= "
                      f"{min_events_per_subwindow} above-Mc events",
        }

    # Least-squares slope of b vs time
    slope, intercept = np.polyfit(centers[valid], bs[valid], deg=1)
    return {
        "ok": True,
        "drift_slope_per_s": float(slope),
        "drift_slope_per_day": float(slope * 86400.0),
        "intercept": float(intercept),
        "b_series": bs,
        "t_centers": centers,
        "n_above_mc": n_above,
    }
