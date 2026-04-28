"""Benioff strain accumulation: sum sqrt(E) over events in a window.

The seismic-energy estimate per event is the Gutenberg-Richter form
    log10(E_J) = 1.5 M + 4.8     (joules)
giving E_M=4.5 ~ 6.3e10 J, E_M=6.0 ~ 1.0e13 J. Benioff (1951) defined the
"strain release" rate as Σ √E per unit time. Bowman+ 1998 promoted ΣE^(1/2)
as a precursor metric (the "AMR" — Accelerating Moment Release — hypothesis).

We compute three numbers per window:
    benioff_total  = Σ √E_i                (the cumulative scalar)
    benioff_rate   = Σ √E_i / window_seconds
    benioff_curv   = curvature of the cumulative ΣE^(1/2)(t) curve
                     (positive curvature → accelerating release; the AMR signature)
"""
from __future__ import annotations

import math

import numpy as np


def event_energy_joules(magnitude: float | np.ndarray) -> float | np.ndarray:
    """Gutenberg-Richter energy: log10(E) = 1.5 M + 4.8."""
    return 10 ** (1.5 * np.asarray(magnitude) + 4.8)


def benioff_features(
    times_seconds: np.ndarray,
    magnitudes: np.ndarray,
    window_seconds: float,
) -> dict:
    """Compute Benioff strain features for a window of events.

    times_seconds: monotonically increasing event times (any origin).
    magnitudes:    same length, event magnitudes.
    window_seconds: total length of the window for rate normalization.

    Returns dict with: benioff_total, benioff_rate (per s), benioff_curv,
    cum_strain (cumulative time series), n_events.
    """
    if len(times_seconds) != len(magnitudes):
        raise ValueError("times_seconds and magnitudes must have same length")
    if window_seconds <= 0:
        raise ValueError("window_seconds must be positive")

    n = int(len(times_seconds))
    if n == 0:
        return {
            "benioff_total": 0.0,
            "benioff_rate": 0.0,
            "benioff_curv": 0.0,
            "cum_strain": np.array([]),
            "n_events": 0,
        }

    sqrt_E = np.sqrt(event_energy_joules(magnitudes))
    cum = np.cumsum(sqrt_E)

    # Curvature: fit a quadratic c2 t^2 + c1 t + c0 to (t, cum) and report c2.
    # Positive c2 = accelerating; negative = decelerating.
    t = np.asarray(times_seconds, dtype=float)
    if n >= 3 and t[-1] > t[0]:
        # Normalize time so the fit is well-conditioned
        t_norm = (t - t[0]) / (t[-1] - t[0])
        coefs = np.polyfit(t_norm, cum, deg=2)
        # coefs are in ORDER [c2, c1, c0]; rescale c2 from per-(unit-norm-t)^2 back
        # to per-second^2 to be consistent with rate units.
        scale_factor = (t[-1] - t[0]) ** 2
        curv_per_s2 = float(coefs[0] / scale_factor)
    else:
        curv_per_s2 = 0.0

    return {
        "benioff_total": float(cum[-1]),
        "benioff_rate": float(cum[-1] / window_seconds),
        "benioff_curv": curv_per_s2,
        "cum_strain": cum,
        "n_events": n,
    }
