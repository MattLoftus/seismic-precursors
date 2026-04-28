"""Spectral slope of a waveform — Brune 1970 high-frequency falloff.

For a single ground-velocity record, compute the power spectral density and
fit a log-log slope in a chosen frequency band (default 0.5-20 Hz, capped by
fs/3 to dodge anti-alias-filter artifacts at the upper edge).

Brune (1970) proposed that earthquake source spectra fall off as f^(-2)
above the corner frequency f_c; ambient noise typically follows ~f^(-1) to
f^(-1.5) flicker noise. Drift in the spectral slope thus reflects shifting
balance between source-radiation, scattering, and noise content.

This is a per-WINDOW scalar; for a 30-day pre-event window it should be
computed on a representative subsample of the continuous waveform (e.g.,
hourly Welch estimates averaged) — full continuous compute is expensive.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.signal import welch


def spectral_slope(
    data: np.ndarray,
    fs: float,
    fmin: float = 0.5,
    fmax: float = 20.0,
    nperseg: int | None = None,
) -> dict:
    """Fit log10 P(f) = a + slope * log10(f) over [fmin, min(fmax, fs/3)].

    Returns dict with: slope, intercept, r_squared, fmin_used, fmax_used,
    n_freq_bins.
    """
    fmax_used = min(fmax, fs / 3.0)
    if fmax_used <= fmin:
        return {
            "slope": float("nan"),
            "intercept": float("nan"),
            "r_squared": float("nan"),
            "fmin_used": fmin,
            "fmax_used": fmax_used,
            "n_freq_bins": 0,
            "reason": f"fmax_used {fmax_used} <= fmin {fmin}",
        }

    if nperseg is None:
        # Aim for ~0.05 Hz frequency resolution if signal allows
        target = int(round(fs / 0.05))
        nperseg = min(len(data), max(64, target))

    f, Pxx = welch(np.asarray(data, dtype=float), fs=fs, nperseg=nperseg,
                   noverlap=nperseg // 2, scaling="density")
    band = (f >= fmin) & (f <= fmax_used) & (Pxx > 0)
    if int(band.sum()) < 5:
        return {
            "slope": float("nan"),
            "intercept": float("nan"),
            "r_squared": float("nan"),
            "fmin_used": fmin,
            "fmax_used": fmax_used,
            "n_freq_bins": int(band.sum()),
            "reason": "fewer than 5 positive PSD bins in the chosen band",
        }

    log_f = np.log10(f[band])
    log_P = np.log10(Pxx[band])
    slope, intercept = np.polyfit(log_f, log_P, deg=1)
    pred = slope * log_f + intercept
    ss_res = float(((log_P - pred) ** 2).sum())
    ss_tot = float(((log_P - log_P.mean()) ** 2).sum())
    r_sq = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "r_squared": float(r_sq),
        "fmin_used": float(fmin),
        "fmax_used": float(fmax_used),
        "n_freq_bins": int(band.sum()),
    }
