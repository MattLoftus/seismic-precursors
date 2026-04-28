"""Spectral (Shannon) entropy of whitened-spectrogram time slices.

For each short-time window, compute power spectrum; whiten by dividing the
spectrum by its mean (kills amplitude variation, leaves spectral *shape*);
then Shannon entropy on the normalized power distribution H(t) = -sum_f
p(f, t) log2 p(f, t).

A signal whose energy is uniformly spread across frequencies has high
entropy (~log2(N_freq_bins)); a signal whose energy is concentrated in
narrow frequency bands has low entropy. Earthquakes shift entropy DOWN
relative to ambient noise (which is closer to white).

Adapted from the ambient-noise tomography preprocessing chain
(Bensen et al. 2007, GJI 169:1239) but applied as a precursor feature
following the framing in PLAYBOOK §5.2.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.signal import spectrogram


def spectral_entropy_series(
    data: np.ndarray,
    fs: float,
    window_seconds: float = 4.0,
    overlap: float = 0.5,
    fmin: float = 0.5,
    fmax: float = 20.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Time series of Shannon entropy over a whitened spectrogram.

    Returns (t_centers_s, entropy_bits). entropy_bits is in bits (log2).
    """
    nperseg = int(round(window_seconds * fs))
    noverlap = int(round(overlap * nperseg))
    f, t, Sxx = spectrogram(
        np.asarray(data, dtype=float),
        fs=fs,
        nperseg=nperseg,
        noverlap=noverlap,
        scaling="density",
        mode="psd",
    )
    band = (f >= fmin) & (f <= fmax)
    f_band = f[band]
    P = Sxx[band, :]                            # (n_freq, n_time)
    # Whiten: divide each time-slice by its mean spectrum value
    P_mean = P.mean(axis=0, keepdims=True)
    P_mean = np.where(P_mean > 0, P_mean, 1e-30)
    Pw = P / P_mean
    # Convert to a probability distribution per time-slice
    Pw_sum = Pw.sum(axis=0, keepdims=True)
    Pw_sum = np.where(Pw_sum > 0, Pw_sum, 1e-30)
    p = Pw / Pw_sum
    # Shannon entropy (in bits)
    p_safe = np.clip(p, 1e-30, 1.0)
    H = -(p_safe * np.log2(p_safe)).sum(axis=0)
    H_max = math.log2(max(2, len(f_band)))     # max entropy at this binning
    return t, H, H_max


def windowed_mean_entropy(
    data: np.ndarray,
    fs: float,
    *,
    window_seconds: float = 4.0,
    overlap: float = 0.5,
    fmin: float = 0.5,
    fmax: float = 20.0,
) -> tuple[float, float]:
    """Convenience: scalar (mean, std) of the entropy time series for a window."""
    _, H, _ = spectral_entropy_series(
        data, fs=fs,
        window_seconds=window_seconds,
        overlap=overlap,
        fmin=fmin, fmax=fmax,
    )
    return float(np.mean(H)), float(np.std(H))
