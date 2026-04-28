"""Hilbert-Huang Transform: EMD decomposition + instantaneous frequency.

Pipeline (Huang+ 1998):
    1. Empirical Mode Decomposition (EMD): decompose signal into IMFs.
    2. Hilbert transform each IMF -> analytic signal.
    3. Instantaneous frequency = d phase / dt / (2 pi).

For the precursor feature (PLAYBOOK §5.2, "HHT peak drift"), we track the
mean instantaneous frequency of the FIRST IMF over a sliding window.
IMF1 carries the highest-frequency component of the signal; its IF drift
reflects the dominant high-frequency mode shifting (e.g., from ~10 Hz
microseismic noise to ~5 Hz body-wave energy near a transit).
"""
from __future__ import annotations

import numpy as np
from scipy.signal import hilbert


def empirical_mode_decomposition(
    data: np.ndarray,
    max_imfs: int = 6,
) -> np.ndarray:
    """Run EMD via the `emd` package; return IMFs as (n_imfs, n_samples) array.

    The `emd` package (Quinn et al.) is the well-maintained Python EMD
    implementation. It expects a 1D array.
    """
    import emd as emd_pkg

    imfs = emd_pkg.sift.sift(np.asarray(data, dtype=float), max_imfs=max_imfs)
    # emd returns (n_samples, n_imfs); transpose to (n_imfs, n_samples)
    return imfs.T


def instantaneous_frequency(imf: np.ndarray, fs: float) -> np.ndarray:
    """Instantaneous frequency of an IMF via analytic-signal phase derivative.

    Returns a (n_samples - 1,) array of IF values in Hz.
    """
    analytic = hilbert(np.asarray(imf, dtype=float))
    phase = np.unwrap(np.angle(analytic))
    return np.diff(phase) * fs / (2 * np.pi)


def imf1_if_series(
    data: np.ndarray,
    fs: float,
    max_imfs: int = 6,
) -> tuple[np.ndarray, np.ndarray]:
    """Mean IF of the first IMF, using all of `data` at once.

    Returns (imfs, if1_hz). imfs has shape (n_imfs, n_samples). if1_hz has
    shape (n_samples - 1,).
    """
    imfs = empirical_mode_decomposition(data, max_imfs=max_imfs)
    if len(imfs) == 0:
        return imfs, np.array([])
    if1 = instantaneous_frequency(imfs[0], fs=fs)
    # Clip pathological IF values (Hilbert can produce negative or huge IFs at
    # signal edges or near zero-crossings of poorly-defined IMFs).
    if1 = np.clip(if1, 0.0, fs / 2.0)
    return imfs, if1
