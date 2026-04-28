"""Repeating-event detection via waveform cross-correlation.

Following Nadeau & McEvilly (1999), repeaters are pairs of events at the
same location with near-identical waveforms (cross-correlation > 0.95 on
broadband). The RATE of repeaters per unit time is the precursor feature.

This module contains the cross-correlation kernel + a synthetic test path.
The full catalog application — for each event in a window, fetch its
waveform, cross-correlate against all earlier events of similar magnitude
within ~5 km — requires waveform-pull infrastructure that we will add in
Session 5+ when we move to Round C cross-regional evaluation.

The synthetic path lets us validate that the kernel correctly identifies
near-duplicate waveforms in noise.
"""
from __future__ import annotations

import numpy as np
from scipy.signal import correlate


def normalized_xcorr(a: np.ndarray, b: np.ndarray) -> float:
    """Maximum normalized cross-correlation between two equal-length signals.

    Returns the peak value of the (full-mode) cross-correlation, normalized
    by the product of L2 norms. Range: [-1, 1].
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch: {a.shape} vs {b.shape}")
    a = a - a.mean()
    b = b - b.mean()
    norm = np.linalg.norm(a) * np.linalg.norm(b)
    if norm == 0:
        return 0.0
    return float(np.max(correlate(a, b, mode="full")) / norm)


def count_repeaters(
    waveforms: list[np.ndarray],
    threshold: float = 0.95,
) -> dict:
    """Count event pairs (i, j) with i < j whose normalized cross-correlation
    is >= threshold. Returns count, indices of repeater pairs, and rate
    (pairs / N events).

    waveforms: list of N equal-length 1-D arrays. For Round C cross-regional
    application this will be a list per window of pre-aligned event waveforms.
    """
    n = len(waveforms)
    pair_indices: list[tuple[int, int]] = []
    if n < 2:
        return {
            "n_pairs": 0,
            "pair_indices": pair_indices,
            "rate_per_event": 0.0,
            "n_events": n,
            "threshold": threshold,
        }
    # O(N^2) — acceptable for ~50-200 events per window
    for i in range(n):
        for j in range(i + 1, n):
            cc = normalized_xcorr(waveforms[i], waveforms[j])
            if cc >= threshold:
                pair_indices.append((i, j))
    return {
        "n_pairs": len(pair_indices),
        "pair_indices": pair_indices,
        "rate_per_event": len(pair_indices) / n,
        "n_events": n,
        "threshold": threshold,
    }
