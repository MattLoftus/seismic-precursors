"""Per-window waveform feature compute.

For each precursor / null window, fetch N short snapshots from the region's
primary broadband station, preprocess, and compute three pre-registered
waveform-derived features per snapshot:

    spectral_slope          Brune log-log slope on Welch PSD
    waveform_entropy        Shannon entropy of whitened spectrogram
    hht_imf1_if_median      EMD IMF1 instantaneous frequency, post-bandpass

The per-window scalars are computed by aggregating across snapshots:
mean for the spectrum-shape / entropy features, median for IMF1 IF
(robust to edge artifacts).

Snapshots are evenly spaced across the 30-day window (default 6 × 10-min).
This is a within-window subsampling strategy NOT specified by the pre-reg;
we document the choice and it's deterministic per (region, t_start, t_end).

Fetch failures are NaN-filled at the snapshot level; per-window aggregation
uses np.nanmean / np.nanmedian. If all snapshots fail, the per-window
features are NaN.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
import pandas as pd
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

from src.features.entropy import windowed_mean_entropy
from src.features.hht import imf1_if_series
from src.features.preprocessing import preprocess_waveform
from src.features.spectral import spectral_slope
from src.regions import Region, Station


@dataclass
class WaveformFeatureParams:
    """Knobs for the waveform-feature pipeline. Defaults mirror exp03 sanity.

    `fdsn_client` here is the FALLBACK if region.fdsn_client is not set;
    Region.fdsn_client takes precedence in `compute_window_features`.
    """
    n_snapshots_per_window: int = 6
    snapshot_seconds: int = 600                # 10-min snapshots
    fdsn_client: str = "IRIS"                  # default if region.fdsn_client missing
    fdsn_timeout: int = 60
    # Restrict channel to broadband Z (skip LHZ 1 Hz long-period and SHZ short-period
    # which would break the bandpass at fs < 24 Hz):
    channel_pattern: str = "BHZ,HHZ,EHZ"
    min_sampling_rate: float = 20.0            # reject traces with fs < this
    location: str = "*"
    bandpass_freqmin: float = 0.5
    bandpass_freqmax: float = 12.0             # post fs/3 anti-alias logic in preprocessing
    emd_max_imfs: int = 4
    entropy_window_seconds: float = 4.0


@dataclass
class SnapshotResult:
    t_center_iso: str
    fs: float
    n_samples: int
    success: bool
    spectral_slope: float = float("nan")
    spectral_r2: float = float("nan")
    entropy_mean_bits: float = float("nan")
    imf1_if_median_hz: float = float("nan")
    error: str = ""


@dataclass
class WindowFeatures:
    region: str
    station: str
    t_start_iso: str
    t_end_iso: str
    n_snapshots: int
    n_successful: int
    spectral_slope: float = float("nan")
    spectral_r2: float = float("nan")
    waveform_entropy: float = float("nan")
    hht_imf1_if_median_hz: float = float("nan")


def snapshot_centers(t_start: UTCDateTime, t_end: UTCDateTime,
                     n: int) -> list[UTCDateTime]:
    """Evenly-spaced snapshot centers across [t_start, t_end)."""
    if n <= 0:
        return []
    duration = float(t_end - t_start)
    # Place snapshots at fractions (k+0.5) / n
    centers = [t_start + duration * (k + 0.5) / n for k in range(n)]
    return centers


def _fetch_one_snapshot(client: Client, station: Station, t_center: UTCDateTime,
                        params: WaveformFeatureParams):
    half = params.snapshot_seconds / 2
    t1 = t_center - half
    t2 = t_center + half
    inv = client.get_stations(
        network=station.network, station=station.station,
        channel=params.channel_pattern,
        starttime=t1, endtime=t2, level="response",
    )
    st = client.get_waveforms(
        network=station.network, station=station.station,
        location=params.location, channel=params.channel_pattern,
        starttime=t1, endtime=t2,
    )
    return st, inv


def compute_window_features(region: Region, t_start_iso: str, t_end_iso: str,
                            params: WaveformFeatureParams | None = None,
                            log: callable | None = None) -> WindowFeatures:
    """Compute the three pre-registered waveform-derived features for one window."""
    if params is None:
        params = WaveformFeatureParams()
    t_start = UTCDateTime(t_start_iso)
    t_end = UTCDateTime(t_end_iso)
    primary = region.primary
    out = WindowFeatures(
        region=region.name, station=primary.code(),
        t_start_iso=t_start_iso, t_end_iso=t_end_iso,
        n_snapshots=params.n_snapshots_per_window, n_successful=0,
    )

    fdsn_short = getattr(region, "fdsn_client", None) or params.fdsn_client
    client = Client(fdsn_short, timeout=params.fdsn_timeout)
    centers = snapshot_centers(t_start, t_end, params.n_snapshots_per_window)

    snapshot_results: list[SnapshotResult] = []
    for tc in centers:
        snap = SnapshotResult(t_center_iso=tc.isoformat(), fs=float("nan"),
                              n_samples=0, success=False)
        try:
            st, inv = _fetch_one_snapshot(client, primary, tc, params)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                st_proc = preprocess_waveform(
                    st, inventory=inv,
                    freqmin=params.bandpass_freqmin,
                    freqmax=params.bandpass_freqmax,
                    log=None,
                )
            z = st_proc.select(channel="?HZ")
            if not z:
                snap.error = "no Z channel after preprocessing"
                snapshot_results.append(snap)
                continue
            tr = z[0]
            data = tr.data.astype(float)
            fs = float(tr.stats.sampling_rate)
            if fs < params.min_sampling_rate:
                snap.error = (f"sampling rate {fs} Hz below minimum "
                              f"{params.min_sampling_rate} (rejecting LHZ/SHZ-style traces)")
                snapshot_results.append(snap)
                continue
            if len(data) < int(0.5 * params.snapshot_seconds * fs):
                snap.error = f"insufficient samples: n={len(data)} for fs={fs}"
                snapshot_results.append(snap)
                continue
            snap.fs = fs
            snap.n_samples = int(len(data))

            # Spectral slope
            sp = spectral_slope(data, fs=fs, fmin=0.5, fmax=params.bandpass_freqmax)
            snap.spectral_slope = float(sp.get("slope", float("nan")))
            snap.spectral_r2 = float(sp.get("r_squared", float("nan")))

            # Entropy
            ent_mean, _ = windowed_mean_entropy(
                data, fs=fs,
                window_seconds=params.entropy_window_seconds,
                fmin=0.5, fmax=params.bandpass_freqmax,
            )
            snap.entropy_mean_bits = float(ent_mean)

            # HHT IMF1 IF. Bandpass freqmax must be < fs/2 (Nyquist); we clamp
            # to fs/3 to dodge any anti-alias-filter ringing near Nyquist.
            fmax_emd = min(params.bandpass_freqmax, fs / 3.0)
            if fmax_emd <= params.bandpass_freqmin:
                snap.error = f"fs={fs} too low for bandpass {params.bandpass_freqmin}-{fmax_emd}"
                snapshot_results.append(snap)
                continue
            from scipy.signal import butter, sosfiltfilt
            sos = butter(4, [params.bandpass_freqmin, fmax_emd],
                         btype="bandpass", fs=fs, output="sos")
            data_bp = sosfiltfilt(sos, data)
            try:
                _, if1 = imf1_if_series(data_bp, fs=fs,
                                         max_imfs=params.emd_max_imfs)
                if len(if1) > 0:
                    snap.imf1_if_median_hz = float(np.median(if1))
            except Exception as e:
                snap.error = f"hht: {type(e).__name__}: {e}"

            snap.success = True
        except Exception as e:
            snap.error = f"fetch: {type(e).__name__}: {e}"
            if log:
                log(f"[wf] snapshot {tc.isoformat()} failed: {snap.error}")
        snapshot_results.append(snap)

    succ = [s for s in snapshot_results if s.success]
    out.n_successful = len(succ)
    if succ:
        out.spectral_slope = float(np.nanmean([s.spectral_slope for s in succ]))
        out.spectral_r2 = float(np.nanmean([s.spectral_r2 for s in succ]))
        out.waveform_entropy = float(np.nanmean([s.entropy_mean_bits for s in succ]))
        out.hht_imf1_if_median_hz = float(np.nanmedian(
            [s.imf1_if_median_hz for s in succ if np.isfinite(s.imf1_if_median_hz)]
        )) if any(np.isfinite(s.imf1_if_median_hz) for s in succ) else float("nan")
    return out


def compute_features_for_dataframe(
    region: Region,
    df: pd.DataFrame,
    *,
    params: WaveformFeatureParams | None = None,
    progress_every: int = 5,
    log: callable | None = print,
    n_workers: int = 1,
) -> pd.DataFrame:
    """Apply `compute_window_features` to each row of df.

    df must contain `t_start` and `t_end` columns (ISO datetime). Returns a new
    DataFrame with original columns + waveform feature columns.

    `n_workers > 1` parallelizes via ThreadPoolExecutor. Each thread creates
    its own FDSN client (cheap; thread-safe under ObsPy). 4 workers gives
    roughly 4x speedup on FDSN-bound workloads; beyond ~6 workers the rate
    limits dominate.
    """
    if n_workers <= 1:
        return _compute_serial(region, df, params=params,
                               progress_every=progress_every, log=log)

    from concurrent.futures import ThreadPoolExecutor, as_completed

    items = []
    for i, row in enumerate(df.itertuples(index=False)):
        t_start_iso = pd.Timestamp(row.t_start).isoformat()
        t_end_iso = pd.Timestamp(row.t_end).isoformat()
        rec = dict(row._asdict()) if hasattr(row, "_asdict") else df.iloc[i].to_dict()
        items.append((i, t_start_iso, t_end_iso, rec))

    results: dict[int, dict] = {}

    def worker(item):
        i, t_start_iso, t_end_iso, rec = item
        wf = compute_window_features(region, t_start_iso, t_end_iso,
                                     params=params, log=None)
        rec = dict(rec)
        rec.update({
            "wf_station": wf.station,
            "wf_n_snapshots": wf.n_snapshots,
            "wf_n_successful": wf.n_successful,
            "wf_spectral_slope": wf.spectral_slope,
            "wf_spectral_r2": wf.spectral_r2,
            "wf_waveform_entropy": wf.waveform_entropy,
            "wf_hht_imf1_if_median_hz": wf.hht_imf1_if_median_hz,
        })
        return i, rec, wf

    completed = 0
    with ThreadPoolExecutor(max_workers=n_workers) as ex:
        futures = [ex.submit(worker, it) for it in items]
        for fut in as_completed(futures):
            i, rec, wf = fut.result()
            results[i] = rec
            completed += 1
            if log and completed % progress_every == 0:
                log(f"[wf] {completed}/{len(df)}  "
                    f"({100*completed/len(df):.0f}%)  "
                    f"last: succ={wf.n_successful}/{wf.n_snapshots} "
                    f"slope={wf.spectral_slope:.3f} ent={wf.waveform_entropy:.2f} "
                    f"if1={wf.hht_imf1_if_median_hz:.2f}Hz")

    rows = [results[i] for i in range(len(items))]
    return pd.DataFrame(rows)


def _compute_serial(region: Region, df: pd.DataFrame, *,
                    params: WaveformFeatureParams | None,
                    progress_every: int, log) -> pd.DataFrame:
    rows = []
    for i, row in enumerate(df.itertuples(index=False)):
        t_start_iso = pd.Timestamp(row.t_start).isoformat()
        t_end_iso = pd.Timestamp(row.t_end).isoformat()
        wf = compute_window_features(region, t_start_iso, t_end_iso,
                                     params=params, log=log)
        rec = dict(row._asdict()) if hasattr(row, "_asdict") else df.iloc[i].to_dict()
        rec.update({
            "wf_station": wf.station,
            "wf_n_snapshots": wf.n_snapshots,
            "wf_n_successful": wf.n_successful,
            "wf_spectral_slope": wf.spectral_slope,
            "wf_spectral_r2": wf.spectral_r2,
            "wf_waveform_entropy": wf.waveform_entropy,
            "wf_hht_imf1_if_median_hz": wf.hht_imf1_if_median_hz,
        })
        rows.append(rec)
        if log and (i + 1) % progress_every == 0:
            log(f"[wf] {i+1}/{len(df)}  ({100*(i+1)/len(df):.0f}%)  "
                f"last: succ={wf.n_successful}/{wf.n_snapshots} "
                f"slope={wf.spectral_slope:.3f} ent={wf.waveform_entropy:.2f} "
                f"if1={wf.hht_imf1_if_median_hz:.2f}Hz")
    return pd.DataFrame(rows)
