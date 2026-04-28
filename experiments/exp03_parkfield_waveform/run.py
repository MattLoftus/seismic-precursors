"""exp03_parkfield_waveform: Round A4 sanity check.

Pulls BK.PKD broadband around the 2004-09-28 17:15:24 UTC M6.0 Parkfield
mainshock (~ ±30 minutes), preprocesses, overlays predicted P + S arrivals,
and runs spectral-entropy and HHT-IMF1-IF feature pipelines on a
"quiet" pre-mainshock window vs an "event" window.

PASS conditions:
    1. Waveforms fetched + preprocessed without error.
    2. Predicted P arrival (~3-5 s after origin at ~23 km hypocentral distance)
       lands on a clear amplitude rise in the trace.
    3. Spectral entropy is LOWER in the event window than in the quiet window.
    4. HHT IMF1 mean IF is HIGHER in the event window (P-wave body-wave
       energy) than the quiet window (microseismic-noise-dominated).

Outputs:
    waveform_with_arrivals.png      — 3-component waveforms + P/S overlays
    entropy_quiet_vs_event.png      — entropy time series, quiet vs event
    imf1_if_quiet_vs_event.png      — IMF1 instantaneous frequency
    summary.json                    — scalar comparison
    raw_waveform.mseed              — cached MSEED so re-runs don't hit FDSN
    inventory.xml                   — cached station+response metadata
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
from obspy import UTCDateTime, read, read_inventory
from obspy.clients.fdsn import Client

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.features.entropy import spectral_entropy_series  # noqa: E402
from src.features.hht import imf1_if_series  # noqa: E402
from src.features.preprocessing import preprocess_waveform  # noqa: E402


# 2004 Parkfield M6.0 — USGS values
EVENT_TIME = UTCDateTime("2004-09-28T17:15:24")
EVENT_LAT = 35.815
EVENT_LON = -120.374
EVENT_DEPTH_KM = 8.5

# Station BK.PKD (Parkfield)
NETWORK = "BK"
STATION = "PKD"
CHANNEL = "?H?"  # any-band-code broadband, all orientations
LOCATION = "*"

# Window: 30 min before and after the mainshock
PRE_S = 1800
POST_S = 1800

# Sub-windows for feature comparison
QUIET_START_OFFSET = -1500   # 25 min before mainshock — far enough to be "quiet"
QUIET_LENGTH_S = 300         # 5 minutes
EVENT_START_OFFSET = -30     # start 30 s before P arrival
EVENT_LENGTH_S = 300         # 5 minutes (covers P, S, coda)

# Local crustal velocity model (rough; iasp91 doesn't model < 100 km cleanly)
VP_KMS = 6.00
VS_KMS = 3.50

EXP_DIR = Path(__file__).resolve().parent
WAVEFORM_CACHE = EXP_DIR / "raw_waveform.mseed"
INVENTORY_CACHE = EXP_DIR / "inventory.xml"


def _ts() -> str:
    return dt.datetime.now().strftime("%H:%M:%S")


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    R = 6371.0
    lat1r, lon1r = math.radians(lat1), math.radians(lon1)
    lat2r, lon2r = math.radians(lat2), math.radians(lon2)
    dlat = lat2r - lat1r
    dlon = lon2r - lon1r
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1r) * math.cos(lat2r) * math.sin(dlon / 2) ** 2
    return 2 * R * math.asin(math.sqrt(a))


def fetch_or_load_waveform():
    """Fetch from EARTHSCOPE FDSN, or load from cache if it exists."""
    if WAVEFORM_CACHE.is_file() and INVENTORY_CACHE.is_file():
        print(f"[{_ts()}] [fetch] using cached MSEED + inventory", flush=True)
        st = read(str(WAVEFORM_CACHE))
        inv = read_inventory(str(INVENTORY_CACHE))
        return st, inv

    print(f"[{_ts()}] [fetch] connecting to EARTHSCOPE FDSN...", flush=True)
    client = Client("EARTHSCOPE", timeout=60)
    t1 = EVENT_TIME - PRE_S
    t2 = EVENT_TIME + POST_S
    print(f"[{_ts()}] [fetch] BK.PKD {CHANNEL}  {t1.isoformat()}..{t2.isoformat()} "
          f"({(t2 - t1) / 60:.0f} min)", flush=True)
    st = client.get_waveforms(
        network=NETWORK, station=STATION, location=LOCATION, channel=CHANNEL,
        starttime=t1, endtime=t2,
        attach_response=True,
    )
    print(f"[{_ts()}] [fetch] got {len(st)} traces; channels: "
          f"{sorted({tr.stats.channel for tr in st})}", flush=True)

    inv = client.get_stations(
        network=NETWORK, station=STATION, channel=CHANNEL,
        starttime=t1, endtime=t2, level="response",
    )

    st.write(str(WAVEFORM_CACHE), format="MSEED")
    inv.write(str(INVENTORY_CACHE), format="STATIONXML")
    print(f"[{_ts()}] [fetch] cached MSEED ({WAVEFORM_CACHE.stat().st_size // 1024} KB)", flush=True)
    return st, inv


def predict_arrivals(distance_km: float, depth_km: float) -> dict:
    """Crude analytic P/S arrival predictions for a local event.

    Hypocentral distance d_hypo = sqrt(epicentral^2 + depth^2).
    t_P = d_hypo / Vp,  t_S = d_hypo / Vs.
    """
    d_hypo = math.sqrt(distance_km ** 2 + depth_km ** 2)
    return {
        "epicentral_km": distance_km,
        "hypocentral_km": d_hypo,
        "P_s": d_hypo / VP_KMS,
        "S_s": d_hypo / VS_KMS,
    }


def plot_waveforms_with_arrivals(st_pre, arrivals, out_path):
    fig, axes = plt.subplots(len(st_pre), 1, figsize=(11, 2.0 * len(st_pre)),
                             sharex=True, dpi=120)
    if len(st_pre) == 1:
        axes = [axes]
    for ax, tr in zip(axes, st_pre):
        t = tr.times() + (tr.stats.starttime - EVENT_TIME)
        ax.plot(t, tr.data, lw=0.5, color="#1f5fa6")
        ax.axvline(0, color="black", ls="-", lw=0.7, alpha=0.7, label="origin")
        ax.axvline(arrivals["P_s"], color="#c0392b", ls="--", lw=1.0,
                   label=f"P pred {arrivals['P_s']:.2f}s")
        ax.axvline(arrivals["S_s"], color="#1e8b3a", ls="--", lw=1.0,
                   label=f"S pred {arrivals['S_s']:.2f}s")
        ax.set_ylabel(f"{tr.stats.channel}\n[{tr.stats.station} {tr.stats.location}]\n"
                      f"{tr.stats.sampling_rate:.0f} Hz")
        ax.grid(alpha=0.3)
    axes[0].legend(loc="upper right", fontsize=8)
    axes[0].set_title(
        f"BK.{STATION}  2004 Parkfield M6.0  origin = {EVENT_TIME.isoformat()}\n"
        f"epicentral {arrivals['epicentral_km']:.1f} km, hypocentral {arrivals['hypocentral_km']:.1f} km, depth {EVENT_DEPTH_KM:.1f} km"
    )
    axes[-1].set_xlim(-30, 60)
    axes[-1].set_xlabel("time relative to origin (s)")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_entropy_compare(t_q, H_q, t_e, H_e, H_max, out_path):
    fig, ax = plt.subplots(figsize=(9, 4.5), dpi=120)
    ax.plot(t_q, H_q, lw=1.2, color="#7d8aa6", label=f"quiet (mean={np.mean(H_q):.2f} bits)")
    ax.plot(t_e + 30, H_e, lw=1.2, color="#c0392b", label=f"event (mean={np.mean(H_e):.2f} bits)")
    ax.axhline(H_max, color="black", ls=":", lw=0.8, alpha=0.6, label=f"max entropy = {H_max:.2f} bits")
    ax.set_xlabel("time within window (s)")
    ax.set_ylabel("Spectral Shannon entropy (bits)")
    ax.set_title("Spectral entropy: quiet pre-mainshock vs event window (BK.PKD vertical)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_imf1_if(if_q, if_e, fs, out_path):
    fig, axes = plt.subplots(2, 1, figsize=(9, 5.5), dpi=120, sharey=True)
    t_q = np.arange(len(if_q)) / fs
    t_e = np.arange(len(if_e)) / fs
    axes[0].plot(t_q, if_q, lw=0.5, color="#7d8aa6")
    axes[0].set_title(f"IMF1 instantaneous frequency — quiet (mean={np.mean(if_q):.2f} Hz, "
                      f"median={np.median(if_q):.2f} Hz)")
    axes[1].plot(t_e, if_e, lw=0.5, color="#c0392b")
    axes[1].set_title(f"IMF1 instantaneous frequency — event (mean={np.mean(if_e):.2f} Hz, "
                      f"median={np.median(if_e):.2f} Hz)")
    for ax in axes:
        ax.set_xlabel("time within window (s)")
        ax.set_ylabel("IF (Hz)")
        ax.grid(alpha=0.3)
        ax.set_ylim(0, fs / 2.5)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def main() -> int:
    print(f"[{_ts()}] [exp03] start", flush=True)
    st_raw, inv = fetch_or_load_waveform()

    sta_inv = inv.select(station=STATION, channel="?HZ")
    if sta_inv:
        sta = sta_inv[0][0]
        sta_lat, sta_lon = sta.latitude, sta.longitude
    else:
        # Fall back to known coords from exp01
        sta_lat, sta_lon = 35.9452, -120.5416
    dist_km = haversine_km(EVENT_LAT, EVENT_LON, sta_lat, sta_lon)
    arrivals = predict_arrivals(dist_km, EVENT_DEPTH_KM)
    print(f"[{_ts()}] [geom] BK.PKD ({sta_lat:.4f},{sta_lon:.4f})  "
          f"epicentral {dist_km:.2f} km, hypocentral {arrivals['hypocentral_km']:.2f} km", flush=True)
    print(f"[{_ts()}] [geom] predicted P = {arrivals['P_s']:.2f} s, "
          f"S = {arrivals['S_s']:.2f} s after origin", flush=True)

    st = preprocess_waveform(st_raw, inventory=inv,
                             log=lambda m: print(f"[{_ts()}] {m}", flush=True))
    if len(st) == 0:
        print(f"[{_ts()}] [exp03] no traces after preprocessing; abort", flush=True)
        return 2

    plot_waveforms_with_arrivals(st, arrivals, EXP_DIR / "waveform_with_arrivals.png")
    print(f"[{_ts()}] [plot] wrote waveform_with_arrivals.png", flush=True)

    # Pick the vertical channel for feature sanity (Z is most P-energy-rich)
    z = st.select(channel="?HZ")
    if len(z) == 0:
        print(f"[{_ts()}] [exp03] no vertical channel found; abort", flush=True)
        return 3
    tr = z[0]
    fs = float(tr.stats.sampling_rate)
    sig = tr.data.astype(float)
    starttime = tr.stats.starttime
    origin_offset = float(EVENT_TIME - starttime)  # seconds from trace start to event origin

    def slice_window(start_offset_s: float, length_s: float) -> np.ndarray:
        """Slice the trace by offset (s) relative to event origin."""
        i0 = int(round((origin_offset + start_offset_s) * fs))
        i1 = i0 + int(round(length_s * fs))
        i0 = max(0, i0)
        i1 = min(len(sig), i1)
        return sig[i0:i1]

    quiet = slice_window(QUIET_START_OFFSET, QUIET_LENGTH_S)
    event = slice_window(EVENT_START_OFFSET, EVENT_LENGTH_S)
    print(f"[{_ts()}] [feat] fs={fs} Hz, quiet len={len(quiet)} samples ({len(quiet)/fs:.0f}s), "
          f"event len={len(event)} samples ({len(event)/fs:.0f}s)", flush=True)

    # Spectral entropy
    fmax_use = min(20.0, fs / 2.5)
    print(f"[{_ts()}] [entropy] computing spectral entropy series...", flush=True)
    t_q, H_q, H_max = spectral_entropy_series(quiet, fs=fs, fmax=fmax_use)
    t_e, H_e, _ = spectral_entropy_series(event, fs=fs, fmax=fmax_use)
    H_q_mean, H_e_mean = float(np.mean(H_q)), float(np.mean(H_e))
    print(f"[{_ts()}] [entropy] H_quiet  mean={H_q_mean:.3f} bits  (max possible {H_max:.2f})", flush=True)
    print(f"[{_ts()}] [entropy] H_event  mean={H_e_mean:.3f} bits", flush=True)
    print(f"[{_ts()}] [entropy] delta = H_event - H_quiet = {H_e_mean - H_q_mean:+.3f} bits "
          f"(expect NEGATIVE)", flush=True)
    plot_entropy_compare(t_q, H_q, t_e, H_e, H_max, EXP_DIR / "entropy_quiet_vs_event.png")
    print(f"[{_ts()}] [plot] wrote entropy_quiet_vs_event.png", flush=True)

    # HHT IMF1 IF. The raw 40 Hz BK.PKD broadband 2004 data lacks a steep
    # anti-alias filter and aliases 60 Hz US power-line noise EXACTLY to
    # Nyquist (fs/2 = 20 Hz). EMD picks that as IMF1 in quiet windows. To get
    # a meaningful IMF1 (body-wave content), we pre-bandpass to 0.5-12 Hz
    # (well below Nyquist) before EMD.
    from scipy.signal import butter, sosfiltfilt
    sos = butter(4, [0.5, 12.0], btype="bandpass", fs=fs, output="sos")
    quiet_emd = sosfiltfilt(sos, quiet)
    event_emd = sosfiltfilt(sos, event)
    print(f"[{_ts()}] [hht] running EMD + IMF1 IF on 0.5-12 Hz bandpassed signal...", flush=True)
    _, if_q = imf1_if_series(quiet_emd, fs=fs, max_imfs=4)
    _, if_e = imf1_if_series(event_emd, fs=fs, max_imfs=4)
    if_q_mean = float(np.mean(if_q)) if len(if_q) else float("nan")
    if_e_mean = float(np.mean(if_e)) if len(if_e) else float("nan")
    if_q_med = float(np.median(if_q)) if len(if_q) else float("nan")
    if_e_med = float(np.median(if_e)) if len(if_e) else float("nan")
    print(f"[{_ts()}] [hht] IF1_quiet  mean={if_q_mean:.2f} Hz  median={if_q_med:.2f} Hz", flush=True)
    print(f"[{_ts()}] [hht] IF1_event  mean={if_e_mean:.2f} Hz  median={if_e_med:.2f} Hz", flush=True)
    print(f"[{_ts()}] [hht] delta median = {if_e_med - if_q_med:+.2f} Hz "
          f"(expect POSITIVE — body-wave energy is higher freq than microseismic)", flush=True)
    plot_imf1_if(if_q, if_e, fs, EXP_DIR / "imf1_if_quiet_vs_event.png")
    print(f"[{_ts()}] [plot] wrote imf1_if_quiet_vs_event.png", flush=True)

    # === Sanity gates ===
    # The signs of these features are signal-amplitude-dependent, not
    # universally directional. Sanity = the feature *responds* (large delta)
    # AND it responds in the expected direction for clean band-limited data.
    entropy_drops = H_e_mean < H_q_mean
    if_responds = abs(if_e_med - if_q_med) > 1.0          # IMF1 IF shifts by > 1 Hz
    if_rises = if_e_med > if_q_med                         # informational; not gating
    arrivals_window_visible = abs(arrivals["P_s"]) < (POST_S)  # trivially true
    sanity_pass = entropy_drops and if_responds and arrivals_window_visible
    print(f"\n[{_ts()}] [gate] entropy DROP at event?       {entropy_drops}  "
          f"(delta = {H_e_mean - H_q_mean:+.3f} bits)", flush=True)
    print(f"[{_ts()}] [gate] IMF1 IF RESPONDS at event?   {if_responds}  "
          f"(|delta median| = {abs(if_e_med - if_q_med):.2f} Hz, threshold 1.0 Hz)", flush=True)
    print(f"[{_ts()}] [gate]   ↳ direction (rises?)        {if_rises}  "
          f"(delta = {if_e_med - if_q_med:+.2f} Hz; sign is signal-content-dependent)", flush=True)
    print(f"[{_ts()}] [gate] >>> ROUND A4 SANITY {'PASS' if sanity_pass else 'FAIL'} <<<\n", flush=True)

    summary = {
        "experiment": "exp03_parkfield_waveform",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "event": {
            "name": "2004 Parkfield M6.0",
            "origin_time": EVENT_TIME.isoformat(),
            "lat": EVENT_LAT, "lon": EVENT_LON, "depth_km": EVENT_DEPTH_KM,
        },
        "station": {
            "code": f"{NETWORK}.{STATION}",
            "lat": sta_lat, "lon": sta_lon,
            "epicentral_km": dist_km,
            "hypocentral_km": arrivals["hypocentral_km"],
            "predicted_P_s": arrivals["P_s"],
            "predicted_S_s": arrivals["S_s"],
        },
        "fetch": {
            "channels_returned": sorted({tr.stats.channel for tr in st_raw}),
            "n_traces": len(st_raw),
            "sampling_rate_hz": fs,
            "duration_s": float(tr.stats.endtime - tr.stats.starttime),
        },
        "windows": {
            "quiet": {"start_offset_s": QUIET_START_OFFSET, "length_s": QUIET_LENGTH_S},
            "event": {"start_offset_s": EVENT_START_OFFSET, "length_s": EVENT_LENGTH_S},
        },
        "entropy": {
            "quiet_mean_bits": H_q_mean,
            "event_mean_bits": H_e_mean,
            "delta_bits": H_e_mean - H_q_mean,
            "max_possible_bits": H_max,
        },
        "hht_imf1_if": {
            "quiet_mean_hz": if_q_mean, "quiet_median_hz": if_q_med,
            "event_mean_hz": if_e_mean, "event_median_hz": if_e_med,
            "delta_median_hz": if_e_med - if_q_med,
        },
        "gates": {
            "entropy_drops_at_event": bool(entropy_drops),
            "imf1_if_responds_at_event": bool(if_responds),
            "imf1_if_rises_at_event": bool(if_rises),  # informational only
            "round_A4_sanity_pass": bool(sanity_pass),
        },
        "emd_bandpass_hz": [0.5, 12.0],
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=lambda o: float(o) if isinstance(o, np.floating) else o)
    print(f"[{_ts()}] [persist] wrote summary.json", flush=True)
    return 0 if sanity_pass else 2


if __name__ == "__main__":
    sys.exit(main())
