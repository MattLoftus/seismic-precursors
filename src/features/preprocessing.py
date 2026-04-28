"""Standard preprocessing for broadband seismic waveforms.

Pipeline (in order):
    1. Merge gaps with linear interpolation, drop empty traces.
    2. Detrend linear + demean.
    3. Cosine taper (5% on each end).
    4. Remove instrument response to ground velocity (m/s).
    5. Bandpass filter 0.5-20 Hz (covers most local-event energy band).

Returns a new ObsPy Stream; does not mutate the input.
"""
from __future__ import annotations

from obspy import Stream


DEFAULT_FREQMIN = 0.5
DEFAULT_FREQMAX = 20.0
DEFAULT_TAPER = 0.05
ANTIALIAS_FRACTION = 1.0 / 3.0   # cap freqmax at fs/3 to dodge Nyquist aliasing artifacts


def preprocess_waveform(
    stream: Stream,
    inventory=None,
    freqmin: float = DEFAULT_FREQMIN,
    freqmax: float = DEFAULT_FREQMAX,
    taper: float = DEFAULT_TAPER,
    output: str = "VEL",
    log: callable | None = print,
) -> Stream:
    """Standard preprocessing. Returns a new Stream.

    inventory: ObsPy Inventory with the station response. Required for response
    removal (the older `attach_response=True` get_waveforms idiom is deprecated).

    The bandpass freqmax is automatically clamped to min(freqmax, fs/3) per
    trace — Session 3 (exp03) found 60 Hz line noise aliasing exactly to
    Nyquist on BK.PKD's 2004 40 Hz broadband sampler.
    """
    st = stream.copy()
    n_in = len(st)
    st.merge(method=1, fill_value="interpolate")
    st = st.select(channel="?H?")  # broad/high-broadband; drop LH (long-period 1 Hz) and SH (short-period)
    if len(st) == 0:
        if log:
            log(f"[preprocess] WARNING: no traces remain after merge/select (in: {n_in})")
        return st

    st.detrend("linear")
    st.detrend("demean")
    st.taper(max_percentage=taper, type="cosine")

    if inventory is not None:
        try:
            st.remove_response(inventory=inventory, output=output, water_level=60)
            if log:
                log(f"[preprocess] removed response (output={output})")
        except Exception as e:
            if log:
                log(f"[preprocess] WARNING: response removal failed ({type(e).__name__}: {e}); "
                    f"keeping raw counts")
    else:
        # Fall back to the deprecated trace-attached response (still functional)
        if any(getattr(tr.stats, "response", None) is not None for tr in st):
            try:
                st.remove_response(output=output, water_level=60)
                if log:
                    log(f"[preprocess] removed trace-attached response (output={output})")
            except Exception as e:
                if log:
                    log(f"[preprocess] WARNING: trace-attached response removal failed "
                        f"({type(e).__name__}: {e})")
        else:
            if log:
                log("[preprocess] WARNING: no response available; keeping raw counts")

    # Per-trace bandpass with anti-alias cap
    for tr in st:
        fs = float(tr.stats.sampling_rate)
        fmax_eff = min(freqmax, fs * ANTIALIAS_FRACTION)
        if fmax_eff <= freqmin:
            if log:
                log(f"[preprocess] WARNING: {tr.id} fs={fs} too low for {freqmin}-{freqmax} band; "
                    f"applying highpass at {freqmin} Hz only")
            tr.filter("highpass", freq=freqmin, corners=4, zerophase=True)
        else:
            tr.filter("bandpass", freqmin=freqmin, freqmax=fmax_eff, corners=4, zerophase=True)
            if log and abs(fmax_eff - freqmax) > 1e-3:
                log(f"[preprocess] {tr.id} fs={fs}Hz: clamped fmax {freqmax}->{fmax_eff:.2f} (fs/3)")
    if log:
        log(f"[preprocess] processed {len(st)} traces")
    return st
