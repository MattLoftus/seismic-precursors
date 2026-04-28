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

    inventory: ObsPy Inventory with the station response. If None and traces have
    response attached (`attach_response=True` on get_waveforms), uses that. If neither
    is available, the response-removal step is skipped and a warning is logged.
    """
    st = stream.copy()
    n_in = len(st)
    st.merge(method=1, fill_value="interpolate")
    st = st.select(channel="BH?")  # broadband only; drop short-period if any
    n_after_merge = len(st)
    if n_after_merge == 0:
        if log:
            log(f"[preprocess] WARNING: no traces remain after merge/select (in: {n_in})")
        return st

    st.detrend("linear")
    st.detrend("demean")
    st.taper(max_percentage=taper, type="cosine")

    have_response = inventory is not None or any(
        getattr(tr.stats, "response", None) is not None for tr in st
    )
    if have_response:
        try:
            st.remove_response(inventory=inventory, output=output, water_level=60)
            if log:
                log(f"[preprocess] removed response (output={output})")
        except Exception as e:
            if log:
                log(f"[preprocess] WARNING: response removal failed ({type(e).__name__}: {e}); "
                    f"keeping raw counts")
    else:
        if log:
            log("[preprocess] WARNING: no response available; keeping raw counts")

    st.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)
    if log:
        log(f"[preprocess] bandpass {freqmin:.2f}-{freqmax:.1f} Hz, {len(st)} traces")
    return st
