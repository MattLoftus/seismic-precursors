"""exp01_parkfield_bvalue: calibration gate for the seismic-precursors pipeline.

Two checks executed back-to-back:

(1) IRIS hello-world: fetch BK.PKD station metadata via FDSN. Confirms the seismic
    waveform stack is alive end-to-end (network -> ObsPy -> inventory parse).

(2) Parkfield b-value: pull ANSS ComCat catalog for the Parkfield region
    (50 km radius around 35.85 N / 120.40 W, M >= 1.0, 2000-2024), estimate
    Mc via Wiemer & Wyss 2000 max-curvature + Woessner & Wiemer 2005 +0.2
    correction, then compute Aki 1965 MLE b-value with bootstrap CI.

PASS GATE: |b - 0.9| / 0.9 <= 0.10  (Bakun et al. 2005 reports b ~ 0.85-0.95
for the Parkfield earthquake-prediction segment of the San Andreas).

Outputs: magnitude_frequency.png, summary.json, console log.
"""
from __future__ import annotations

import datetime as dt
import json
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Make src/ importable when run as a script.
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.features.bvalue import (  # noqa: E402
    bvalue_with_bootstrap,
    magnitude_of_completeness_max_curvature,
)


PARKFIELD_LAT = 35.85
PARKFIELD_LON = -120.40
PARKFIELD_RADIUS_KM = 50.0
M_MIN = 1.0
START = dt.datetime(2000, 1, 1)
END = dt.datetime(2025, 1, 1)
CHUNK_YEARS = 2

BAKUN_2005_B = 0.9
B_TOLERANCE = 0.10  # +/- 10%

EXP_DIR = Path(__file__).resolve().parent


def hello_world_iris() -> dict:
    """Fetch BK.PKD station metadata. Returns inventory summary as dict."""
    print(f"[{_ts()}] [hello-world] connecting to IRIS FDSN...", flush=True)
    from obspy.clients.fdsn import Client

    client = Client("IRIS", timeout=30)
    inv = client.get_stations(network="BK", station="PKD", level="station")
    nets = [n.code for n in inv]
    stas = [(n.code, s.code, s.latitude, s.longitude, s.elevation, str(s.start_date)[:10]) for n in inv for s in n]
    print(f"[{_ts()}] [hello-world] OK -- networks={nets} stations={[s[1] for s in stas]}", flush=True)
    for code_n, code_s, lat, lon, elev, start in stas:
        print(f"               {code_n}.{code_s}  {lat:.4f},{lon:.4f}  elev={elev:.0f}m  start={start}", flush=True)
    return {"networks": nets, "stations": stas}


def fetch_parkfield_catalog() -> np.ndarray:
    """Pull Parkfield M>=M_MIN events 2000-2024 in chunks. Returns array of magnitudes."""
    from libcomcat.search import search

    print(
        f"[{_ts()}] [catalog] Parkfield M>={M_MIN}, "
        f"radius={PARKFIELD_RADIUS_KM}km @ ({PARKFIELD_LAT}, {PARKFIELD_LON}), "
        f"{START.date()} to {END.date()}",
        flush=True,
    )
    all_mags: list[float] = []
    chunk_count = 0
    n_chunks = (END.year - START.year) // CHUNK_YEARS
    cur = START
    while cur < END:
        nxt = dt.datetime(min(cur.year + CHUNK_YEARS, END.year), 1, 1)
        chunk_count += 1
        for attempt in (1, 2, 3):
            try:
                events = search(
                    starttime=cur,
                    endtime=nxt,
                    minmagnitude=M_MIN,
                    latitude=PARKFIELD_LAT,
                    longitude=PARKFIELD_LON,
                    maxradiuskm=PARKFIELD_RADIUS_KM,
                )
                break
            except Exception as e:  # libcomcat raises ConnectionError, JSONDecodeError, etc.
                if attempt == 3:
                    raise
                wait = 5 * attempt
                print(f"[{_ts()}] [catalog]   chunk {chunk_count}/{n_chunks} {cur.year}-{nxt.year} "
                      f"attempt {attempt}/3 failed ({type(e).__name__}); retrying in {wait}s", flush=True)
                time.sleep(wait)
        mags = [float(ev.magnitude) for ev in events if ev.magnitude is not None]
        all_mags.extend(mags)
        print(
            f"[{_ts()}] [catalog]   chunk {chunk_count}/{n_chunks} "
            f"{cur.year}-{nxt.year}: {len(events):>5d} events  "
            f"(running total: {len(all_mags)})",
            flush=True,
        )
        cur = nxt
    arr = np.array(all_mags, dtype=float)
    print(f"[{_ts()}] [catalog] done -- {len(arr)} events, M range "
          f"[{arr.min():.2f}, {arr.max():.2f}]", flush=True)
    return arr


def plot_fmd(magnitudes: np.ndarray, mc: float, b: float, a: float, out_path: Path) -> None:
    """Frequency-magnitude distribution: cumulative log10 N vs M, with G-R fit overlay."""
    dm = 0.1
    edges = np.arange(np.floor(magnitudes.min() * 10) / 10, np.ceil(magnitudes.max() * 10) / 10 + dm, dm)
    counts, _ = np.histogram(magnitudes, bins=edges)
    centers = edges[:-1] + dm / 2
    cum_counts = np.cumsum(counts[::-1])[::-1]

    fig, ax = plt.subplots(figsize=(7.0, 4.5), dpi=120)
    ax.semilogy(centers, np.where(counts > 0, counts, np.nan), "o", ms=3.5, color="#999",
                label="binned N(M)")
    ax.semilogy(centers, np.where(cum_counts > 0, cum_counts, np.nan), "s", ms=3.5, color="#1f5fa6",
                label="cumulative N(M >= m)")
    m_fit = np.linspace(mc, magnitudes.max() + 0.2, 50)
    ax.semilogy(m_fit, 10 ** (a - b * m_fit), "-", lw=1.5, color="#c0392b",
                label=f"G-R fit:  log10 N = {a:.2f} - {b:.3f} M")
    ax.axvline(mc, color="k", ls="--", lw=0.8, alpha=0.6, label=f"Mc = {mc:.2f}")
    ax.set_xlabel("Magnitude")
    ax.set_ylabel("N")
    ax.set_title(f"Parkfield FMD ({START.year}-{END.year - 1}, r={PARKFIELD_RADIUS_KM:.0f} km, "
                 f"M>={M_MIN})  N={len(magnitudes)}")
    ax.legend(loc="lower left", fontsize=8, framealpha=0.95)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def _ts() -> str:
    return dt.datetime.now().strftime("%H:%M:%S")


def main() -> int:
    EXP_DIR.mkdir(parents=True, exist_ok=True)

    # === Hello world ===
    inv_summary = hello_world_iris()

    # === Catalog + b-value ===
    mags = fetch_parkfield_catalog()
    np.savez_compressed(EXP_DIR / "catalog_magnitudes.npz", magnitudes=mags)
    print(f"[{_ts()}] [persist] saved catalog_magnitudes.npz ({mags.nbytes // 1024} KB)", flush=True)

    mc_raw = magnitude_of_completeness_max_curvature(mags)
    print(f"[{_ts()}] [bvalue] Mc raw (max-curvature) = {mc_raw:.2f}", flush=True)

    res = bvalue_with_bootstrap(mags, mc=None, n_boot=1000, mc_correction=0.2)
    print(f"[{_ts()}] [bvalue] Mc adjusted = {res.mc:.2f}  "
          f"({res.method_mc}, n>=Mc = {res.n_above_mc}/{res.n_total})", flush=True)
    print(f"[{_ts()}] [bvalue] b = {res.b:.4f} +/- {res.b_se:.4f} (Shi-Bolt 1982)", flush=True)
    print(f"[{_ts()}] [bvalue] b 95% bootstrap CI = [{res.b_boot_lo:.4f}, {res.b_boot_hi:.4f}]", flush=True)
    print(f"[{_ts()}] [bvalue] a (log10 N at M=0) = {res.a:.4f}", flush=True)

    plot_path = EXP_DIR / "magnitude_frequency.png"
    plot_fmd(mags, mc=res.mc, b=res.b, a=res.a, out_path=plot_path)
    print(f"[{_ts()}] [plot] wrote {plot_path.name}", flush=True)

    # === Calibration gate ===
    rel_err = abs(res.b - BAKUN_2005_B) / BAKUN_2005_B
    bakun_in_ci = res.b_boot_lo <= BAKUN_2005_B <= res.b_boot_hi
    gate_passed = rel_err <= B_TOLERANCE or bakun_in_ci
    verdict = "PASS" if gate_passed else "FAIL"

    print(f"\n[{_ts()}] [gate] Bakun+ 2005 reference b = {BAKUN_2005_B:.2f}", flush=True)
    print(f"[{_ts()}] [gate] |our_b - 0.9| / 0.9 = {rel_err:.3%}  (tolerance {B_TOLERANCE:.0%})", flush=True)
    print(f"[{_ts()}] [gate] Bakun b inside our 95% CI? {bakun_in_ci}", flush=True)
    print(f"[{_ts()}] [gate] >>> CALIBRATION {verdict} <<<", flush=True)

    summary = {
        "experiment": "exp01_parkfield_bvalue",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "region": "Parkfield",
        "params": {
            "lat": PARKFIELD_LAT, "lon": PARKFIELD_LON,
            "radius_km": PARKFIELD_RADIUS_KM,
            "M_min": M_MIN,
            "start": START.isoformat(),
            "end": END.isoformat(),
            "chunk_years": CHUNK_YEARS,
        },
        "iris_hello_world": {
            "networks": inv_summary["networks"],
            "n_stations": len(inv_summary["stations"]),
        },
        "catalog": {
            "n_events": int(len(mags)),
            "m_min_observed": float(mags.min()) if len(mags) else None,
            "m_max_observed": float(mags.max()) if len(mags) else None,
        },
        "bvalue": {
            "Mc_raw_max_curvature": float(mc_raw),
            "Mc_used": float(res.mc),
            "Mc_method": res.method_mc,
            "n_above_Mc": int(res.n_above_mc),
            "b": float(res.b),
            "b_se_shi_bolt_1982": float(res.b_se),
            "b_boot_ci95": [float(res.b_boot_lo), float(res.b_boot_hi)],
            "a_intercept": float(res.a),
        },
        "calibration_gate": {
            "reference_paper": "Bakun et al. 2005",
            "reference_b": BAKUN_2005_B,
            "tolerance_frac": B_TOLERANCE,
            "rel_error": float(rel_err),
            "reference_in_bootstrap_ci": bool(bakun_in_ci),
            "passed": bool(gate_passed),
        },
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[{_ts()}] [persist] wrote summary.json", flush=True)
    return 0 if gate_passed else 2


if __name__ == "__main__":
    sys.exit(main())
