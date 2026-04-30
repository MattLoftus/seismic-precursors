"""exp10 — California waveform-feature validation (precursor windows only).

Validates the waveform-feature pipeline on California precursor windows from
exp07. 37 windows × 6 snapshots × 10 min ≈ 222 FDSN fetches, ~25-50 min wall.

If this completes successfully and gives sensible feature distributions:
    - timing extrapolation for the full 4-region run is established
    - the per-window aggregation strategy (mean/median across 6 snapshots)
      is validated
    - we can move to Session 13 full run with confidence

If it fails (e.g., FDSN rate-limit, BK.PKD coverage gaps):
    - station-availability fallback to backups becomes a Session 13 priority
    - or we reduce snapshot count / duration

This is NOT pre-registered as a separate experiment; it is the first step of
implementing the Round C waveform sub-protocol that the pre-registration
deferred.

Outputs:
    california_precursor_waveform_features.csv
    summary.json
    run.log
"""
from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.regions import CALIFORNIA  # noqa: E402
from src.waveform_pipeline import (  # noqa: E402
    WaveformFeatureParams,
    compute_features_for_dataframe,
)

EXP_DIR = Path(__file__).resolve().parent
SOURCE_CSV = ROOT / "experiments" / "exp07_macro_pra2" / "feature_summary.csv"


def _ts() -> str:
    return dt.datetime.now().strftime("%H:%M:%S")


def main() -> int:
    print(f"[{_ts()}] [exp10] California waveform-feature validation", flush=True)
    if not SOURCE_CSV.is_file():
        print(f"[{_ts()}] missing {SOURCE_CSV}", flush=True)
        return 2

    df_all = pd.read_csv(SOURCE_CSV)
    df_all["t_start"] = pd.to_datetime(df_all["t_start"], utc=True, format="ISO8601")
    df_all["t_end"] = pd.to_datetime(df_all["t_end"], utc=True, format="ISO8601")
    df = df_all[(df_all["region"] == "California") &
                (df_all["window_kind"] == "precursor")].reset_index(drop=True)
    print(f"[{_ts()}] [exp10] {len(df)} California precursor windows", flush=True)
    print(f"[{_ts()}] [exp10] window range: {df['t_start'].min()} .. {df['t_end'].max()}",
          flush=True)
    print(f"[{_ts()}] [exp10] primary station: {CALIFORNIA.primary.code()}", flush=True)

    params = WaveformFeatureParams(
        n_snapshots_per_window=6,
        snapshot_seconds=600,    # 10 min
        fdsn_client="EARTHSCOPE",
        bandpass_freqmin=0.5,
        bandpass_freqmax=12.0,
    )
    print(f"[{_ts()}] [exp10] params: {params.n_snapshots_per_window} snapshots × "
          f"{params.snapshot_seconds}s, bandpass [{params.bandpass_freqmin}, "
          f"{params.bandpass_freqmax}] Hz", flush=True)

    t0 = dt.datetime.now()
    out_df = compute_features_for_dataframe(
        CALIFORNIA, df, params=params, progress_every=5,
        log=lambda m: print(f"[{_ts()}] {m}", flush=True),
    )
    elapsed = (dt.datetime.now() - t0).total_seconds()
    print(f"[{_ts()}] [exp10] done in {elapsed:.0f}s "
          f"({elapsed/max(1,len(df)):.1f}s per window)", flush=True)

    out_path = EXP_DIR / "california_precursor_waveform_features.csv"
    out_df.to_csv(out_path, index=False)
    print(f"[{_ts()}] [persist] {out_path.name} ({len(out_df)} rows)", flush=True)

    n_succ = int((out_df["wf_n_successful"] > 0).sum())
    n_full = int((out_df["wf_n_successful"] == params.n_snapshots_per_window).sum())
    print(f"\n[{_ts()}] [summary] windows with >=1 successful snapshot: "
          f"{n_succ}/{len(out_df)} ({100*n_succ/max(1,len(out_df)):.0f}%)", flush=True)
    print(f"[{_ts()}] [summary] windows with all {params.n_snapshots_per_window} "
          f"successful: {n_full}/{len(out_df)}", flush=True)
    print(f"[{_ts()}] [summary] feature distributions:", flush=True)
    for col in ("wf_spectral_slope", "wf_spectral_r2",
                "wf_waveform_entropy", "wf_hht_imf1_if_median_hz"):
        vals = out_df[col].to_numpy(dtype=float)
        finite = vals[np.isfinite(vals)]
        if len(finite):
            print(f"[{_ts()}]   {col:<32s} N={len(finite):>3d}  "
                  f"med={np.median(finite):.3f}  IQR=[{np.percentile(finite,25):.3f}, "
                  f"{np.percentile(finite,75):.3f}]", flush=True)
        else:
            print(f"[{_ts()}]   {col:<32s} N=0 finite", flush=True)

    summary = {
        "experiment": "exp10_california_waveform_validation",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "elapsed_seconds": elapsed,
        "elapsed_per_window_s": elapsed / max(1, len(df)),
        "params": {
            "n_snapshots_per_window": params.n_snapshots_per_window,
            "snapshot_seconds": params.snapshot_seconds,
            "fdsn_client": params.fdsn_client,
            "bandpass": [params.bandpass_freqmin, params.bandpass_freqmax],
        },
        "n_windows": int(len(df)),
        "primary_station": CALIFORNIA.primary.code(),
        "n_windows_any_snapshot_ok": n_succ,
        "n_windows_all_snapshots_ok": n_full,
        "feature_summaries": {
            col: {
                "n_finite": int(np.isfinite(out_df[col]).sum()),
                "median": float(np.nanmedian(out_df[col])) if int(np.isfinite(out_df[col]).sum()) else None,
                "p25": float(np.nanpercentile(out_df[col], 25)) if int(np.isfinite(out_df[col]).sum()) else None,
                "p75": float(np.nanpercentile(out_df[col], 75)) if int(np.isfinite(out_df[col]).sum()) else None,
            }
            for col in ("wf_spectral_slope", "wf_spectral_r2",
                        "wf_waveform_entropy", "wf_hht_imf1_if_median_hz")
        },
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[{_ts()}] [persist] summary.json", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
