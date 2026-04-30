"""exp11 — Full 4-region waveform-feature compute.

Applies the validated pipeline (`src/waveform_pipeline.py`) to all qualifying
training regions × all window types (precursor + null A + null B), in
parallel with 4 worker threads.

Per-region FDSN routing (set in `src/regions.py`):
    California → NCEDC
    Cascadia   → IRIS  (UW.LON)
    Turkey     → IRIS  (IU.ANTO)
    Italy      → IRIS  (MN.AQU; pre-reg listed IV.AQU but MedNet is
                        the post-2009 operational continuation)

Window count: ~1,739 across 4 regions. Estimated wall time ~1.5 hours
with 4 workers (vs ~5.6 hours serial).

Outputs:
    <region>_waveform_features.csv  per region
    summary.json
    run.log
"""
from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.regions import (  # noqa: E402
    CALIFORNIA, CASCADIA, ITALY, TURKEY, Region,
)
from src.waveform_pipeline import (  # noqa: E402
    WaveformFeatureParams,
    compute_features_for_dataframe,
)

EXP_DIR = Path(__file__).resolve().parent
SOURCE_CSV = ROOT / "experiments" / "exp07_macro_pra2" / "feature_summary.csv"

REGIONS_TO_RUN: list[Region] = [CALIFORNIA, CASCADIA, TURKEY, ITALY]
N_WORKERS = 4


def _ts() -> str:
    return dt.datetime.now().strftime("%H:%M:%S")


def main() -> int:
    print(f"[{_ts()}] [exp11] start — {N_WORKERS} workers, "
          f"{len(REGIONS_TO_RUN)} regions: "
          f"{[r.name for r in REGIONS_TO_RUN]}", flush=True)

    df_all = pd.read_csv(SOURCE_CSV)
    df_all["t_start"] = pd.to_datetime(df_all["t_start"], utc=True, format="ISO8601")
    df_all["t_end"] = pd.to_datetime(df_all["t_end"], utc=True, format="ISO8601")

    params = WaveformFeatureParams(
        n_snapshots_per_window=6,
        snapshot_seconds=600,
    )

    overall_t0 = dt.datetime.now()
    region_summaries = []
    for region in REGIONS_TO_RUN:
        df_r = df_all[df_all["region"] == region.name].reset_index(drop=True)
        if df_r.empty:
            print(f"[{_ts()}] [exp11] {region.name} no source rows — skip", flush=True)
            continue
        print(f"\n[{_ts()}] [exp11] === {region.name} "
              f"(N={len(df_r)} windows, FDSN={region.fdsn_client}, "
              f"primary={region.primary.code()}) ===", flush=True)

        t0 = dt.datetime.now()
        try:
            out_df = compute_features_for_dataframe(
                region, df_r, params=params, n_workers=N_WORKERS,
                progress_every=20,
                log=lambda m: print(f"[{_ts()}] {m}", flush=True),
            )
        except Exception as e:
            print(f"[{_ts()}] [exp11] {region.name} FAILED: {type(e).__name__}: {e}",
                  flush=True)
            continue
        elapsed = (dt.datetime.now() - t0).total_seconds()
        out_path = EXP_DIR / f"{region.name}_waveform_features.csv"
        out_df.to_csv(out_path, index=False)
        n_succ = int((out_df["wf_n_successful"] > 0).sum())
        n_full = int((out_df["wf_n_successful"] == params.n_snapshots_per_window).sum())
        print(f"[{_ts()}] [exp11] {region.name} done in {elapsed:.0f}s "
              f"({elapsed/max(1,len(df_r)):.1f}s/win)", flush=True)
        print(f"[{_ts()}] [exp11] {region.name}: any-snapshot {n_succ}/{len(df_r)} "
              f"({100*n_succ/max(1,len(df_r)):.0f}%); all-6 {n_full}/{len(df_r)}", flush=True)
        region_summaries.append({
            "region": region.name,
            "fdsn_client": region.fdsn_client,
            "primary_station": region.primary.code(),
            "n_windows": int(len(df_r)),
            "n_any_snapshot_ok": n_succ,
            "n_all_snapshots_ok": n_full,
            "elapsed_seconds": elapsed,
            "elapsed_per_window_s": elapsed / max(1, len(df_r)),
        })

    overall_elapsed = (dt.datetime.now() - overall_t0).total_seconds()
    print(f"\n[{_ts()}] [exp11] all regions done in {overall_elapsed:.0f}s "
          f"({overall_elapsed/60:.1f} min)", flush=True)

    summary = {
        "experiment": "exp11_full_waveform_features",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "n_workers": N_WORKERS,
        "params": {
            "n_snapshots_per_window": params.n_snapshots_per_window,
            "snapshot_seconds": params.snapshot_seconds,
        },
        "regions": region_summaries,
        "total_elapsed_seconds": overall_elapsed,
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[{_ts()}] [persist] summary.json", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
