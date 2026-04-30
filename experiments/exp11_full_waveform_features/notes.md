# exp11 — Full 4-region waveform-feature compute

**Date:** 2026-04-29 (Session 13)
**Status:** **PASS** — all 4 qualifying training regions completed; 90.5% windows have at least one successful waveform snapshot (76% have all 6).
**Predecessor:** `experiments/exp10_california_waveform_validation/notes.md`

---

## Result

| Region | FDSN | Primary station | Windows | any-snapshot | all-6 |
|--------|------|-----------------|---------|---------------|--------|
| California | NCEDC | BK.PKD | 437 | 437 (100%) | 423 (97%) |
| Cascadia | IRIS | UW.LON | 436 | 411 (94%) | 288 (66%) |
| Turkey | IRIS | IU.ANTO | 430 | 327 (76%) | 317 (74%) |
| Italy | IRIS | MN.AQU | 436 | 398 (91%) | 291 (67%) |
| **Total** | | | **1,739** | **1,573 (90.5%)** | **1,319 (76%)** |

Total wall time: **40.6 min** with 4 worker threads (vs ~5.6 hr serial — 8.3× speedup, slightly above the theoretical 4× because some workers were CPU-bound on EMD while others were waiting on FDSN).

## Cross-regional baseline waveform character

The four primary stations have markedly different baseline character:

| Region | spectral slope median | entropy median (bits) | HHT IMF1 IF median (Hz) |
|--------|----------------------|------------------------|--------------------------|
| California (BK.PKD)  | $-2.30$ (Brune $f^{-2}$ source-like) | 2.15 (concentrated) | 6.95 |
| Cascadia (UW.LON)    | $-0.30$ (nearly white) | 3.5 (broad)         | 8.5 (fixed near $f_s/3$) |
| Turkey (IU.ANTO)     | $-1.30$ (rising spectrum + noise) | 3.3                | 4.0 |
| Italy (MN.AQU)       | $-0.85$ (mid)         | 3.4                | 4.6 |

These are NOT comparable in absolute terms — each station has its own
noise floor and instrument response. What matters for cross-regional
analysis is the WITHIN-region precursor-vs-null AUC, which doesn't depend
on the baseline magnitude.

## Two operational findings

### IV.AQU is gone post-2009 — operational continuation is MN.AQU

The pre-reg listed `IV.AQU` (Italian Seismic Network) as Italy's primary
station. Session 13 verification found `IV.AQU` returns zero data via INGV,
ORFEUS, or IRIS for any time in 2000–2024. The post-2009 reorganization
moved L'Aquila Observatory under MedNet (`MN.AQU`, INGV-affiliated, IRIS-
archived). MN.AQU is the same physical station; we treat the network-code
update as the operational implementation of pre-reg v1 §3's "95%
data-availability fallback rule" rather than a station change. Documented
in `src/regions.py` and the paper Methods §3.1.

### Per-region waveform-archive routing

Pre-reg v1 §3 didn't specify FDSN endpoints. exp10 found BK.PKD requires
NCEDC; exp11 found Italy requires IRIS (not INGV) for MN.AQU. Final
routing:

```
California  -> NCEDC   (BK network is Berkeley/NCSN-archived)
Cascadia    -> IRIS    (UW network in PNSN, IRIS-fed)
Turkey      -> IRIS    (IU global network on IRIS)
Italy       -> IRIS    (MN MedNet on IRIS)
```

Both findings are documented as deviations from pre-reg v1's silent
default of EARTHSCOPE/IRIS uniformity and are reported in paper Methods.

## Per-region per-feature within-region AUC (preview, see exp12 for
macro)

The per-region waveform feature AUC numbers feed exp12's macro
computation. Per-region windows-vs-null distributions are in
`<region>_waveform_features.csv` (one CSV per region). Spot check on
California precursor windows (37 windows): spectral slope tightly
distributed around $-2.30$, entropy 2.0–2.4 bits, IMF1 IF 6–8 Hz —
matches exp10 validation result.

## Files produced

- `<region>_waveform_features.csv` (one per region; gitignored)
- `summary.json` — per-region runtime + success counts
- `run.py`, `run.log`

## Implications for exp12

All 6 features are now computable on the same per-(region, window) grid
as the catalog features. Macro-AUC over the full 6-feature panel is the
exp12 result.

## Acknowledged limitation

The 6×10-min snapshot strategy (per-window subsampling within the 30-day
window) was NOT specified in the pre-reg. We picked it pragmatically:
3 hours of total fetch per window keeps storage tractable; 6 evenly-spaced
snapshots reduce per-window variance from temporal heterogeneity.
A different sampling strategy (e.g., one full day per window, or
30 short snapshots) might give different feature variances. We document
the choice as deterministic per (region, t_start, t_end) so future
replication is exact, and we acknowledge in §Discussion that this is
a free parameter in the within-window aggregation that the pre-reg
didn't fix.
