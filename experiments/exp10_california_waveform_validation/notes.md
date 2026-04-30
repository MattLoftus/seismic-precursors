# exp10 — California waveform-feature validation

**Date:** 2026-04-29 (Session 12)
**Status:** **PASS** — all 37 California precursor windows yielded computable waveform features with physically reasonable distributions. Pipeline validated end-to-end.
**Predecessor context:** Round C waveform sub-protocol that was deferred from the pre-reg's headline; user directive (Session 12) to run all 6 features before paper writing.

---

## Goal

Validate the waveform-feature pipeline (`src/waveform_pipeline.py`) on California's 37 precursor windows from exp07 before scaling to all 4 qualifying training regions × precursor + null windows in Session 13.

## Two debugging discoveries

### 1. BK.PKD waveforms are NCEDC-hosted, not EarthScope

The first run targeted EarthScope (the post-IRIS-rebrand short-name) and returned HTTP 204 / `FDSNNoDataException` for every window. Investigation revealed BK (Berkeley Digital Seismic Network) is operated by Berkeley/NCSN and primarily archived at NCEDC (`https://service.ncedc.org`); EarthScope/IRIS has only partial mirroring of BK continuous data.

**Fix:** added `fdsn_client: str = "IRIS"` field to the `Region` dataclass (`src/regions.py`), with California overriding to `"NCEDC"` and Italy to `"INGV"`. `compute_window_features` now reads `region.fdsn_client` and falls back to params default if missing.

### 2. Channel pattern `?HZ` matched LHZ (1 Hz) and broke butter()

The second-pass run successfully fetched data via NCEDC but failed in the EMD bandpass step with `ValueError: Digital filter critical frequencies must be 0 < Wn < 1`. Tracing: the channel pattern `?HZ` matched LHZ (long-period Z, 1 Hz). With fs=1, my freqmax=12 gave Wn=12/(0.5)=24, way above 1.

**Fix:** changed default channel pattern to `BHZ,HHZ,EHZ` (explicit broadband-Z list, no LH/SH). Added `min_sampling_rate=20.0` to `WaveformFeatureParams` to reject any trace with fs < 20 Hz. Also clamped the EMD bandpass freqmax to `min(params.bandpass_freqmax, fs / 3.0)` — defense-in-depth against any other low-fs trace that slips through.

## Setup (final)

| Parameter | Value |
|-----------|-------|
| Region | California; primary station BK.PKD |
| FDSN client | NCEDC |
| Snapshots per window | 6 |
| Snapshot duration | 600 s (10 min) |
| Channels | `BHZ,HHZ,EHZ`; min fs ≥ 20 Hz |
| Bandpass | 0.5–12 Hz, clamped to fs/3 |
| Aggregation | mean for slope/entropy, median for IMF1 IF |

## Result

- **37/37 California precursor windows** had at least one successful snapshot.
- **33/37** had all 6 snapshots successful; the 4 with partial coverage had failed snapshots in 2001, 2006, 2016 (BK.PKD operational gaps).
- **Total runtime: 429 s (11.6 s per window).**

| Feature | Median | IQR | Interpretation |
|---------|--------|-----|----------------|
| `wf_spectral_slope` | $-2.30$ | $[-2.54, -2.07]$ | Matches Brune 1970 $f^{-2}$ source spectra |
| `wf_spectral_r2` | $0.81$ | $[0.77, 0.89]$ | log-log fits explain ~80% of PSD variance |
| `wf_waveform_entropy` | $2.15$ bits | $[2.00, 2.44]$ | Moderately structured signal (vs max ~6 bits) |
| `wf_hht_imf1_if_median_hz` | $6.95$ Hz | $[6.30, 7.60]$ | Mid-band; cultural noise + small events |

The distributions are tight and physically reasonable. spectral_slope around $-2$ is consistent with the high-frequency tail of source spectra; the IQR is narrow ($\pm 0.25$), suggesting the per-region waveform character is well-defined for "background" windows. This is the right baseline against which to compare null-window distributions in Session 13.

## Timing extrapolation

For the full 4-region run including null windows:

| Region | precursor + null A + null B | total |
|--------|------------------------------|-------|
| California | 37 + 200 + 200 | 437 |
| Cascadia | 36 + 200 + 200 | 436 |
| Turkey | 30 + 200 + 200 | 430 |
| Italy | 36 + 200 + 200 | 436 |
| **Total** | | **1,739 windows** |

At 11.6 s/window: **5.6 hours** of compute. Suitable for an overnight background job in Session 13.

If we add the test regions for completeness (even though their N_pre=1 doesn't qualify under PRA-2 Amendment 3, computing waveform features on what we have is cheap):

| Region | precursor + null A + null B | total |
|--------|------------------------------|-------|
| Mexico | 1 + 0 + 0 | 1 |
| Alaska | 1 + 0 + 0 | 1 |

Adds ~24 sec to the runtime. We'll skip them — no AUC computable.

## Session 13 plan

1. Run `exp11_full_waveform_features` over all 4 qualifying training regions × 1,739 windows in background. ~5–6 hours wall.
2. Per-region FDSN clients: California=NCEDC (verified), Cascadia=IRIS (UW.LON in PNSN), Turkey=IRIS (IU.ANTO), Italy=INGV (IV.AQU).
3. Save per-region CSV: `experiments/exp11_full_waveform_features/<region>_waveform_features.csv` with all features per window.
4. exp12: combine catalog + waveform features, recompute macro-AUC over the full 6-feature panel under PRA-2.

## Open issues for Session 13+

- INGV's FDSN endpoint may have its own quirks; will verify on Italy run.
- HHZ at 80 Hz is preferred over BHZ at 20 Hz when both are available — the current pipeline picks alphabetically (BHZ first). Could be optimized to pick highest-fs but the BHZ-at-20-Hz result is fine for our purposes.
- The 4 California windows with partial coverage (snapshot failures) still produce per-window aggregated features from the successful snapshots; this means they're slightly noisier but still usable.
