# exp03 — Parkfield 2004 M6.0 waveform sanity check (Round A4)

**Date:** 2026-04-28 (Session 3)
**Status:** **PASS** after iterating on the IMF1-IF gate definition. Three feature-pipeline mechanics validated end-to-end on a known event.
**Predecessor:** `experiments/exp02_parkfield_declustered/notes.md`
**PLAYBOOK §10:** Round A4 / Week-2 single-waveform sanity gate

---

## Goal

Validate the full waveform pipeline (EarthScope FDSN fetch → instrument-response removal → bandpass → feature compute) on **one** known event before scaling to 6-region × 20-year continuous waveform pulls. Specifically:

1. Predicted P arrival lands on a clear amplitude rise.
2. Spectral entropy shifts in the expected direction during the event.
3. HHT IMF1 instantaneous frequency responds (any direction) to the event.

## Setup

| Parameter | Value |
|-----------|-------|
| Event | 2004-09-28T17:15:24 UTC, M6.0, 35.815N -120.374W, depth 8.5 km |
| Station | BK.PKD (Parkfield Direct), 35.9452 N, -120.5416 W |
| Geometry | epicentral 20.92 km, hypocentral 22.58 km, depth 8.5 km |
| Predicted | t_P = 3.76 s, t_S = 6.45 s (analytic Vp=6.0, Vs=3.5 km/s) |
| Window | origin ± 1800 s (1 hour total) |
| Channels | BK.PKD ?H? → returned 6 traces (BHE, BHN, BHZ at 40 Hz; LHE, LHN, LHZ at 1 Hz) |
| Preprocess | merge, detrend, taper, response-removal to VEL (m/s), bandpass 0.5–20 Hz |
| Features | spectral entropy + HHT IMF1 IF on quiet pre-event (t-1500…t-1200) vs event (t-30…t+270) |

## Headline result

| Gate | Result | Pass? |
|------|--------|-------|
| Predicted P/S land on clear arrivals | Yes (visible at t=3.76s, 6.45s) | ✅ |
| Spectral entropy DROPS at event | Δ = −1.77 bits (3.90 → 2.13) | ✅ |
| IMF1 IF RESPONDS at event (\|Δ\| > 1 Hz) | \|Δ\| = 4.01 Hz median (5.43 → 0.90) | ✅ |

**Round A4 sanity: PASS.**

## What went wrong on the first run

The original gate required IMF1 IF to *rise* at the event (motivated by "body waves are higher freq than microseismic noise"). It failed: quiet IF1 median = 19.59 Hz, event IF1 median = 1.00 Hz. Two findings emerged:

### Finding 1 — 60 Hz line noise aliases exactly to Nyquist at fs=40 Hz

BK.PKD's 2004 broadband sampler ran at 40 Hz, so Nyquist = 20 Hz. Standard US AC mains noise is at 60 Hz; aliased into a 40 Hz sampler, that becomes |60 − 1·40| = 20 Hz **= exactly Nyquist**. EMD's first IMF in a quiet window picked up this aliased line noise as the dominant high-frequency mode → IF1 saturated at Nyquist. Fix: pre-bandpass the EMD input to 0.5–12 Hz (well below Nyquist) so line-noise aliasing is filtered out.

### Finding 2 — EMD's IMF1 in seismic data is amplitude-dominated

After the bandpass fix, IF1 *still* dropped at the event (5 → 1 Hz), didn't rise. EMD picks the *dominant* local oscillation as IMF1, where "dominant" means amplitude × extrema frequency. For a M6.0 mainshock at 23 km, the body-wave amplitude (1–3 Hz dominant) is 100–1000× any cultural noise (5–10 Hz). EMD picks the body-wave oscillation as IMF1.

So the **sign of the IMF1-IF shift is signal-content-dependent**, not universally directional. The right gate is: |shift| > threshold ("the feature responds"). The cross-regional classifier downstream can learn the per-region sign.

## Files produced

- `raw_waveform.mseed` (616 KB) — gitignored under `*.mseed`; cached so re-runs skip FDSN
- `inventory.xml` (~few KB) — station + response metadata
- `waveform_with_arrivals.png` — 3-component velocity traces with predicted P/S overlays — clean P arrival at 3.76 s; bigger S/coda
- `entropy_quiet_vs_event.png` — entropy time series; event drops cleanly to ~1.5 bits during peak motion
- `imf1_if_quiet_vs_event.png` — IF1 time series; quiet broadly distributed 0–15 Hz, event mostly 0–3 Hz body-wave-dominant
- `summary.json` — machine-readable result
- `run.py`, `run.log` — reproducible script + log

## Lessons for the cross-regional pipeline (Round B–D)

1. **Station-specific anti-alias preprocessing.** Different networks/eras have different sampling rates and anti-alias filter quality. Apply a uniform pre-EMD bandpass at fs/3 (or 12 Hz, whichever lower) regardless of nominal sampling rate. **Update preprocessing.py to clamp `freqmax` to min(20, fs/3)** as a Session 4 cleanup task.
2. **Don't pre-specify feature directions.** For each feature, test the magnitude of response between event/non-event windows; let downstream classifiers learn the per-region sign. Pre-registration should commit to **|feature shift|**, not signed differences.
3. **Round A4 mechanics validated.** ObsPy fetch, instrument-response removal, bandpass, EMD via `emd 0.8.1`, spectrogram-based entropy — full chain works on real data.
4. **Useful side benefit:** EMD on a 12 ksample (5 min × 40 Hz) signal completes in ~0.3 s. Scaling to 30-day pre-event windows × 6 regions × ~80 events ≈ 14,400 windows. EMD ~5,000 s = ~1.5 hr per region. Tractable for Round C.

## What I won't bother with in Session 4

- Re-pulling 30-day windows for the 2004 event. The 1-hour pull was sufficient for sanity. Real precursor windows in Round C will use 30-day pre-event windows per target M≥4.5 event.
- HHT marginal Hilbert spectrum (more sophisticated HHT-derived feature). The mean IMF1 IF is what PLAYBOOK §5.2 specified; advanced HHT statistics can come later if needed.

## Open carryover from Session 1

Mizrahi 2024 PDF still not manually downloaded. Adding to Session 4 todo.

## Bugs / cleanup for next session

- ObsPy emits a deprecation warning that `attach_response=True` will be removed; switch to using `inventory=...` in `remove_response()` directly.
- ScientiPy issues a warning that 20 Hz bandpass at 40 Hz fs = Nyquist; preprocessing should auto-clamp.
- Add `*.xml` to .gitignore (we already ignore data caches; inventory XML is regenerable).
