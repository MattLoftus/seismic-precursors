# exp04 — Parkfield catalog feature distributions (Round B sanity)

**Date:** 2026-04-28 (Session 4)
**Status:** Mechanical sanity **PASS**, scientific sanity **CRITICAL FINDING** — target-set leakage caught early.
**Predecessor:** `experiments/exp03_parkfield_waveform/notes.md`

---

## Goal

Compute catalog-derived precursor features (b-value, b-value drift, Benioff strain total/curvature, n_above_Mc) on real 30-day pre-event windows for Parkfield M≥4.5 targets, compared against random aftershock-free null windows. Validate that:

1. Each feature is *computable* on every precursor window (no NaN-explosion).
2. Distributions are *plottable* — i.e., the pipeline doesn't silently degenerate.
3. Spot-check whether any feature trivially separates precursor from null at AUC=1.0 (which would suggest data leakage).

This is **mechanical sanity, not science**. AUC numbers reported here are NOT pre-registered.

## Method

| Parameter | Value |
|-----------|-------|
| Catalog | exp02 cache, 15,732 events 2000–2024, Parkfield 50 km radius |
| Window length | 30 days |
| Mc for features | 1.0 (relaxed from calibration Mc=1.5; needed for ≥10 events per drift sub-window) |
| Targets | all M ≥ 4.5 in catalog |
| Null sampling | 200 random 30-day windows, requiring ≥30 days before AND ≥60 days after every target ("Null A" semantics) |
| Drift sub-windows | 3 equal-time chunks per 30-day window |
| Benioff energy | log10 E_J = 1.5 M + 4.8; total = Σ √E |
| Per-feature AUC | Mann-Whitney U |

## Headline result

```
Targets found:  4    (all in late September 2004)

  2004-09-28 17:15:24  M5.97  (35.818, -120.366)  depth 8.1 km   ← MAINSHOCK
  2004-09-28 17:24:15  M4.71  (35.804, -120.350)  depth 6.0 km   ← aftershock +9 min
  2004-09-29 17:10:04  M5.00  (35.954, -120.502)  depth 10.8 km  ← aftershock +24 hr
  2004-09-30 18:54:28  M4.88  (35.988, -120.538)  depth 9.9 km   ← aftershock +49 hr
```

Preview AUCs:

| Feature | AUC | Interpretation |
|---------|-----|----------------|
| b-value (Mc=1.0) | 0.212 | Precursor windows have *lower* b. Suggestive but N=4 is too small. |
| b-value drift / day | 0.430 | Indistinguishable from null. |
| log10 Benioff total | **0.995** | **Aftershock contamination — see below** |
| Benioff curvature | **1.000** | **Aftershock contamination** |
| n events ≥ Mc | **0.990** | **Aftershock contamination** |

## Critical scientific finding — target leakage

**Parkfield 50 km radius has only ONE independent M≥4.5 target in 25 years.** The four apparent targets are:

- 1 mainshock (2004-09-28 M5.97)
- 3 aftershocks (within 49 hours of the mainshock, all in the immediate ZBZ-clustered sequence)

For the three aftershock targets, the 30-day pre-event window includes the M5.97 mainshock + thousands of aftershocks. That's why those windows have:

- n_events ≈ 600–850 (vs null median ~50)
- Benioff total ≈ 10^7.2 J^0.5 (vs null median ~10^5.3)
- Benioff curvature ≈ 1.5×10^6 to 8.5×10^6 (vs null near zero)

These are **the mainshock's own seismicity counted as a "precursor" to its aftershocks** — not a real precursor signal. The AUC ≈ 1.0 is data leakage by construction, not skill.

The b-value AUC=0.212 (precursors have lower b) is *suggestive* — Gulia & Wiemer 2019 documented b-decrease before mainshocks — but with N=4 precursor points the statistic is uninterpretable. Even pure noise produces AUC ranging anywhere in [~0.05, 0.95] at this sample size.

## Implications for Round D pre-registration

This is exactly the kind of finding worth catching at exp04, not at peer review. Three pre-registration constraints emerge:

1. **Target events must be declustered.** Use the ZBZ output from exp02-style declustering on the target set itself: keep only background-mode (parent-less) M≥4.5 events as targets. The 3 Sept 2004 aftershocks should be dropped, leaving the M5.97 mainshock as the single target. (Pre-reg should commit to ZBZ thresholding on targets exactly as on the catalog at large.)
2. **Regional scope must be broadened.** Parkfield-only is hopeless: 1 independent target in 25 years. The pre-registered "California" training region should encompass the entire San Andreas + Long Valley + Eastern California shear zone — likely O(50–100) independent declustered M≥4.5 in 2000–2024 across that scope. (PLAYBOOK §5.4 already lists "California" as one of the 6 training regions; we just need to be precise about its bounds.)
3. **30-day pre-event windows must exclude post-mainshock periods of OTHER recent events.** Even after target declustering, a precursor window may overlap with an unrelated nearby mainshock's aftershock sequence. The exclusion buffer (currently 60 days post-target in Null A) needs to be applied to PRECURSOR windows too, not just null windows. Otherwise we get the same leakage by a less obvious path.

These three are now non-negotiable pre-registration commitments.

## Mechanical sanity — pipeline status

| Check | Result |
|-------|--------|
| Catalog re-loads from cache via ISO8601 parsing | ✅ (after fixing mixed-microsecond format issue in src/data.py + this exp's loader) |
| All 4 catalog features computable on every precursor window | ✅ 4/4 |
| All 4 features computable on null windows | ✅ b: 130/200, b_drift: 114/200, Benioff: 200/200 (b feature requires ≥30 events≥Mc; sparser null windows give NaN, this is correct) |
| b-value drift requires ≥10 events per sub-window | ✅ enforced; sparser windows return NaN |
| Plots and CSV produced | ✅ |

The pipeline mechanics work. The science we have to throw away on this region is informative.

## Files produced

- `feature_summary.csv` — 204 rows × ~12 columns; per-window features
- `feature_distributions.png` — 2x3 grid of histograms (precursor red, null gray) + AUC table panel
- `summary.json` — machine-readable, includes the preview AUCs (clearly tagged NOT pre-registered)
- `run.py`, `run.log`

## Bugs / cleanup

- Fixed: mixed-microsecond ISO8601 parsing in `src/data.py` and this experiment's loader.
- Done: `src/features/preprocessing.py` now clamps `freqmax` to `fs/3` per trace and uses `inventory=` rather than the deprecated `attach_response=True`.
- Open: the AUC computation here uses scipy-free Mann-Whitney; will move to `scipy.stats.mannwhitneyu` once we want exact p-values.

## What's left for Round B

| Feature | Module | Status |
|---------|--------|--------|
| b-value (Aki MLE) | `src/features/bvalue.py` | ✅ exp01/02 |
| b-value drift | `src/features/bvalue.py` | ✅ exp04 |
| Benioff strain | `src/features/benioff.py` | ✅ exp04 |
| Spectral slope | `src/features/spectral.py` | ✅ implemented; not yet exercised on real data (needs waveforms) |
| Repeating-event rate | `src/features/repeating.py` | scaffold only — full catalog application requires waveform-pull infrastructure (Round C) |
| Waveform entropy | `src/features/entropy.py` | ✅ exp03 |
| HHT IMF1 IF | `src/features/hht.py` | ✅ exp03 |
| Declustering | `src/features/declustering.py` | ✅ exp02 |

## Next session (Session 5)

1. **Pre-registration draft.** Now is the right time. We have all 6 features implemented; Mc choice is calibrated; target leakage is identified; null sampling protocol is exercised. Draft `papers/pre_registration.md` committing: feature stack, Mc=1.5 background calibration / Mc=1.0 feature compute, declustered targets, 6 train + 2 test regions with explicit lat/lon polygons, 3-null-A/B/C protocol, ROC-AUC primary metric, 3σ bootstrap threshold. Hash the commit SHA.
2. **Broaden region definition.** Replace "Parkfield 50 km" with "California" (e.g., 32–42 N, 125–114 W) for the cross-regional pipeline. exp05: re-run feature distributions with declustered targets across California to see what real signal remains after the Parkfield-only artifacts are removed.
3. Carryover: download Mizrahi 2024 PDF manually before Round E.
