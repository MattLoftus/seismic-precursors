# Pre-Registration — Cross-Regional Null-Corrected Earthquake Precursor Search

**Project:** seismic-precursors (Cedar Loop LLC)
**Date committed:** 2026-04-28
**Author:** Matthew Loftus
**Repo:** https://github.com/MattLoftus/seismic-precursors
**This document:** `papers/pre_registration.md`
**Pre-reg commit SHA:** *(written into PLAYBOOK and any later paper draft after this file is committed; this SHA is the canonical reference for the pre-registered protocol)*

---

## 0. Why this document exists

The earthquake-precursor literature has been burned repeatedly by post-hoc feature tuning, region cherry-picking, and silent Mc adjustments. The pre-registration protocol used by CSEP (Schorlemmer+ 2007; Mizrahi+ 2024) and adopted in adjacent-domain work (e.g., the MDPI 2025 ionospheric-precursor study) is the cleanest defense against these failure modes.

This file commits, **before any cross-regional evaluation has been run**, every choice that could be tuned post-hoc to inflate apparent skill. After this commit, the protocol is locked. The commit SHA will be cited in any subsequent paper draft as the reference point for "what we said we would do."

The Session 4 (`experiments/exp04_parkfield_feature_distributions/notes.md`) finding — that naive precursor-window definition gives AUC ≈ 1.0 on Benioff/n_events through aftershock self-contamination — is the immediate motivation. Three commitments emerge as non-negotiable and are incorporated below.

---

## 1. Catalog protocol

### 1.1 Source
- **Catalog:** ANSS ComCat via `usgs-libcomcat` Python client (`libcomcat.search.search`).
- **Magnitude type:** prefer M_w; fall back to whatever ComCat returns as the headline magnitude (this is a known limitation of multi-decade global catalogs and we cannot fix it post-hoc).
- **Time range:** 2000-01-01 00:00:00 UTC to 2024-12-31 23:59:59 UTC.
- **Lower magnitude cutoff for catalog ingestion:** M ≥ 1.5 (deliberately below typical regional Mc so we can estimate Mc ourselves per region).

### 1.2 Magnitude of completeness (Mc) per region
Mc is estimated **per region** using:

1. Wiemer & Wyss 2000 maximum-curvature on a 0.1-magnitude binning.
2. Woessner & Wiemer 2005 +0.20 correction.
3. **Mandatory diagnostic:** plot the b-vs-Mc curve over Mc ∈ {1.0, 1.15, 1.25, 1.35, 1.50, 1.65, 1.80}. If b is still rising at the chosen Mc (no plateau), increment Mc by 0.15 until the plateau is reached or 5 sub-windows of 30-day data have fewer than 30 events ≥ Mc — whichever happens first.
4. The chosen Mc per region is recorded as an artifact, not tuned downstream.

### 1.3 Catalog-level declustering
- **Method:** Zaliapin & Ben-Zion 2013 nearest-neighbor (`src/features/declustering.py`).
- **Parameters:** `b = 1.0` (ZBZ 2013 published default), `d_f = 1.6`, `p = q = 0.5`.
- **Threshold:** `log10(η_0) = −5.0` (ZBZ 2013 published default).
  - Auto-threshold via KDE-trough is computed and reported as a sensitivity check; it is **NOT** the headline.
- Background events (NN proximity ≥ η_0) are kept; clustered events are dropped from the catalog used for catalog-derived features.

### 1.4 Target events
- **Magnitude cutoff:** M ≥ 4.5 (PLAYBOOK §2).
- **Target declustering (NEW Session 4 commitment):** the ZBZ output above is also applied to the target set. Only background-mode (parent-less) M ≥ 4.5 events are kept. Aftershocks of other targets are dropped.
- **Minimum target count per region:** ≥ 8 declustered targets. Regions with fewer are dropped from the panel and noted as a limitation.

---

## 2. Region polygons (LOCKED)

### 2.1 Training regions (6)
Bounding boxes (lat_min, lat_max, lon_min, lon_max), inclusive:

| Region | lat_min | lat_max | lon_min | lon_max | Tectonic regime |
|--------|---------|---------|---------|---------|-----------------|
| California | 32.0 | 42.0 | −125.0 | −114.0 | Strike-slip (San Andreas) + thrust margins |
| Cascadia | 40.0 | 50.0 | −130.0 | −121.0 | Subduction (Juan de Fuca + Pacific) |
| Japan | 30.0 | 46.0 | 128.0 | 148.0 | Subduction (Pacific + Philippine) |
| Chile | −45.0 | −18.0 | −76.0 | −68.0 | Subduction (Nazca) |
| Turkey | 36.0 | 42.0 | 26.0 | 44.0 | Strike-slip (NAF + EAF) + thrust |
| Italy | 36.0 | 47.0 | 6.0 | 19.0 | Continental thrust + extension |

### 2.2 Test regions (held out) (2)
| Region | lat_min | lat_max | lon_min | lon_max | Tectonic regime |
|--------|---------|---------|---------|---------|-----------------|
| Mexico | 14.0 | 32.0 | −118.0 | −86.0 | Subduction (Cocos) + Gulf of California rift |
| Alaska | 51.0 | 72.0 | −180.0 | −130.0 | Subduction (Pacific) + transform |

### 2.3 Pre-registered caveat about test set tectonics
Mexico and Alaska are tectonically **subduction-similar** to Chile / Japan / Cascadia in the training set. This is a softer cross-regime extrapolation than e.g. Iran (continental thrust) + Indonesia (volcanic arc) would have been. The result will be reported as "cross-regional within plate-boundary settings," **not** as "universal across all tectonic regimes." This caveat is committed here so it cannot be quietly omitted from the paper.

### 2.4 Region locking
The above 6 + 2 set is **frozen**. We will NOT swap regions if a particular train/test split looks unfavorable post-evaluation.

---

## 3. Per-region station selection (waveform-derived features)

The waveform-derived features (spectral slope, waveform entropy, HHT IMF1 IF) require a representative broadband station per region. The primary station per region is fixed below; one or two backups are listed for the case of insufficient data availability.

| Region | Primary | Backup 1 | Backup 2 |
|--------|---------|----------|----------|
| California | BK.PKD | BK.CMB | IU.COR |
| Cascadia | UW.LON | IU.COR | UW.RATT |
| Japan | II.MAJO | IU.MAJO | II.INU |
| Chile | IU.LCO | C1.GO01 | G.PEL |
| Turkey | IU.ANTO | KO.ANTO | II.KIV |
| Italy | IV.AQU | IV.MGAB | G.SSB |
| Mexico | IU.TEIG | IU.TUC | G.UNM |
| Alaska | IU.COLA | IU.KDAK | AT.KIAG |

**Fallback rule:** if a primary station has < 95% data availability over 2000-2024 within ±30 days of any target, use the highest-availability backup. This is decided **per region, once**, before features are computed; no swapping during evaluation.

---

## 4. Feature stack (LOCKED)

Six features per 30-day window. Implementations are pinned to commit `b188ce3` (Round B closure):

| # | Feature | Module | Per-window scalar(s) |
|---|---------|--------|----------------------|
| 1 | **b-value** | `src/features/bvalue.py::aki_bvalue` | b at the per-region Mc |
| 2 | **b-value drift** | `src/features/bvalue.py::bvalue_drift` | linear-fit slope of b across 3 equal-time sub-windows (per day) |
| 3 | **Benioff strain** | `src/features/benioff.py::benioff_features` | log10(Σ √E) and quadratic-fit cumulative curvature |
| 4 | **Spectral slope** | `src/features/spectral.py::spectral_slope` | log-log Brune slope on Welch-PSD aggregated over the 30-day window's hourly subsamples (band 0.5 Hz to min(20, fs/3)) |
| 5 | **Waveform entropy** | `src/features/entropy.py::windowed_mean_entropy` | mean Shannon entropy (bits) over hourly subsampled segments |
| 6 | **HHT IMF1 IF** | `src/features/hht.py::imf1_if_series` | median IMF1 instantaneous frequency over hourly subsampled segments, after 0.5 Hz to fs/3 bandpass |

**No feature additions post-hoc.** If, after this commit, an additional candidate feature is identified, it cannot be added to the headline evaluation. It can be reported in a separate "exploratory" section that is explicitly NOT subject to the 3σ test.

**No feature definition tweaks post-hoc.** The Mc used for b-value is the region's Mc from §1.2. The drift sub-windows are 3, not 5. The Benioff energy formula is `log10 E_J = 1.5 M + 4.8`. These are locked.

### 4.1 Repeating-event rate — deferred
A seventh feature (Nadeau-McEvilly repeating-event rate) was scaffolded in `src/features/repeating.py` but its full at-scale catalog application requires waveform fetch for every event in every window — out of scope for this evaluation. **It is excluded from the pre-registered headline evaluation.** If we later add a Round-C waveform-pull infrastructure that makes it tractable, repeating-event rate will be reported as an exploratory feature, not in the headline 6-feature panel.

---

## 5. Window definitions (LOCKED)

### 5.1 Precursor windows
- **Length:** 30 days.
- **Definition:** for each declustered target with origin time `t`, the precursor window is `[t − 30 days, t)`.
- **Rejection rule (NEW Session 4 commitment):** the precursor window is rejected if any portion of it overlaps another target's `[t' − 30 days, t' + 60 days]` interval. (Both ends of the buffer apply: 30-day pre-window and 60-day post-event aftershock-tail buffer.)

### 5.2 Null A — aftershock-free random
- **Method:** uniformly random 30-day windows over the catalog time-span.
- **Acceptance rule:** window must be at least 30 days BEFORE every target and at least 60 days AFTER every target.
- **Per-region count:** 200 null windows.
- **Same rejection rule as precursor windows** (no overlap with another target's [−30, +60] zone).

### 5.3 Null B — random window with minimal exclusion
- **Method:** uniformly random 30-day windows over the catalog time-span.
- **Acceptance rule:** window must not contain any M ≥ 4.5 event AND must be at least 7 days after such an event (to remove the immediate strong-motion aftershock contamination on waveform features).
- **Per-region count:** 200 null windows.
- This is a less-strict null than A; it admits steady-state aftershock-rich periods that A excludes.

### 5.4 Null C — colored-noise synthetic
- **Method:** for each waveform-derived feature, generate synthetic 30-day records with the same power spectrum (1/f^α best-fit α from a year of station noise during 2010, a tectonically quiet year for that region) but randomized phases.
- **Per-region count:** 200 null synthetic records.
- **Catalog features (b, b-drift, Benioff)** are not computable on synthetic noise; for those, Null C is **not applied** (exclusion noted).

---

## 6. Evaluation protocol (LOCKED)

### 6.1 Primary metric
**Per-feature ROC-AUC** computed via Mann-Whitney U statistic.

For each (training region × feature × null type), compute AUC over precursor vs null windows. Then aggregate across the 6 training regions: macro-average and pool. The cross-regional headline AUC is the **macro-averaged** value per feature.

### 6.2 Held-out test
After training-region AUCs are computed and before any tuning, the same protocol is applied to Mexico + Alaska. The test-region AUC is the **headline result**.

### 6.3 Significance threshold
Bootstrap (`n_resamples = 1000`) over (precursor windows, null windows). The **3σ criterion** for "feature works" is:

   **bootstrap CI95 lower bound > 0.5  AND  rel-effect (|AUC − 0.5|) > 0.05  AND  z-score > 3 vs the null-of-the-null**

The "null-of-the-null" is computed by permuting precursor / null labels within each region and recomputing macro-AUC; the z-score is computed against the permutation distribution (1000 permutations).

This is stricter than CI95-only because it adds an effect-size floor (precursors with AUC = 0.51 + tight CI shouldn't count as "working") and a permutation safeguard.

### 6.4 Multiple-comparisons correction
Six features × three null types × two test regions = 36 tests in the panel. **Bonferroni-corrected α = 0.05 / 36 ≈ 0.00139**, two-sided z-threshold 3.20.

The 3σ headline (one-sided) and the Bonferroni-corrected threshold are both reported. A feature that passes the one-sided 3σ but not the Bonferroni-corrected test is reported as "suggestive, not confirmed."

### 6.5 What we will NOT do
- We will **not** tune Mc per region after seeing AUCs.
- We will **not** swap stations after seeing AUCs.
- We will **not** swap regions if a particular train/test split looks unfavorable.
- We will **not** add features after seeing AUCs.
- We will **not** widen the precursor window beyond 30 days, or shrink it.
- We will **not** drop null types from the headline if one of them gives less favorable AUCs.
- We will **not** exclude specific events from the target panel after seeing per-event contributions to AUC.
- If any of the above seems necessary in retrospect, the deviation is documented and the result is downgraded from "headline" to "exploratory."

---

## 7. Decision tree

After running the locked protocol on the test regions:

| Outcome | Result framing | Expected paper score |
|---------|---------------|----------------------|
| ≥ 1 feature passes 3σ AND Bonferroni-corrected, on the test regions | "Cross-regional precursor signal at AUC > 0.55, test-region 3σ" | 7.0–7.5 (or 7.5–8.0 if a physical mechanism is established) |
| ≥ 1 feature passes 3σ but not Bonferroni | "Suggestive cross-regional signal" | 6.5–7.0 |
| 0 features pass 3σ on test regions | "Tightened cross-regional null-corrected upper bound" | 6.0–6.5 |
| Pipeline failure (insufficient data, station gaps, etc.) | "Methodological contribution; result inconclusive" | 5.0–5.5 |

Rounding: subtract 0.5 across the board if the cold-read subagent score (Round E) lands ≥ 0.5 below the in-context score; this is per `~/CLAUDE.md` cold-read protocol.

---

## 8. Reporting commitments

The paper will report:

1. The pre-reg commit SHA (this file at the time of commit).
2. Per-region Mc, target counts, station availability — as a table.
3. Per-region per-feature per-null AUC + bootstrap CI + permutation z — as a table.
4. The macro-average and pooled AUCs.
5. The test-region AUCs side-by-side with training-region AUCs.
6. **Every protocol deviation**, if any, with explicit acknowledgement that it moves the affected statistic from "headline" to "exploratory."
7. The b-vs-Mc plot per region as a supplementary figure.
8. ZBZ histogram per region as a supplementary figure (the auto-threshold vs the published −5.0 default — sensitivity check).

---

## 9. Computational protocol

- Python 3.9.6 (system), venv per `requirements.txt` at commit `b188ce3`.
- Random seed for null-window sampling and bootstrapping: 42.
- Catalog cache and station metadata cache are persisted to disk per experiment; full re-runs from raw FDSN/ComCat are reproducible.
- All per-region waveform pulls executed with `nice -n 10` per `~/CLAUDE.md`.

---

## 10. Acknowledged limitations of this pre-registration

1. **Mexico and Alaska are tectonically subduction-similar to several training regions.** True cross-tectonic-regime extrapolation would require Iran or Indonesia in the test set; we did not commit to those.
2. **Magnitude type heterogeneity.** ComCat returns Mw, Mb, Md, ML mixed; we use whatever the headline magnitude is. This adds noise to b and Benioff features that we cannot control.
3. **Single station per region** for waveform-derived features. A multi-station median would be more robust but requires substantially more waveform-pull infrastructure.
4. **No explicit ETAS baseline.** Mignan & Broccardo 2019's "trivial-baseline" lesson recommends comparing to a 1-2 parameter baseline. We will add this comparison in the paper but it is not pre-registered as a primary result; treat as exploratory.

---

## 11. Pre-registration commitment

**By committing this file to git, I (Matt Loftus) commit to running the cross-regional evaluation under the protocol specified above, with no post-hoc tuning of features, regions, Mc, stations, windows, nulls, or significance thresholds. Deviations are reported and downgraded.**

After Round D evaluation, this commit's SHA is the canonical reference point for "what was pre-registered" in any later paper draft.

— end of pre-registration —
