# Paper outline (Round E1)

**Working title:** *Pre-Registered Cross-Regional Earthquake Precursor Search: A Four-Failure-Mode Catalogue*

**Alternative:** *When CSEP-Style Pre-Registration Meets Feature-Based Cross-Regional Precursor Work: A Failure-Mode Catalogue*

**Author:** Matthew Loftus (Cedar Loop LLC, Bremerton WA)
**Target venue:** Geophysical Research Letters (primary); Seismological Research Letters Statistical Seismology section (backup)
**arXiv categories:** physics.geo-ph (primary), stat.ML (cross-list)
**Word target:** GRL ~2,500 words main text, 4 figures; SRL up to ~5,000.
**Pre-reg v1 SHA:** `a4f1c6fa093ecea8bd00077fb8dbcab008c70eaa`
**PRA-2 SHA:** `05a4b0f4f7d26b076fc5169c5cda493e9f343652`
**Code + data:** https://github.com/MattLoftus/seismic-precursors (public)

---

## 1. Pitch in one sentence

We pre-registered a CSEP-style cross-regional earthquake-precursor evaluation using six features (b-value, b-drift, Benioff strain total + curvature, n events ≥ Mc, plus three waveform features deferred to Round C waveforms) on six training regions and two test regions, ran it, found four structural failure modes that the field hasn't characterized, and report the lessons.

## 2. Why this is a paper, not a non-result

Most precursor null results don't make it to the literature because the field has historically required a positive claim to publish. But the **specific failure modes we identified are reusable knowledge**:

1. The CSEP overlap-rejection rule scales badly to event-dense subduction zones.
2. ANSS ComCat is teleseismic-only outside the US, making per-region Mc inconsistent.
3. The b-value MLE's N≥30 floor combined with 30-day windows at moderate Mc renders b-feature precursors uncomputable in most regions, defeating Gulia & Wiemer 2019's signal at our window length.
4. The pre-registered test regions (Mexico, Alaska) — both subduction zones — also fail #1, leaving the cross-regional test-region step structurally inaccessible.

Each is documented with diagnostic evidence + a pre-registered amendment + an empirical fix outcome. The chain of custody (v1 SHA → exp06 evidence → PRA-2 SHA → exp07–09 evidence) is publicly verifiable. **This is the kind of methods contribution that prevents the next group from rediscovering the same failure modes.**

## 3. Section structure

### Abstract (~150 words)
- One-sentence question: do feature-based precursors generalize cross-regionally under CSEP-style pre-registration?
- Methods: pre-registered protocol; six regions trained, two test; ISC catalog; ZBZ declustering; A/B null windows; ROC-AUC at 3σ + Bonferroni; PRA-2 amendment after v1 surfaced two failures.
- Result: 0/10 features pass 3σ. Three structural failure modes characterized.
- One-sentence implication: pre-registration discipline catches its own failure modes early; subduction-zone aftershock density and short-window b-value both have specific, reusable mitigations or known-limitation flags for future work.

### 1. Introduction (~400 words)
- Earthquake precursor literature: Bowman 1998 AMR, Gulia & Wiemer 2019 b-drop, VAN, DeVries 2018 (with Mignan & Broccardo 2019 critique). Generally per-region or post-hoc; rare to see cross-regional, pre-registered, null-corrected.
- CSEP/RELM (Mizrahi 2024) solved pre-registration for *rate* forecasts; ionospheric-precursor work (MDPI 2025) ported to a different physical domain. **No equivalent for feature-based seismic precursors.**
- Our contribution: port CSEP-style pre-registration to feature-based seismic precursor evaluation across six tectonic regimes; report what fails and why.
- Cite TLS / matched-filter background (Hippke & Heller 2019; Senobari 2024) as related work; they're not direct competitors because they target event detection, not precursors.

### 2. Pre-registration protocol (~500 words)
- The pre-reg v1 commitments. Cite SHA. Reproduce the locked feature set, region polygons, ZBZ parameters, null-window definitions, significance thresholds.
- Pre-reg amendments mechanism: PRA-2 grounded in exp06 evidence, committed BEFORE re-evaluation. Cite SHA. Three amendments: ISC catalog source, overlap rule shrinks to [t', t'+60], per-region N≥8 minimum.
- Why we report PRA-2 results as the headline rather than v1: chain-of-custody discipline; v1 was a feasibility experiment; PRA-2 is our headline protocol.
- The "will NOT" list (no Mc tuning, no station/region swaps, no feature additions, etc.) — explicit.

### 3. Methods (~600 words)
- Catalog source (ISC), time range, region polygons, station selection.
- Mc estimation: Wiemer-Wyss + Woessner +0.20 + b-vs-Mc plateau check.
- ZBZ 2013 declustering (b=1.0, d_f=1.6, log10(η_0)=−5.0) applied to both catalog and target set.
- Six features: b (Aki MLE), b-drift (3 sub-windows), Benioff total + curvature, n_above_mc.
  - Three waveform features (spectral slope, entropy, HHT IF) defined but deferred to Round C waveform sub-protocol — deferred for compute reasons; flag explicitly.
  - Repeating-event rate excluded from headline per pre-reg v1 §4.1.
- Null windows: A (aftershock-free, ≥30 d before AND ≥60 d after every target), B (random with 7-d post-event exclusion). Null C (colored noise) deferred (waveform-only).
- Per-region ROC-AUC via Mann-Whitney U; cross-region bootstrap (resample regions); permutation z (shuffle within region); macro-AUC across qualifying regions; Bonferroni-corrected α.

### 4. Results (~700 words)

#### 4.1 Per-region panel under PRA-2 (Table 1, Figure 1)
- Mc per region, declustered targets, kept precursor windows, null counts.
- 4/6 training regions qualify (California 37, Cascadia 36, Turkey 30, Italy 36 kept precursor windows). Japan and Chile drop (1 each) by Amendment 3 N≥8 minimum.
- **Failure mode #1 (overlap rule):** v1 overlap zone [t'-30, t'+60] collides with itself when target density > ~10/yr. Specific evidence: Japan 4,563 declustered M≥4.5 → 0 precursor windows kept under v1, 1 under PRA-2 (Amendment 2). Subduction-zone aftershock density is structurally incompatible.
- **Failure mode #2 (catalog × Mc):** ComCat under v1 gave Japan Mc=4.55, Chile Mc=4.45 (≈ target M≥4.5 threshold). ISC under PRA-2 gave Japan 3.50, Chile 3.50 (ISC merges 130+ regional networks). Tradeoff: California Mc 2.75 → 3.50 because ISC has US-coverage gap relative to ComCat.
- **Failure mode #3 (b-feature × Mc-window-incompatibility):** b-value MLE's N≥30 floor combined with Mc=3.5 + 30-day windows yields N_above_mc median = 4–5 in California/Cascadia/Italy; b is uncomputable in 97% of windows. Only Turkey (Mc=3.10, n_above_mc median 76) has b reliably computable. Gulia-Wiemer 2019 signal does not transfer to short windows.

#### 4.2 Macro-AUC on training regions (Figure 2)
- Computable features (Benioff total, Benioff curvature, n_above_mc): macro AUC = 0.51, 0.51, 0.51 vs Null A; 0.53, 0.51, 0.54 vs Null B. **Clean null.**
- Spurious-pass features: b-vs-null_A AUC=0.722 with z=+3.41 fails CI gate (CI95=[0.37, 1.00]) because per-region: California 0.80 (1 finite per side), Cascadia 1.00 (1×2 finite values), Turkey 0.37 (146 vs 22 finite), Italy NaN. Cascadia "1.00" is computed on 3 data points and is the artifact this paper documents.
- Single-region Turkey result: b vs null_A = 0.37 (z=−1.20, sub-3σ). Direction matches Gulia-Wiemer 2019 foreshock prediction. Reported as exploratory.

#### 4.3 Test-region results — failure mode #4 (no figure, table only)
- Mexico (Mc=3.95): 843 declustered targets → 1 precursor window kept (drops below N≥8).
- Alaska (Mc=3.10): 847 declustered targets → 1 precursor window kept (drops below N≥8).
- **Zero of two test regions qualify** under PRA-2 Amendment 3. The cross-regional test-region claim is structurally inaccessible.
- This is failure mode #4: the subduction-density structural issue from training (#1) extends to BOTH pre-registered test regions, demonstrating that #1 is not a training-curation artifact but a regime-level incompatibility between the precursor-window framework and high-rate subduction-zone seismicity.
- Implication: the strongest cross-regional precursor upper-bound we can defend is the training-panel-only macro AUC ≈ 0.51 (Benioff total, Benioff curvature, n_above_mc). Test-region strengthening of this bound is unattainable without PRA-3 amendments addressing subduction-zone density.

### 5. Discussion (~500 words)

#### 5.1 Each failure mode in context
- #1 (overlap): the rule is inherited from clinical-trials risk-set logic but the analogy fails when post-event aftershock cascades are months long. PRA-2 Amendment 2 (drop pre-event side) is principled but only fixes moderate-density regions. For subduction zones, the precursor-window framework needs a different definition (e.g., conditional precursor windows that EXCLUDE recent-aftershock periods individually).
- #2 (ComCat × Mc): naturally addressed by switching to ISC, but ISC has its own US-coverage limitation. A truly clean approach would be per-region local catalogs (JMA, IPOC, AFAD, INGV, SSN, AEC) — this is what national operational forecasting systems do but is not currently a single-source pipeline. PRA-3 if a follow-up project pursues this.
- #3 (b-feature × short-window): Gulia-Wiemer 2019 worked because their windows were months, their Mc was lower, and their setting was high-rate (volcanic). For 30-day windows at moderate Mc, b is structurally uncomputable. If b is to be a per-window feature, the windows or Mc need to change — but those are headline pre-reg parameters.

#### 5.2 What survived
- The Benioff features, n_above_mc, and (in Turkey) b: clean macro nulls. With 4 qualifying training regions and effectively ~139 pooled precursor windows, the upper bound on cross-regional precursor AUC is informative (CI95 ≈ ±0.05 around 0.5). This is the "tight upper bound" framing.

#### 5.3 What this means for the field
- CSEP-style pre-registration is an asset for honest cross-regional work, but **the protocol writers must internalize three concrete failure modes** identified here. We provide the diagnostic evidence so the next group doesn't have to re-discover them.
- The exploratory Turkey b-direction match with Gulia-Wiemer 2019 is suggestive but unconfirmed; flagged as a target for region-specific follow-up at a different (larger) window length.

### 6. Conclusion (~150 words)
- Three structural failure modes catalogued + amendment-grounded narrative.
- Honest null on testable features; tight upper bound.
- One exploratory direction (Turkey b-drop) replicating Gulia-Wiemer 2019 sub-3σ.
- Recommend protocol amendments for any future cross-regional feature-based precursor study.
- Public code + data + commit history at https://github.com/MattLoftus/seismic-precursors

### Acknowledgements
- ANSS ComCat, ISC bulletin, ObsPy + libcomcat.
- Anthropic Claude for protocol drafting + diagnostic + paper-writing assistance (per arXiv co-authorship policy: Claude is acknowledged but not a co-author).

### Data + code availability
- All code at GitHub repo (public).
- All catalog caches regenerable from ISC FDSN endpoints.
- Pre-reg v1 + PRA-2 SHAs cited as canonical commitment points.

### References (target ~25 citations)
- Aki 1965, Wiemer & Wyss 2000, Shi & Bolt 1982 — b-value methodology
- Bowman et al. 1998 — AMR
- Gulia & Wiemer 2019 — b-drop (the closest direct comparison)
- Mignan 2014 — precursor review (skeptical)
- Mignan & Broccardo 2019 — DeVries critique (trivial baseline)
- Bakun et al. 2005 — Parkfield experiment
- Mizrahi et al. 2024 — CSEP review (the framework we ported from)
- MDPI 2025 (Tao et al.) — ionospheric precursor cross-station test (methodological analog)
- Senobari et al. 2024 — matrix profile in seismology (matched-filter context)
- Reasenberg 1985, Zaliapin & Ben-Zion 2013 — declustering
- Hippke & Heller 2019 — TLS source paper (port lineage)
- Bensen et al. 2007 — ambient-noise processing (entropy feature)
- Huang et al. 1998 — HHT
- Brune 1970 — spectral slope
- Beyreuther et al. 2010 — ObsPy

---

## 4. Figures

| # | Title | Source |
|---|-------|--------|
| 1 | Per-region panel: targets → declustering → kept precursor windows under v1 vs PRA-2 | exp07 `pra2_comparison_v1_vs_v2.png` |
| 2 | Macro-AUC across training regions: 5 features × 2 nulls, with bootstrap CI95 | exp07 `pra2_macro_auc_plot.png` |
| 3 | Test-region per-feature AUC: Mexico + Alaska | exp09 `test_region_plot.png` |
| 4 | The Cascadia AUC=1.00 artifact: per-region b-distribution histograms | exp08 `audit.png` |

Plus 1 supplementary table per region with Mc, declustered targets, kept windows, AUC + CI per feature × null.

## 5. Acknowledged caveats (consolidated for the paper)

1. **Catalog M_min=2.5 is a preview deviation** from pre-reg v1's M≥1.5. Justification: ISC over 2000–2024 at M≥1.5 is intractably large. For the paper, we report M≥2.5 results and explicitly flag this as the carried deviation.
2. **Mexico/Alaska test set is tectonically subduction-similar** to several training regions. Acknowledged in pre-reg v1 §2.3.
3. **Single station per region** for waveform-derived features (deferred from headline).
4. **No ETAS baseline.** Mignan & Broccardo 2019's trivial-baseline critique remains an obligation we partially discharge by reporting Benioff/n_events macro AUCs near 0.5; a fitted ETAS model would be a stronger comparison.
5. **Magnitude type heterogeneity** in ISC bulletin (Mw / mb / Ms / ML) — acknowledged.

## 6. What's NOT in scope of this paper

- Operational earthquake forecasting (we don't claim any). This is statistical methodology.
- Test of Gulia-Wiemer 2019 b-drop at their original window length / Mc — the Turkey direction-match is suggestive but underpowered; we recommend a follow-up.
- The waveform-derived features (spectral slope, entropy, HHT IMF1 IF) — implemented in code but not exercised cross-regionally for this paper; deferred to a follow-up Round C waveform sub-protocol.

## 7. Submission timeline (target)

- Session 11: Draft Methods + Results sections, populate Figures 1–4
- Session 12: Draft Introduction + Discussion + Conclusion, peer-review simulation (3 subagents per CLAUDE.md)
- Session 13: Cold-read score subagent (per CLAUDE.md, blind subagent reading only the paper)
- Session 14: Final polish + arXiv submission + GRL submission cover letter
- Total estimated effort: ~4 sessions of work

## 8. Score expectation

Per pre-reg + PRA-2 decision tree:
- Outcome: 0 features pass 3σ on test regions but the methods catalogue is documented → **6.5–7.0**
- Cold-read score adjustment per CLAUDE.md: −0.5 if subagent disagrees by ≥0.5

Realistic expected: **6.5** at GRL or SRL Statistical Seismology with ~3 reviews.
