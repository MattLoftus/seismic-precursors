# Novelty check — Session 1

**Date:** 2026-04-27
**Status:** PASS — project remains novel; specific positioning has shifted
**Required by:** `PLAYBOOK.md §3.3`, `MEMORY:feedback_novelty_check`,
                  RESEARCH_LEARNINGS lesson on the DNA-storage near-miss

---

## Queries executed

- `earthquake precursor null model cross-regional generalization pre-registered 2024 2025`
- `matched filter earthquake precursor template feature space seismology cross-region`
- `Mignan Broccardo 2019 Nature DeVries critique earthquake aftershock neural network`
- `Gulia Wiemer 2019 Nature b-value foreshock aftershock traffic light system`
- `earthquake precursor pre-registration prospective forecast experiment Stark 2020 Schorlemmer`
- `Mizrahi 2024 Reviews Geophysics earthquake forecasts current practices CSEP cross-region`
- `"earthquake forecast" cross-region "held-out" generalization CSEP test 2023 2024`
- `earthquake precursor "Benioff strain" "waveform entropy" "Hilbert-Huang" combined feature pipeline`
- `Senobari "matrix profile" seismology 2024 template matching everything`
- `libcomcat USGS python pip install package name PyPI`  *(infrastructure, not novelty)*

## Closest priors and their distance from this project

### 1. CSEP / RELM (Collaboratory for the Study of Earthquake Predictability)
**Touchstone:** Mizrahi et al. 2024 *Reviews of Geophysics* and the
2025 *Sci Data* benchmark database of 10 years of prospective next-day
California forecasts.

- **What they do:** pre-registered, prospective evaluation of earthquake
  forecast models. Models are submitted in advance, tested on future data,
  scored under fixed metrics. Italy, NZ, US, Japan, China, Europe testing
  centers; 442 models under current evaluation across 30-min, 1-day,
  3-month, 1-year, and 5-year horizons.
- **What they don't do:** the models are almost exclusively rate-based
  (ETAS variants, Reasenberg, smoothed seismicity). Feature-based precursor
  models (b-value drift, spectral slope, entropy) are not the workhorse.
- **Implication for our novelty:** **pre-registration itself is NOT novel
  in seismology.** CSEP solved that ~15 years ago. Our pre-reg artifact
  should explicitly cite CSEP and frame our work as a CSEP-style
  prospective protocol applied to a new feature stack. **Score impact:
  −0.5 from the original PLAYBOOK estimate** (we lose the "first to do
  pre-reg" angle; we keep the "first to apply this exact protocol to this
  exact feature stack" angle).

### 2. SafeNet (Sci Reports 2025) — multimodal NN for intermediate-term forecasting
- **What they do:** ML earthquake forecast trained on global catalog
  excluding US, fine-tuned on early US data, tested on later US data —
  i.e., a real cross-region held-out test.
- **What they don't do:** uses rate features (ETAS-flavored), not the b-value
  / spectral / Benioff / entropy / HHT stack. Doesn't run a TLS-style matched
  filter over feature time series. Doesn't test the A/B/C null triplet.
- **Implication:** confirms cross-region held-out testing IS happening in 2025
  (so we lose the "first cross-region test" angle), but the feature space and
  null discipline differ. **Score impact: −0.25** (cross-region test isn't
  the stand-out novelty either; the *combination* with the feature stack and
  null triplet is what's left).

### 3. Gulia & Wiemer 2019 *Nature* — b-value foreshock-traffic-light
- **What they do:** real-time monitoring of b-value drift around large
  events. Uses a 10–30% drop in localized b-value relative to background as
  a foreshock indicator vs. an aftershock-tail indicator.
- **What they don't do:** combine b-value with other features. Don't run
  cross-regional held-out tests at scale. The follow-up SRL 2020 paper
  ("Two Foreshock Sequences Post Gulia and Wiemer 2019") shows the system's
  actual performance is highly sensitive to expert judgment and Monte-Carlo
  ambiguous in many cases — i.e., the b-value-only signal is real but
  marginal.
- **Implication:** **b-value drift as a precursor IS being studied with
  rigor.** Our angle has to be (a) we test b-value alongside 5 other
  features in a head-to-head AUC bake-off so we can quantify b's marginal
  contribution, and (b) we apply the A/B/C null triplet that the
  Gulia–Wiemer line of work has not. **Score impact: neutral** as long as
  we cite GW heavily and frame our work as the rigorous-null companion.

### 4. Mignan & Broccardo 2019 *Nature* — "one neuron beats deep learning"
- **What they do:** showed that DeVries et al. 2018's deep neural network
  for aftershock pattern prediction (AUC=0.85) is reproduced by a 2-parameter
  logistic regression. Lesson: trivial-baseline comparison matters.
- **Implication:** we MUST include trivial baselines (raw rate; log-linear
  ETAS) in our AUC tables, not just compare to chance. Otherwise reviewers
  with this paper in mind will discount the result. **Methodological
  obligation, not a novelty hit.**

### 5. Matrix Profile in Seismology (Senobari et al. 2024 *JGR Solid Earth*)
- **What they do:** "correlate everything with everything" to find new
  earthquake events in continuous seismograms. Demonstrated on Parkfield.
  Detection task, not precursor task.
- **What they don't do:** apply matched-filter logic to feature time series
  for precursor identification. Their templates are seismograms; ours would
  be feature trajectories.
- **Implication:** this is the closest "matched filter in seismology" work.
  Our TLS port to feature-space precursor scanning is genuinely distinct.
  Cite as related-work for matched-filter context. **No novelty hit.**

### 6. MDPI 2025 — ionospheric precursors with cross-station validation
- **What they do:** 38 years × 100+ ionosonde stations. Multi-station ML with
  strict temporal non-overlap, cross-station independence. Best classifier
  weighted F1 = 0.70 on 10% temporal split, degrading to 0.56 on 50% split.
  Reports that prior ionospheric-precursor claims were inflated by
  spatiotemporal autocorrelation artifacts.
- **What they don't do:** seismic features. They use ionospheric TEC.
- **Implication:** **methodologically the closest analog to our protocol —
  but in a different physical domain.** Our paper should explicitly frame
  itself as "the seismic-feature counterpart to the rigor that the
  ionospheric-precursor community has now achieved." Cite prominently.
  **Score impact: −0.25** (we're applying a known methodological framework
  to a new domain, not inventing it).

### 7. CSEP California 10-year benchmark (*Sci Data* 2025)
- A public 10-year prospective benchmark of next-day California forecasts.
  Useful as a calibration target if we want to add a "we beat the
  reference forecast at AUC for M≥4.5" claim later. **Not a novelty hit;
  potentially a useful baseline.**

## Net positioning after the novelty check

What's left of our novelty after these results land:

1. **The specific 6-feature stack** (b-value drift + spectral slope drift +
   repeating-event rate + Benioff strain accumulation + waveform entropy +
   HHT peak drift) computed on equal footing in a single pipeline.
   Each feature is published individually; the head-to-head comparison
   under one protocol is not.
2. **TLS-style matched filter applied to feature time series for precursor
   templates.** Direct port from `~/workspace/exoplanet-search`. Genuinely
   new in the seismic-precursor literature based on these searches.
3. **A/B/C null triplet on equal footing.** Aftershock-free + random-window
   + colored-noise. The Gulia-Wiemer line has not done this; CSEP does
   different (likelihood-based) statistics; nothing we found applies all
   three.
4. **Honest cross-regional negative result framing.** If we get AUC ≤ 0.55
   the paper is "tightened upper bound on cross-regional precursor
   detectability under a rigorous null triplet" — which is a contribution
   the field needs and which the GW / Mignan threads validate.

## Revised score expectation

| Outcome | PLAYBOOK estimate | Post-novelty-check estimate |
|---------|------------------|----------------------------|
| Tight upper bound, all features fail cross-region | 6.5–7.0 | **6.0–6.5** |
| One feature works cross-region (AUC > 0.55 at 3σ) | 7.5–8.0 | **7.0–7.5** |
| Multiple features + physical mechanism | 8.5 | **7.5–8.0** |

The downward revision reflects: (a) CSEP solved pre-registration, (b) cross-
region testing exists in SafeNet 2025, (c) the MDPI 2025 ionospheric paper is
the methodological precedent in an adjacent domain. The 0.5–1.0 point
reduction is across the board.

The project remains worth doing. **Publishable at GRL even at the lower end.**

## Citations to include in any draft

**Must cite:**
- Bakun et al. 2005 (Parkfield experiment summary — calibration target)
- Mignan 2014 (precursor review, skeptical baseline)
- Mignan & Broccardo 2019 (one-neuron critique — trivial-baseline obligation)
- Gulia & Wiemer 2019 (b-value foreshock traffic light)
- Mizrahi et al. 2024 (CSEP review — pre-reg framing)
- CSEP California 10-year benchmark (Sci Data 2025)
- SafeNet (Sci Reports 2025) — cross-region forecasting analog
- MDPI 2025 ionospheric — methodological analog
- Senobari et al. 2024 (matrix profile in seismology) — matched-filter context
- Aki 1965, Wiemer & Wyss 2000, Shi & Bolt 1982 — b-value methodology

**Should cite:**
- Reasenberg 1985, Zaliapin & Ben-Zion 2013 — declustering
- Bowman et al. 1998 — accelerating moment release
- Hippke & Heller 2019 — TLS source paper (port lineage)
- Beyreuther et al. 2010 — ObsPy

## Verdict

**Novelty gate: PASS** with a one-tier downward revision of expected paper
score. Project not killed. Repositioning from "we invented cross-region pre-
registered null-corrected precursor testing" to "we applied the rigor that
adjacent domains have now achieved to the seismic feature stack for the
first time, with a TLS port that is genuinely new."

Re-run this check at end of Round D before paper writing.
