# Novelty check v2 — post-exp14 foreshock-detection finding

**Date:** 2026-04-30 (Session 15)
**Trigger:** exp14 (Session 14) found Benioff strain trajectory passes 3σ + Bonferroni cross-regionally at AUC=0.70/z=+8. Mandatory novelty re-check before claiming this as a positive finding.
**Predecessor:** `papers/novelty_check.md` (Session 1, project-start novelty audit)

---

## Queries executed

- "cross-regional foreshock detection Benioff strain template matching earthquake precursor"
- "'foreshock detection' 'cross-regional' OR 'multiple regions' template matching pre-registered"
- "Trugman Ross 2019 high-resolution foreshock California catalog pre-mainshock Benioff"
- "global earthquake foreshock survey statistical signal Helmstetter Sornette 2003 2010 2020"
- "'Pervasive foreshock activity' Trugman Ross 2019 cross-regional AUC pre-registered"
- "earthquake foreshock detection AUC ROC cross-regional precursor magnitude scaling 2024 2025"

## Closest priors

### Trugman & Ross 2019 GRL — *"Pervasive Foreshock Activity Across Southern California"*
**The most direct prior.** 1.81M-event high-resolution template-matched catalog 2008-2017. Found **72% of M≥4 mainshocks have foreshock activity** elevated above local background. Foreshock sequences range from days to weeks (median 16.6 days).

- **vs ours**: They report fraction-of-mainshocks-with-foreshocks; we report cross-regional template-correlation AUC. They are southern California only; we are 4 regions. They use point-counts; we use Benioff $\Sigma\sqrt{E_J}$ trajectory.
- **Novelty distinction**: same phenomenon (foreshocks exist before most mainshocks), different metric and scope.

### Helmstetter, Sornette, Grasso 2003 JGR — *"Mainshocks are aftershocks of conditional foreshocks"* + *"Foreshocks Explained by Cascades of Triggered Seismicity"*
**The theoretical baseline.** ETAS-derived foreshock properties: power-law acceleration, increased seismicity, modified b-value. Foreshocks emerge from cascade triggering, not requiring distinct precursor physics.

- **vs ours**: ETAS predicts what we measured. The Benioff trajectory's last-sub-window spike IS a power-law foreshock acceleration in M-weighted form.
- **Novelty distinction**: we did not propose a mechanism; we measured a previously-predicted phenomenon under a pre-registered cross-regional protocol.

### Khan et al 2025 JGR — *"Effect of Mainshock Selection, Earthquake Catalog and Definition on Foreshock Rate Estimates in Southern California"*
**Very recent**, methodological work on foreshock rate sensitivity to definition. SCSN catalog only.

### Lippiello 2025 GRL — *"Toward Recognizing the Waveform of Foreshocks"*
**Very recent (2025)**. Waveform-based foreshock recognition. 68 M≥6 events since 2011. Reports foreshocks have characteristic rapid-fluctuation / elevated-energy patterns differing from isolated events.

- **vs ours**: They look at individual-foreshock waveforms; we look at strain-release trajectory in pre-event window.

### Hirose et al 2021 JGR — *"Characteristics of Foreshocks Revealed by an Earthquake Forecasting Method Based on Precursory Swarm Activity"*
JMA-catalog Japan-region foreshock characterization.

### Sci. Reports 2024 — *"Deep learning forecasting of large induced earthquakes via precursory signals"*
**Most direct AUC-comparable prior.** Cross-regional deep-learning forecast on 246 induced earthquakes. Headline AUC = 0.657. Notes: "as more events or stations are added to the data pool, the Q value generally decreases" — i.e., generalization is hard.

- **vs ours**: induced earthquakes (different physics from natural tectonic), AUC=0.66 (lower than our 0.70-0.81), no pre-registration.

## What's novel about our work after this check

Honestly:

1. **The phenomenon (foreshocks exist before mainshocks; their strain release is M-weighted-elevated)** — NOT NOVEL. Trugman 2019, Helmstetter 2003, many others.
2. **Cross-regional AUC quantification of foreshock detection at AUC=0.70 with Bonferroni-corrected significance** — modestly novel. Most prior work is per-region or uses different metrics.
3. **Pre-registered cross-regional protocol with explicit deviations log** — genuinely novel for this question. CSEP did pre-registration for rate forecasts; this is the first feature-based cross-regional pre-reg in our literature search.
4. **Quantitative target-M scaling of cross-regional AUC** (exp17): 0.624 for M=4.5–5.0, 0.813 for M≥5.0 — modestly novel. Helmstetter 2003 ETAS predicts M-scaling; we quantify it cross-regionally.
5. **Failure-mode catalogue (5 modes: overlap rule × density, ComCat × non-US, b-feature × Mc-window, test-region inaccessibility, foreshock-included-in-precursor-window-definition)** — methodologically substantial; reusable knowledge for future cross-regional pre-reg work.

The HONEST framing for the paper is now:

> **Pre-registered cross-regional confirmation of pervasive foreshock activity (Trugman & Ross 2019) with quantitative AUC, M-scaling, and a five-failure-mode methodological catalogue for CSEP-style cross-regional precursor pre-registration.**

NOT:

> "First positive cross-regional precursor finding."

## Score implication

Original ceiling 6.0–8.0; previous Session 14 raised to 7.0–7.5 on the assumption that exp14 was a novel positive. Post-novelty-check, the realistic landing zone is **6.5–7.0**:

- 6.5: methods paper centered on the failure-mode catalogue + foreshock-confirmation
- 7.0: same plus the M-scaling quantification + pre-reg discipline angle
- 7.5: only if cold-read subagent finds the framing very compelling

Avoiding score inflation. The DNA storage lesson (17 experiments before discovering heavy prior work) doesn't apply here — we caught the prior work after exp14, before paper writing. Score adjusts down honestly.

## Recommended paper framing

Move away from "positive finding" language. Instead:

1. Headline: "We pre-registered a CSEP-style cross-regional feature-based precursor evaluation across 4+2 regions, ran it under v1 + PRA-2 amendments, identified five concrete pre-reg failure modes, and **quantitatively confirmed** the Trugman 2019 foreshock-pervasiveness phenomenon cross-regionally at AUC=0.70 with M-scaling 0.62→0.81."
2. The contribution is: **rigor + breadth + quantification + failure-mode documentation**, not a new discovery.
3. Cite Trugman 2019 + Helmstetter 2003 prominently. Position as confirmation + methodological.
