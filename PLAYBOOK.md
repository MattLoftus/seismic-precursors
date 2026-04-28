# Seismic Precursors — Null-Corrected Cross-Regional Earthquake Precursor Search

**Project:** Seismic Precursor Hunt via IRIS FDSN Open Data
**Dir:** `~/workspace/seismic-precursors/`
**Master playbook:** `~/workspace/MASTER_PLAYBOOK.md`
**Status:** INITIALIZED — not yet started
**Created:** 2026-04-23
**Owner:** Matt Loftus / Cedar Loop LLC

---

## 1. What Is This?

A research project that applies the TLS-style matched-filter + null-model methodology from Cedar Loop's exoplanet-search pipeline (which has produced 4 CTOIs submitted to ExoFOP) to 20 years of IRIS FDSN continuous seismic waveform data, searching for statistical precursors to M≥4.5 earthquakes across 6-10 active regions. The field has been burned by false precursor claims (VAN, Parkfield) specifically because null models were absent — an honest null-corrected analysis is both methodologically rare and scientifically valuable.

**Deliverable:** A research paper — either (a) a null-corrected detection of a statistical precursor signal, or (b) a rigorous null-corrected upper bound that sharpens the field's understanding of what's actually detectable. Either outcome is publishable at 7.0+ in the seismology methods literature.

**Target venue:** GRL (Geophysical Research Letters), BSSA (Bulletin of the Seismological Society of America), or SRL (Seismological Research Letters).

**Not a deliverable:** operational earthquake prediction. Precursor detection in specific regions, on specific time scales, with honest statistical bounds — yes. "We predict the next California quake" — no, ever.

---

## 2. Scientific Question

The earthquake-precursor literature has been in a reputational hole since Parkfield (1985-2004 experiment, failed to predict the 2004 event) and VAN (Varotsos-Alexopoulos-Nomicos SES-based method, still disputed). Most precursor claims fail because:
1. No rigorous null-model testing (is the "signal" just noise?).
2. No cross-regional generalization (works in one place, fails everywhere else).
3. Post-hoc tuning to match a specific catalog.

**Question:** Does a **cross-regional, pre-registered, null-corrected** search across 8-10 seismic regions, using features drawn from the established seismology literature (b-value drift, spectral slope change, repeating event rate, Benioff strain accumulation, waveform entropy), produce a detectable statistical precursor to M≥4.5 earthquakes?

**Specifically:** train on 6 regions, test on 2 held-out regions. The precursor must generalize (out-of-sample ROC-AUC > 0.55, beat random at 3σ on bootstrap) to count as a positive result.

---

## 3. Novelty Positioning

### 3.1 What's been done

**Precursor methods (broad literature):**
- Bowman+ 1998, *JGR* 103:24359 — accelerating moment release (AMR) / Benioff strain
- Varotsos, Alexopoulos, Nomicos 1981+ — VAN method (disputed)
- Wyss & Martirosyan 1998 — precursor review
- Nanjo+ 2012 — precursor analysis framework
- Rundle+ 2000+ — Pattern Informatics
- Shcherbakov+ 2015 — review of precursor statistical methods
- Mignan 2014 — precursor review (notably skeptical, important baseline)

**b-value drift:**
- Nanjo, Hirata, Obara+ 2012 — b-value decrease before M9 Tohoku
- Gulia & Wiemer 2019, *Nature* 574:193 — b-value as aftershock hazard indicator (the newer hot area)

**Machine learning for earthquake forecasting:**
- DeVries+ 2018, *Nature* 560:632 — neural net for aftershock location (methodologically criticized, Mignan & Broccardo 2019)
- Mignan & Broccardo 2019, *Nature* 574:193 — ML critique; simpler model beats the NN
- Rouet-Leroy+ 2017, *GRL* — slow slip precursor ML
- Johnson+ 2020 — ML precursor review

**Parkfield (the cautionary tale):**
- Bakun+ 2005 — Parkfield experiment summary; precursors expected, not observed
- Ellsworth 2013, *Science* 342:732 — "Earthquake Prediction" — the balanced review

**VAN method (still disputed):**
- Varotsos+ many — SES precursor claims
- Kagan & Jackson 1996, Mulargia & Gasperini 1996 — critiques

**Nulls + statistical rigor in the field:**
- Swiss Seismological Service has published null-model baselines
- Papadopoulos, Drakopoulos, Lemeille — precursor statistical review
- Stark+ 2020 — importance of pre-registration for earthquake forecasting

### 3.2 What's novel about our approach

> **Updated 2026-04-27 after Session 1 novelty check.** Pre-registration is
> NOT novel in seismology — CSEP/RELM solved it ~15 years ago. Cross-region
> held-out testing is also being done (SafeNet 2025). The methodologically
> closest prior is the MDPI 2025 ionospheric-precursor paper, which applies
> a similar rigor framework to a different physical domain. See
> `papers/novelty_check.md` for the full audit. What remains genuinely novel:

1. **The specific 6-feature stack on equal footing.** Each feature
   (b-value drift, spectral slope drift, repeating-event rate, Benioff strain
   accumulation, waveform entropy, HHT peak drift) is published individually,
   often by different groups. A single pipeline that computes them all under
   a unified Mc / preprocessing / null protocol does not appear in the
   literature. This is where the Mignan & Broccardo 2019 lesson bites: any
   single feature has a trivial baseline; the head-to-head AUC bake-off is
   what survives.
2. **TLS-style matched filter applied to feature time series for precursor
   templates.** Direct port from `~/workspace/exoplanet-search`. The
   Senobari et al. 2024 *JGR* matrix-profile work is the closest matched-
   filter analog in seismology, but applied to event detection (raw
   waveforms), not precursor template scanning across feature space.
3. **A/B/C null triplet on equal footing.** Aftershock-free + random-window
   + colored-noise. The Gulia–Wiemer 2019 line uses random windows only;
   CSEP uses likelihood tests (different statistic); nothing combines all
   three. Tightening the null discipline is the main methodological lever
   we have post-Session-1.
4. **CSEP-style protocol applied to feature-based forecasts.** We adopt the
   CSEP pre-registration / hold-out test / metric-fixed framework but apply
   it to feature-based precursors rather than the rate-based ETAS variants
   that dominate CSEP. This is the bridge contribution.
5. **Honest upper-bound reporting.** Even a null result (no precursor) is
   framed as a tightened upper bound on detectability — which is the
   contribution Mignan 2014 said the field needs.

### 3.3 Mandatory novelty check (Week 1 Day 3 + Week 3 end)

Search arXiv + ADS + Google Scholar + USGS Open-File Reports for:
- "earthquake precursor null model"
- "cross-regional precursor generalization"
- "earthquake precursor pre-registration"
- "TLS earthquake" and "matched filter earthquake precursor"
- "earthquake forecasting hold-out test"
- "b-value drift null model"

Also search SSA / AGU abstracts — much seismology work is conference-first.

If a close prior exists, pivot or narrow. Document in `papers/novelty_check.md`.

---

## 4. Hypothesis

**H0 (primary null):** No feature in our set produces a cross-regional out-of-sample ROC-AUC > 0.55 for predicting M≥4.5 events within a 30-day window, after controlling for aftershock sequences.

**H1 (primary, falsifiable):** At least one feature or feature combination achieves AUC > 0.55 on held-out regions at 3σ bootstrap significance.

**H2 (secondary):** If H1 is confirmed, the feature's predictive horizon (window length for which AUC is maximal) is consistent across regions — indicating physics, not statistical artifact.

**H3 (stretch):** The predictive feature can be linked to a physical mechanism (e.g., slow-slip events on fault asperities, pre-seismic stress redistribution) based on literature.

---

## 5. Methodology

### 5.1 Data pipeline

1. **Catalog ingestion:** Pull ANSS ComCat catalog for M≥2.5 events 2004-2024, 8-10 regions. Store in SQLite.
2. **Target event selection:** M≥4.5 events. After aftershock declustering (Reasenberg 1985 or ETAS, see Zaliapin & Ben-Zion 2013).
3. **Continuous waveform data:** IRIS FDSN web service via ObsPy. Fetch 30-day pre-event windows per target.
4. **Feature extraction (per pre-event window):**
   - **b-value trajectory** via maximum likelihood (Aki 1965)
   - **Spectral slope drift** (amplitude-frequency slope change across window)
   - **Repeating-event rate** (template matching in recent events)
   - **Benioff strain accumulation** (∑√E over window)
   - **Waveform entropy** (Shannon entropy of whitened spectrogram)
   - **Hilbert-Huang peak drift** (EMD decomposition + instantaneous frequency drift)
5. **Null windows:** sample N random 30-day windows NOT preceding any M≥4.5 event by < 60 days. Same features.
6. **TLS-style matched filter:** scan candidate precursor templates across feature space; assess detection statistic distribution under real vs null.
7. **Cross-regional evaluation:**
   - Train on 6 regions
   - Test on 2 held-out
   - Metric: ROC-AUC on held-out regions' M≥4.5 events
8. **Bootstrap:** 1000 resamples for 3σ bound.

### 5.2 Feature details

| Feature | Computation | Literature |
|---------|-------------|------------|
| b-value | MLE: b = log10(e) / (⟨M⟩ - M_c). Sliding 30-day window. | Aki 1965 |
| Spectral slope | Least-squares fit of log-amplitude vs log-frequency 0.5-20 Hz. | Brune 1970 |
| Repeating events | Template correlation > 0.95 on broadband. | Nadeau & McEvilly 1999 |
| Benioff strain | Σ 10^(1.5 M + 4.8) over 30-day window (joules). | Benioff 1951 |
| Waveform entropy | Shannon of whitened spectrogram power, 1s × 5Hz bins. | adapted from Bensen+ 2007 |
| HHT peak drift | First IMF's instantaneous frequency change. | Huang+ 1998 |

### 5.3 Null model details

**Null A — Aftershock-free window:** sample windows that are NOT within 60 days after any M≥4.5 event. Controls for aftershock sequences which might spuriously enhance precursor-like signatures.

**Null B — Random window:** uniform random 30-day windows throughout the 20-year record, excluding only a 7-day post-event aftershock window. Stationary-rate assumption.

**Null C — Colored noise:** synthesize synthetic 30-day records with matched power spectra but randomized phases. Tests whether features reduce to spectral coloring alone.

### 5.4 Pre-registration artifact

Before any evaluation, write `papers/pre_registration.md` committing:
- Feature set (no additions post-hoc)
- Null models (fixed)
- Training regions (6): California, Cascadia, Japan, Chile, Turkey, Italy (default)
- Test regions (2): Mexico, Alaska (default)
- Metric: ROC-AUC on held-out test regions
- Significance threshold: 3σ via bootstrap
- Window: 30-day pre-event

Commit this file to git **before** looking at any results. Hash the commit SHA in the final paper.

---

## 6. Data Sources

### 6.1 Earthquake catalogs

| Source | Access | Coverage |
|--------|--------|----------|
| ANSS ComCat | `libcomcat-python` or REST: `https://earthquake.usgs.gov/fdsnws/event/1/` | Global M≥2.5, 1900-present |
| JMA (Japan) | via IRIS mirror | JMA catalog for Japan |
| EMSC | REST API | European + global |

### 6.2 Continuous waveform data

| Source | Access | Coverage |
|--------|--------|----------|
| IRIS FDSN Web Services | `obspy.clients.fdsn.Client('IRIS')` | Global, all major networks |
| Networks of interest | | |
| — BK (Berkeley Digital) | | California |
| — UW (Pacific NW Seismic Network) | | Cascadia |
| — IU + II (Global Seismic) | | Japan, Chile, global |
| — KO (Kandilli) | | Turkey |
| — IV (INGV) | | Italy |
| — C1 + CX | | Chile |

### 6.3 Station selection

For each region, identify 3-5 broadband stations with near-continuous coverage 2004-2024. Use IRIS DMC's Station Information Service to verify data availability.

### 6.4 Access example

```python
from obspy.clients.fdsn import Client
client = Client("IRIS")

# Continuous waveform fetch
from obspy import UTCDateTime
starttime = UTCDateTime("2024-01-01T00:00:00")
endtime = UTCDateTime("2024-01-02T00:00:00")
st = client.get_waveforms(
    network="BK",
    station="PKD",
    location="*",
    channel="BH?",
    starttime=starttime,
    endtime=endtime
)

# ComCat catalog fetch
from libcomcat.search import search
events = search(
    starttime=starttime,
    endtime=endtime,
    minmagnitude=4.5,
    maxradiuskm=500,
    latitude=35.9,  # Parkfield
    longitude=-120.4
)
```

Cache to `data/iris_cache/` (gitignored — data volume is ~50-200 GB for full region × 20-year span).

---

## 7. Tech Stack

| Layer | Tool | Why |
|-------|------|-----|
| Language | Python 3.9 (`/usr/bin/python3`) | Matt standard |
| Seismic I/O | ObsPy | Community standard for seismology Python |
| Catalogs | libcomcat | USGS-supported |
| Numerical | NumPy, SciPy | Standard |
| Signal processing | SciPy + PyEMD (Hilbert-Huang) | Standard |
| Declustering | custom Reasenberg / Zaliapin & Ben-Zion 2013 | Open-source implementations exist |
| ML baseline | scikit-learn (logistic regression, random forest) | For AUC computation |
| Plotting | matplotlib + cartopy | Geographic maps |
| Database | SQLite (WAL mode) | Same as stock-picks |
| Paper | tectonic | Matt LaTeX |

---

## 8. Infrastructure Reuse

Roughly **40% reuse** from Cedar Loop's portfolio:
- **From `~/workspace/exoplanet-search/`:** TLS matched-filter pattern, nightly cron orchestration, calibration-first methodology, null-model harness, pipeline-per-target layout.
- **From `~/workspace/rmt-neural/`:** bootstrap CI computation, cold-read subagent pattern.
- **From `~/workspace/stock-picks/`:** SQLite + SQLAlchemy patterns for catalog storage, backtest / hold-out test framework.
- **From `~/workspace/quantum-gravity/`:** bootstrap null-model conventions, significance testing framework.

**Net-new:** ObsPy install, continuous seismic waveform handling (preprocessing, declustering), b-value / Benioff / HHT feature extraction (all have open-source reference implementations but need integration), geographic mapping with cartopy.

---

## 9. Experimental Plan

### Round A — Calibration (Week 1-2)
**Goal:** reproduce a published b-value for one well-characterized region.

- **A1:** Install ObsPy + libcomcat. IRIS hello-world.  ✅ 2026-04-27
- **A2:** Pull Parkfield ComCat catalog 2000-2024. Compute b-value via MLE. Compare to published Parkfield b-value (≈0.9 from Bakun+ 2005). **Hard gate:** within 10%.  ✅ 2026-04-28 (exp01 + exp02; declustered Mc=1.50 → b=0.881, 2.1% rel-err)
- **A3:** Novelty check #1. Document.  ✅ 2026-04-27 (`papers/novelty_check.md`)
- **A2b (added in Session 2):** ZBZ 2013 declustering + Mc sensitivity sweep. Headline: declustering matters less than Mc choice for Parkfield; the b-vs-Mc curve diagnoses true completeness better than the Wiemer-Wyss + Woessner default.  ✅ 2026-04-28 (exp02)
- **A4:** Pull BK.PKD continuous waveform for 1 month around 2004 Parkfield M6.0 event. Verify P + S arrivals visually against published picks. Compute waveform entropy + HHT features for sanity.  ✅ 2026-04-28 (exp03; pulled ±30 min not ±1 month — 1 month was scope-creep for a sanity check). Lessons: (i) 60 Hz line noise aliases to Nyquist at fs=40 Hz, mandates pre-EMD bandpass; (ii) EMD's IMF1 in seismic data is amplitude-dominated, sign of IF shift is signal-content-dependent.

### Round B — Feature Implementation (Week 2-3)
**Goal:** all 6 features computed cleanly on Parkfield.

- **B1:** b-value trajectory (sliding window).
- **B2:** Spectral slope drift.
- **B3:** Repeating-event rate (template matching).
- **B4:** Benioff strain accumulation.
- **B5:** Waveform entropy.
- **B6:** HHT peak drift.
- **B7:** Declustering pipeline (Reasenberg / Zaliapin-Ben-Zion).
- **B8:** Pre-registration commit to git with SHA.

### Round C — Single-Region Evaluation (Week 3-4)
**Goal:** establish within-region precursor signals on Parkfield.

- **C1:** For all M≥4.5 Parkfield events 2004-2024, extract 30-day pre-event feature windows.
- **C2:** Sample equal-count null windows (Null A, B, C) per feature.
- **C3:** Per-feature distribution comparison (KS test, t-test, Wasserstein distance) real vs each null.
- **C4:** TLS-style matched filter scan over all feature templates.
- **C5:** Significance via 1000-sample bootstrap.

**Decision point:** if within-region AUC is below 0.52 for all features, the cross-regional case is very unlikely to work. Consider broadening to different regions or scope-narrowing to a specific tectonic setting.

### Round D — Cross-Regional (Week 4-6)
**Goal:** pre-registered train/test protocol.

- **D1:** Extend pipeline to 6 training regions (California, Cascadia, Japan, Chile, Turkey, Italy).
- **D2:** Train logistic regression + random forest on features from training regions.
- **D3:** **UNSEAL TEST REGIONS (Mexico, Alaska).** Compute AUC on pre-specified features.
- **D4:** Bootstrap significance (3σ threshold).
- **D5:** If AUC > 0.55 at 3σ: positive result → Round E.
- **D6:** If AUC ≤ 0.55: negative result → frame as upper bound → still Round E.

### Round E — Paper (Week 6-8)
- **E1:** Outline.
- **E2:** Full draft.
- **E3:** Peer-review simulation (3 subagents: seismologist / statistician / methodology).
- **E4:** Cold-read score. Novelty check #2.
- **E5:** Submit to arXiv (`physics.geo-ph`).
- **E6:** Journal submission (GRL, BSSA, or SRL).

---

## 10. Calibration Gates

| Gate | Metric | Pass action | Fail action |
|------|--------|-------------|-------------|
| **Day 5 b-value Parkfield** | Matches Bakun+ 2005 to 10% at the b-vs-Mc plateau | Proceed | Debug catalog handling; most likely an Mc miscalibration |
|   ↳ Result 2026-04-27 (exp01) | b=0.813±0.007 at Mc=1.35, 9.67% rel-err — PASS by margin only | | Bakun outside CI; exp01 hypothesized aftershock contamination, real cause was Mc too low |
|   ↳ Result 2026-04-28 (exp02) | **declustered, Mc=1.50: b=0.881±0.010, 2.1% rel-err — DEFINITIVE PASS** | | ZBZ declustering tagged 11% as clusters (not bimodal — see notes); Mc choice is the bigger lever than declustering for Parkfield |
| **Day 10 novelty** | No close prior paper | Proceed | Pivot narrower or kill |
| **Week 2 single-waveform sanity** | Visually identifiable P+S arrivals on BK.PKD for known event | Proceed | Debug ObsPy fetch / instrument response removal |
| **Week 3 feature compute** | All 6 features produce sensible distributions on 100 random windows | Proceed | Debug individual feature implementations |
| **Week 4 within-region signal** | At least one feature shows > 2σ real-vs-null gap on Parkfield | Proceed to cross-regional | Negative on Parkfield is a warning sign; decide continue vs narrow |
| **Week 6 cross-regional AUC** | AUC > 0.55 at 3σ OR tight upper bound | Proceed to paper | See Round D decision |
| **Week 7 cold-read** | Score ≥ 6.5 from blind subagent | Submit | Revise methodology |

---

## 11. Scoring Rubric (Paper)

> **Revised 2026-04-27 after Session 1 novelty check (−0.5 across the board).**

| Outcome | Original estimate | Post-novelty estimate | Interpretation |
|---------|------------------|----------------------|----------------|
| All features fail cross-regional; tight upper bound only | 6.5–7.0 | **6.0–6.5** | Honest null result — rare in the field, publishable |
| One feature works cross-regional (AUC > 0.55 at 3σ) | 7.5–8.0 | **7.0–7.5** | Genuine contribution; field-important |
| Multiple features work, physical mechanism links to slow-slip / stress redistribution | 8.5 | **7.5–8.0** | Field-changing; Nature Geoscience shot |
| Cold-read score mandatory at Round E3 before commit | | | |

---

## 12. Target Venue + arXiv Categories

- **Primary:** GRL (fast turnaround, broad audience)
- **Alternate:** BSSA or SRL (more methods-focused, seismology-specific)
- **Shot:** Nature Geoscience (only if cross-regional positive with physical mechanism)

arXiv primary: `physics.geo-ph`
arXiv cross-list: `stat.ML` (Mahoney-endorsed — always available) and optionally `physics.data-an`.

---

## 13. Risks + Mitigations

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| ObsPy install issues | Low | Low | Conda fallback exists; ObsPy is well-supported |
| b-value calibration fails | Low-Medium | High — blocks | Most likely magnitude of completeness M_c mis-specified; use Wiemer & Wyss 2000 method |
| Continuous waveform data access rate-limited | Medium | Medium | IRIS is generous but there are limits; cache aggressively |
| Compute blow-up on feature extraction for full 20-year × 6-region | Medium | Medium | Nice -n 10; overnight cron; pre-compute features per day per station |
| Precursor signal turns out to be an aftershock artifact | High | Medium-High | Aftershock declustering (Reasenberg or ZBZ) is the primary defense; test both |
| Signal is region-specific (works in CA, fails Japan) | High (this is the historical failure mode) | Medium | Cross-regional train/test is designed to catch this; honest reporting of AUC per region |
| Publication ethics (specific forecasts) | Medium | Low-Medium | Do NOT make operational forecasts; frame as statistical methodology paper |
| Pre-registration violation temptation | Medium | High — kills credibility | Hard-commit feature set pre-Round C; hash the pre-reg doc SHA in final paper |

---

## 14. Decision Gates (When to Pivot / Kill)

- **Day 10:** novelty check. Close prior paper → pivot or kill.
- **Week 3:** if within-region signal on Parkfield is < 2σ for ALL features, the cross-regional case is very unlikely. Consider narrowing (focus on one feature family) or scope-shift (aftershock forecasting instead of precursor, which is the Gulia-Wiemer 2019 path).
- **Week 6:** if cross-regional AUC < 0.52 for all features, we have a negative result. Paper IS still publishable but angle shifts from "we found X" to "we place an upper bound on Y." Commit to the frame.
- **Week 7:** cold-read < 6.0 → revise.

---

## 15. Milestones + Timeline

| Milestone | Target | Status |
|-----------|--------|--------|
| Project initialized | 2026-04-23 | ✅ Done |
| ObsPy + libcomcat installed | 2026-04-25 | ✅ Done 2026-04-27 (obspy 1.5.0, usgs-libcomcat) |
| Novelty check #1 | 2026-04-28 | ✅ Done 2026-04-27 (`papers/novelty_check.md`) |
| Parkfield b-value calibrated | 2026-05-02 | ✅ Done 2026-04-28 (exp02: declustered Mc=1.50 b=0.881, 2.1% rel-err) |
| ZBZ 2013 declustering implemented | 2026-05-09 | ✅ Done 2026-04-28 (exp02; pre-reg note added re: log10 η_0=−5.0 default + auto as sensitivity) |
| Round A4: BK.PKD waveform + entropy + HHT sanity | 2026-05-09 | ✅ Done 2026-04-28 (exp03; P+S arrivals visible at predicted times; entropy drops 1.77 bits, IMF1 IF responds 4.01 Hz). Lesson: pre-EMD bandpass at fs/3 to avoid Nyquist line-noise aliasing. |
| All 6 features implemented | 2026-05-14 | Pending |
| Pre-registration committed | 2026-05-16 | Pending |
| Within-region evaluation | 2026-05-23 | Pending |
| Cross-regional evaluation | 2026-06-06 | Pending |
| Paper draft | 2026-06-20 | Pending |
| arXiv submission | 2026-07-01 | Pending |
| GRL submission | 2026-07-04 | Pending |

---

## 16. Setup Commands

```bash
cd ~/workspace/seismic-precursors
/usr/bin/python3 -m venv venv
source venv/bin/activate

pip install numpy scipy matplotlib scikit-learn
pip install obspy            # seismic I/O and processing
pip install usgs-libcomcat   # USGS ComCat catalog access (NOT 'libcomcat' on PyPI)
pip install cartopy          # geographic maps (Session 2+)
pip install emd PyEMD        # Hilbert-Huang (Session 3 — Round B)
pip install tqdm sqlalchemy
pip install pytest ipython jupyter

pip freeze > requirements.txt

# Hello world
python -c "
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
client = Client('IRIS')
inventory = client.get_stations(network='BK', station='PKD')
print(inventory)
"
```

---

## 17. First Session Plan (2-4 hours)

1. **Orient** (15 min): PLAYBOOK, SOUL, CLAUDE, RESEARCH_LEARNINGS (pay attention to null models, cross-system universality lessons, pre-registration lesson #30).
2. **Install ObsPy + libcomcat** (15 min).
3. **Hello-world** (15 min): fetch station metadata for BK.PKD (Parkfield). If successful, the stack is alive.
4. **Novelty check** (60 min): arXiv + ADS + Google Scholar. Document in `papers/novelty_check.md`. Key queries from PLAYBOOK §3.3. **Hard requirement.**
5. **Round A2 calibration** (60-90 min): pull Parkfield M≥2.5 catalog 2000-2024 via libcomcat. Compute b-value via MLE. Target: within 10% of Bakun+ 2005's ≈0.9.
6. **Commit + log** (15 min): push to GitHub if repo initialized.

Do NOT in Session 1: download continuous waveforms, extract features beyond b-value, evaluate any precursor claims.

---

## 18. Project Directory Layout

```
seismic-precursors/
├── PLAYBOOK.md
├── KICKOFF_PROMPT.md
├── README.md
├── requirements.txt
├── .gitignore
├── src/
│   ├── __init__.py
│   ├── data.py                    ← IRIS + ComCat fetching
│   ├── declustering.py            ← Reasenberg / Zaliapin-Ben-Zion
│   ├── features/
│   │   ├── __init__.py
│   │   ├── bvalue.py
│   │   ├── spectral.py
│   │   ├── repeating.py
│   │   ├── benioff.py
│   │   ├── entropy.py
│   │   └── hht.py
│   ├── nulls.py                   ← Null A, B, C
│   ├── tls.py                     ← Matched filter scan
│   ├── evaluation.py              ← Train/test, AUC, bootstrap
│   └── scoring.py                 ← Cold-read subagent
├── data/                          ← gitignored
│   ├── catalog.db                 ← SQLite
│   └── iris_cache/
├── experiments/
│   ├── exp01_parkfield_bvalue/
│   ├── exp02_feature_impl/
│   └── ...
├── notebooks/
├── papers/
│   ├── references/
│   ├── novelty_check.md
│   ├── pre_registration.md        ← commit BEFORE Round D
│   └── draft/
└── figures/
```

---

## 19. References

**Precursor methods:**
- Bowman, Ouillon, Sammis, Sornette, Sornette 1998 — AMR
- Nanjo, Hirata, Obara 2012 — b-value pre-M9
- Gulia & Wiemer 2019 *Nature* — b-value / aftershock
- Rundle, Tiampo, Klein, Sá Martins 2000 — Pattern Informatics

**Precursor skeptics / rigor:**
- Mignan 2014 — precursor review (skeptical)
- Mignan & Broccardo 2019 *Nature* — critique of DeVries+ 2018
- Ellsworth 2013 *Science* — balanced prediction review
- Stark+ 2020 — pre-registration in forecasting

**Parkfield experiment:**
- Bakun, Aagaard, Dost+ 2005 — Parkfield 2004 event summary

**VAN:**
- Varotsos, Alexopoulos, Nomicos 1981+ (many)
- Kagan & Jackson 1996; Mulargia & Gasperini 1996 — critiques

**ML for seismology:**
- DeVries+ 2018 *Nature* — aftershock neural net (criticized)
- Rouet-Leroy+ 2017 *GRL* — slow slip ML
- Johnson+ 2020 — ML review

**Methodology:**
- Aki 1965 — b-value MLE
- Wiemer & Wyss 2000 — M_c estimation
- Reasenberg 1985 — declustering
- Zaliapin & Ben-Zion 2013 — improved declustering
- Bensen+ 2007 — ambient-noise processing (for entropy feature)
- Huang+ 1998 — Hilbert-Huang Transform
- Brune 1970 — spectral slope

**Software:**
- Beyreuther+ 2010 — ObsPy
- `libcomcat` — USGS

---

## 20. Upkeep Rules

Per CLAUDE.md:
- Update PLAYBOOK whenever architecture changes.
- Update `~/workspace/MASTER_PLAYBOOK.md` on status change.
- Update `~/workspace/SCIENTIFIC_SUBMISSIONS.md` on paper status change.
- Update `~/workspace/RESEARCH_LEARNINGS.md` with transferable lessons.
- **CRITICAL:** pre-registration doc in git BEFORE Round D. Hash its SHA in final paper.

---

*Initialized 2026-04-23.*
