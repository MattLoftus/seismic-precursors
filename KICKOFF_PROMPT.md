# Seismic Precursors — Kickoff Prompt

*Paste the content below into a new Claude Code chat to begin this project.*

---

## [Begin paste into new chat]

I'm starting a new research project: **a null-corrected, pre-registered, cross-regional search for statistical precursors to M≥4.5 earthquakes using 20 years of public IRIS FDSN seismic data.** The project is fully scoped in `~/workspace/seismic-precursors/PLAYBOOK.md`. Please read that file in full before doing anything else.

Also read, in this order:
1. `~/SOUL.md` — my ethos
2. `~/CLAUDE.md` — project defaults, scoring rubric, upkeep rules
3. `~/.claude/projects/-Users-Loftus/memory/MEMORY.md` — persistent memory
4. `~/workspace/RESEARCH_LEARNINGS.md` — 90+ cross-project lessons. **Required reading** — particularly the null-model lessons, cross-system universality lessons (#52, #55, #59, #68), pre-registration + score-inflation lessons (#30, #62), and lesson #28 on batch-testing many hypotheses.
5. `~/workspace/exoplanet-search/PLAYBOOK.md` — the sibling pipeline pattern; I'm reusing the TLS + nightly pipeline structure here.
6. `~/workspace/seismic-precursors/PLAYBOOK.md` — this project's full plan

## The thesis in one sentence

Apply my null-corrected, TLS-style matched-filter pipeline (from exoplanet-search, which produced 4 CTOIs submitted to ExoFOP) to 20 years of IRIS FDSN continuous seismic waveforms across 8-10 regions — with pre-registered train/test regions, pre-registered feature set, and honest cross-regional generalization evaluation — searching for statistical precursors to M≥4.5 earthquakes. Either a genuine precursor detection at 3σ cross-region (field-changing) or a tightened null-corrected upper bound (still publishable at 7.0 given the field's historical lack of rigor).

## Today's goals (first session, 2-4 hours)

1. **Orient** — read all files listed above.
2. **Install ObsPy + libcomcat** — `pip install obspy libcomcat cartopy emd scikit-learn` in fresh venv.
3. **Hello-world** — IRIS station metadata fetch for BK.PKD (Parkfield). If that returns a station object, the stack is alive.
4. **Novelty check** — search arXiv + ADS + Google Scholar + USGS Open-File Reports. Queries in PLAYBOOK §3.3. Document in `papers/novelty_check.md`. **If a close prior paper exists, stop and report to me; pivot needed.** The field has a LOT of precursor literature — this check matters.
5. **Round A2 calibration** — pull Parkfield ComCat catalog (M≥2.5, 2000-2024) via libcomcat. Compute b-value via Aki 1965 MLE. Compare to Bakun+ 2005's published Parkfield b-value (≈0.9). **Hard gate:** within 10%. If M_c (magnitude of completeness) needs calibration, use Wiemer & Wyss 2000.
6. **Log** — `experiments/exp01_parkfield_bvalue/notes.md`.

Do **not** in Session 1: pull continuous waveforms, implement features beyond b-value, touch anything about evaluation or TLS matching.

## Critical standards I expect you to follow

- **Pre-registration is non-negotiable.** Before Round D (cross-regional evaluation), commit `papers/pre_registration.md` to git with:
  - Fixed feature set
  - Fixed null models (A, B, C per PLAYBOOK)
  - Fixed training regions (6 default)
  - Fixed test regions (2 default)
  - Fixed metric (ROC-AUC)
  - Fixed significance threshold (3σ via bootstrap)
  - Hash the commit SHA in the paper. This is my defense against p-hacking myself.
- **Nice every CPU-heavy command.** Feature extraction on 20 years × 6 regions is CPU-heavy. `nice -n 10`.
- **Run long jobs in background.** Continuous waveform downloads + feature extraction can run for hours. `run_in_background: true`. Log progress every 1000 events or similar.
- **Novelty check mandatory.** The precursor literature is HUGE and hot. A close prior paper could exist. Before scoring above 5, search the web (per CLAUDE.md).
- **Null model first.** Every precursor claim gets all three nulls (A aftershock-free, B random window, C colored noise). See RESEARCH_LEARNINGS #20, #65, #74.
- **Cold-read scoring.** Any result worth publishing gets a zero-context subagent cold-read first.
- **Calibrate before exploring.** b-value must match Bakun+ 2005 within 10% on Parkfield before proceeding to any features. Hard gate.
- **Cross-system universality test.** Lesson #55: test a formula on the model farthest from the known case, not closest. After training on 6 regions, I want the test regions (Mexico, Alaska) to be TECTONICALLY DISTINCT from the training set, not similar.
- **Score honestly.** A null result here is GENUINELY INTERESTING — the field needs honest nulls. Don't feel pressure to find a signal; feel pressure to find the truth.

## Environment quick-reference (from MEMORY.md)

- System Python: `/usr/bin/python3` (3.9.6). Do not use Anaconda.
- LaTeX: `tectonic`.
- GitHub CLI: `gh` authenticated as MattLoftus.
- Project root: `~/workspace/seismic-precursors/`
- My exoplanet-search pipeline is at `~/workspace/exoplanet-search/` — TLS matched filter pattern, nightly cron orchestration.
- My stock-picks project is at `~/workspace/stock-picks/` — has SQLite + SQLAlchemy + backtest framework that can be adapted for catalog storage + train/test methodology.

## Why this project is worth doing (given I haven't touched seismology)

- My methodology (null models + cold-read + rigor) is genuinely rare in precursor literature.
- Data is fully open via IRIS FDSN.
- 4-CTOI exoplanet pipeline gives me a TLS + pipeline template.
- Even a clean null result (no precursor) is publishable at 6.5-7.0 because the field needs honest nulls.
- Downside-bounded: 6-week time budget, with Week 2 decision gate on calibration failure.

## Questions I expect you to ask me (before coding)

1. Target venue: GRL (broad + fast), BSSA (seismology-specific), SRL (methodology-friendly)? I lean GRL.
2. Are the 6 training regions + 2 test regions acceptable (California, Cascadia, Japan, Chile, Turkey, Italy + Mexico, Alaska)? Or do I want to pre-register something different (e.g., include Iran, Indonesia — more test set diversity)?
3. Should GitHub repo be public from day 1 (seismology norm is public) or wait until arXiv?
4. Peer review simulation (3 subagents: seismologist / statistician / methodology) before submission — confirm yes?
5. Do I want to extend the methodology paper into a SaaS product (SwarmSense, from the V4 analysis — real-time PH on USGS seismic feed) later, or is this purely a research paper?

## When you're done with Session 1, report back

- Novelty check findings (with links — especially Mignan 2014, Mignan & Broccardo 2019, Bakun+ 2005, Gulia & Wiemer 2019 — the field's touchstones)
- Whether ObsPy installed cleanly
- Whether the IRIS hello-world succeeded
- Whether the b-value calibration passed (within 10% of Bakun+ 2005 for Parkfield)
- Any surprises — especially anything that suggests the scope needs narrowing
- Plan for Session 2 (likely: declustering + continuous waveform for one event)

The novelty check AND the b-value calibration are both hard gates. Do not skip either.

## [End paste into new chat]

---

## Notes for Matt (not part of the prompt)

- Precursor literature is VAST and CONTENTIOUS. Budget real time for the novelty search. Mignan's work (2014 review, 2019 Nature critique of DeVries+ 2018) is the most rigorous baseline — read at minimum the titles/abstracts of those before committing.
- Pre-registration is a genuine differentiator here. A lot of precursor claims die in review because reviewers can detect post-hoc tuning. Committing your feature list + null models + test regions BEFORE evaluation, with a git SHA hash, makes that critique inapplicable.
- The most likely outcome is a null result — precursors are hard. Budget for it emotionally. The paper is still good: "We applied the most rigorous pre-registered cross-regional framework to the precursor question and found AUC ≤ 0.55 across all features, tightening the upper bound on detectability by X%."
- If the null happens, consider the SwarmSense product angle from V4 analysis (real-time PH on USGS seismic feed, insurance/mining API at $99/mo) — the infrastructure transfers.
- 60-day post-event window for aftershock exclusion in nulls is a default; if that's unusual for the region's aftershock time-constant, adjust (Japan subduction zones often have longer aftershock tails).
