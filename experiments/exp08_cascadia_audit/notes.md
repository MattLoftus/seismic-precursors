# exp08 — Cascadia AUC=1.00 forensic audit + b-value computability finding

**Date:** 2026-04-29 (Session 9)
**Status:** **Cascadia AUC=1.00 confirmed artifact.** A bigger structural finding: b-value features are uncomputable in 3 of 4 qualifying regions at PRA-2's Mc values combined with 30-day windows. This is pre-reg failure mode #3.
**Predecessor:** `experiments/exp07_macro_pra2/notes.md`

---

## Verdict on the b-vs-null_A z=+3.41

**Cascadia AUC = 1.00 was computed on 1 precursor × 2 null_A finite b-values** (3 data points total). Separation 0.007 (precursor 0.672 vs null_A max 0.665). Spurious.

Per-region b-finite rates:

| Region | precursor (b_OK / N) | null_A (b_OK / N) | null_B (b_OK / N) |
|--------|----------------------|-------------------|-------------------|
| California | 1 / 37 (3%) | 5 / 200 (2%) | 0 / 200 (0%) |
| Cascadia | 1 / 36 (3%) | 2 / 200 (1%) | 0 / 200 (0%) |
| Italy | 0 / 36 (0%) | 3 / 200 (2%) | 0 / 200 (0%) |
| Turkey | 22 / 30 (73%) | 146 / 200 (73%) | 157 / 200 (78%) |

Three of four qualifying regions have **<5% of windows with computable b**. Only Turkey has real data.

## Why — the structural incompatibility

The `bvalue_with_bootstrap` floor (and `aki_bvalue` directly) requires **N ≥ 30 events at or above Mc** for a stable estimate. At PRA-2's per-region Mc:

| Region | Mc | n_above_mc median per 30-day window | Probability of N ≥ 30 |
|--------|----|------------------------------------|-----------------------|
| California | 3.50 | 4 | ~3% |
| Cascadia | 3.50 | 5 | ~3% |
| Italy | 3.50 | 4 | ~0% |
| Turkey | 3.10 | 76 | ~73% |

At Mc ≥ 3.5 with 30-day windows, you need a typical event rate of ≥1 event/day above Mc to reliably hit N=30. California, Cascadia, and Italy don't reach that density at the PRA-2 Mc. Turkey does.

**This isn't a tuning problem.** The Mc values from the b-vs-Mc plateau check are correctly placed (PRA-2 §1.2). The issue is the interaction:

- 30-day windows × Mc ≥ 3 → N_above_Mc ≈ 4–7 typically
- b-value MLE requires N ≥ 30 → b is NaN almost everywhere

Gulia & Wiemer 2019, where the b-drop precursor signal was originally documented, used **months-long windows and lower Mc** in much higher-rate volcanic settings. The signal does not transfer to short windows at moderate Mc.

## What this means for PRA-2's overall result

**The macro AUCs that mattered:** Benioff total, Benioff curvature, n_above_mc — these are computable in every window (no Mc≥30 requirement). Their macro AUCs were 0.51, 0.51, 0.51 (vs A) and 0.53, 0.51, 0.54 (vs B). **Clean null on the testable features.**

**The macro AUCs that didn't matter:** b, b_drift — they were "computed" but on a degenerate set (1–22 windows), giving meaningless statistics. The headline z=+3.41 was an artifact of bootstrap on 3 finite per-region AUCs.

**Single-region b on Turkey** (the only place b is computable in volume): AUC vs null_A = 0.37 (z=−1.20). Precursor windows have LOWER b. Direction matches Gulia-Wiemer 2019. Single region, sub-3σ — exploratory.

## Pre-reg failure mode #3

This is now a three-failure-mode catalogue for the methods paper:

| # | Failure | Severity | Found at |
|---|---------|----------|----------|
| 1 | Overlap rule [t'-30, t'+60] collides with itself in event-dense regions | High | exp06 (v1) |
| 2 | ComCat is teleseismic for non-US; Mc lands at target threshold for Japan/Chile | High | exp06 (v1) |
| 3 | b-value features require N≥30 events ≥ Mc per window; 30-day windows at Mc≥3 don't provide it in moderate-rate regions | High | exp07 + exp08 (PRA-2) |

Each failure has:
- **Concrete diagnostic evidence** (per-region rejection counts, Mc tables, n_above_mc distributions)
- **A pre-registered fix attempt** (PRA-2 amendments 2 and 1 respectively, none for #3 because not yet known)
- **Empirical data on the fix outcome** (PRA-2 fixed #1 partially, #2 mostly, didn't fix #3)

The methods paper writes itself: pre-reg → diagnostic → amendment → re-test → next failure → ... → honest framing that **CSEP-style pre-registration ported to feature-based cross-regional precursor work has these specific structural failure modes that practitioners must address ex ante**.

## Decision: methods paper

Going forward to Round E (paper writing) is the right call. The structural-finding-catalogue is a substantive contribution. Test-region run on Mexico + Alaska is the last empirical step (so we can report end-to-end pre-reg + amendments + train + test).

Score expectation: **6.5–7.0** as a methods paper at GRL or SRL Statistical Seismology. Possibly 7.0–7.5 if the test-region Turkey-style b-direction-match-Gulia-Wiemer holds up as a flagged exploratory secondary finding.

## Files produced

- `audit.png` — per-region b-distribution histograms, makes the 3-of-4 emptiness visible
- `audit_table.csv` — per (region, kind) NaN audit + n_above_mc quartiles + b summary
- `summary.json` — machine-readable verdict
- `run.py`, `run.log`

## Next session (Session 10)

1. **Test-region run:** apply PRA-2 protocol to Mexico + Alaska. Same expectation: panel partly recovers, b features mostly uncomputable, Benioff/n_above_mc clean null. ~30 min compute.
2. **Outline methods paper:** abstract + section headings. Round E begins.
3. Carryover: Mizrahi 2024 PDF.
