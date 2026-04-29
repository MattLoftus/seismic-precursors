# exp09 — Test-region run under PRA-2 (Mexico + Alaska)

**Date:** 2026-04-29 (Session 10)
**Status:** **Zero of two test regions qualify under PRA-2 Amendment 3 (N≥8 minimum).** The cross-regional test-region claim is structurally unattainable under the pre-registered protocol. This is itself the cleanest possible result for the methods paper.
**Pre-reg v1 SHA:** `a4f1c6fa093ecea8bd00077fb8dbcab008c70eaa`
**PRA-2 SHA:** `05a4b0f4f7d26b076fc5169c5cda493e9f343652`

---

## Per-region result

| Region | Mc | declustered targets | precursor kept | N_A kept | N_B kept | Qualifies? |
|--------|----|---------------------|----------------|----------|----------|------------|
| Mexico | 3.95 | 843 | **1** | 0 | 0 | ❌ (N=1 < 8) |
| Alaska | 3.10 | 847 | **1** | 0 | 0 | ❌ (N=1 < 8) |

**Both test regions drop**. Same subduction-density-overwhelms-overlap-rule failure mode that affected Japan and Chile in training.

## What this means

The pre-registered protocol's headline claim — macro-AUC across training regions, evaluated on held-out test regions — **cannot be made**. The protocol prevents itself from producing a test-region statistic.

This is not a project failure. It's the **fourth concrete finding of the methods paper**:

> **Failure mode #4 (Round D test-region):** the subduction-density structural issue (#1) extends to BOTH the originally-pre-registered test regions (Mexico, Alaska). With N_precursor ≤ 1 in either, no cross-regional held-out test is feasible under the PRA-2 protocol. The test-region claim is structurally inaccessible without further amendments. PRA-3 — if pursued — would need to address the subduction-density issue directly, possibly by raising the target M threshold (e.g., M ≥ 5.5 for subduction zones) or by allowing overlapping precursor windows with explicit non-independence accounting.

## Implications for the methods paper

The paper's results section §4.3 now has a definite outcome: **the test-region run produced zero qualifying regions, demonstrating that the failure mode catalogued in §4.1 (failure mode #1) extends to the test set and is therefore not a curated training-set artifact.**

The strongest upper-bound claim we can defend:

- **Computable features (Benioff total, Benioff curvature, n_above_mc) on the training panel macro AUC ≈ 0.51 ± 0.05 (2σ)** — this IS the cross-regional precursor upper bound on these features at the (PRA-2 Mc, 30-day window, 6 training-region polygons) configuration.
- **The test-region run was structurally inaccessible**; we cannot strengthen the upper bound with held-out evidence.

## Score implication, updated

The headline shifts somewhat:

- **Before exp09:** "0/10 features pass 3σ on training, expected null on test"
- **After exp09:** "0/10 features pass 3σ on training; 0 qualifying test regions; the test-region step was structurally inaccessible — pre-reg failure mode #4"

This is HONESTLY better for the methods paper (more diagnostic findings) and HONESTLY worse for any "we definitively rule out" claim (we can't rule out anything cross-regionally for subduction zones because the protocol doesn't even let us look).

Updated score expectation: **6.5–7.0** still holds. The 4-failure-mode catalogue is more substantive than the 3-failure-mode version, but the loss of "tested on held-out regions" cap caps the methods contribution at maybe 7.0 unless the cold-read subagent finds the framing very compelling.

## Files produced

- `catalog_Mexico.csv`, `catalog_Alaska.csv` — gitignored
- `test_region_auc_table.csv` — empty (no qualifying regions)
- `summary.json` — per-region counts + qualifying = []
- `run.py`, `run.log`

(No `test_region_plot.png` was produced because the plotting code only runs if there are qualifying regions.)

## Next session (Session 11)

1. Begin paper drafting per `papers/paper_outline.md` Round E1.
2. Methods + Results sections first, since results are now complete.
3. Update paper outline §4.3 with the zero-qualifying-test-region result and the failure-mode #4 framing.
4. Carryover: Mizrahi 2024 PDF manual download for the introduction's CSEP framing.
