# exp13 — Joint multi-feature classification with LORO CV

**Date:** 2026-04-30 (Session 14)
**Status:** Clean null. **0/2 classifiers pass 3σ.** Pooled multi-feature classification does not generalize cross-regionally any better than per-feature scalar AUC.

---

## Goal

Pre-reg v1 §3.2 + §7 mentioned `scikit-learn (logistic regression, random
forest)` as part of the cross-regional ML evaluation. exp01–12 ran
per-feature scalar AUCs (which all came in null). This experiment trains
a classifier over all 6 computable features simultaneously to test
whether a learned combination separates precursor from Null A windows
where no single feature does.

## Setup

- **6 features:** Benioff total (log10), Benioff curvature, n above Mc,
  spectral slope, waveform entropy, HHT IMF1 IF.
- **b and b-drift excluded** (failure mode #3: uncomputable in 3 of 4
  qualifying regions).
- **Per-region z-score normalization** before pooling: each region's
  features are standardized using the union of precursor + null A within
  that region. This isolates the within-region precursor-vs-background
  residual from per-region absolute scale heterogeneity.
- **Leave-one-region-out (LORO) CV:** 4 splits; for each held-out region,
  train on the other 3 and predict on the held-out region; ROC-AUC on
  the test fold.
- **Two classifiers:** Logistic regression (`class_weight="balanced"`) and
  random forest (200 trees, balanced).
- **Macro AUC** across the 4 LORO splits, with cross-region bootstrap CI95
  and within-fold permutation-z significance.
- **Pre-reg gates:** 3σ = (CI95-lower > 0.5 AND |effect| > 0.05 AND |z| > 3).
  Bonferroni α = 0.05 / (2 classifiers) = 0.025.

## Result

| Classifier | Macro AUC | CI95 | $z$ | $p$ | 3σ | Bonf |
|-----------|-----------|------|-----|-----|------|------|
| Logistic regression | 0.526 | [0.502, 0.551] | +0.97 | 0.33 | n | n |
| Random forest | 0.480 | [0.446, 0.518] | $-$0.71 | 0.48 | n | n |

Per-region:

| Classifier | California | Cascadia | Turkey | Italy |
|-----------|-----------|----------|--------|-------|
| LR | 0.510 | 0.548 | 0.554 | 0.493 |
| RF | 0.435 | 0.489 | 0.539 | 0.457 |

LR shows a slight positive macro AUC (0.526 against 0.5) but it doesn't
clear 3σ or Bonferroni. RF is essentially random.

The verdict: **no learned combination of the 6 features generalizes
cross-regionally as a precursor classifier.** This is consistent with the
per-feature scalar AUC results from exp07–12 — the combination doesn't
carry signal that the components don't.

## Why this strengthens the paper

This is the kind of test reviewers will ask for: "did you try a
multi-feature classifier?" The answer is now yes, and the answer is null.
The methods catalog now includes the joint-classification result.

## Files produced

- `loro_auc_table.csv` — per-classifier per-test-region AUC
- `loro_roc_plot.png` — ROC curves
- `feature_importance.png` — LR coefs + RF importances per LORO split
- `summary.json`
- `run.py`, `run.log`
