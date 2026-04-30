# exp16 — Mask-final-5-days control on Benioff TLS

**Date:** 2026-04-30 (Session 15)
**Status:** **Verdict: ENTIRELY foreshock-driven.** Macro AUC drops from 0.704 (exp14) to 0.513 (z=+8.01 → +0.49) when sub-window 5 is masked. Zero residual 25-day-ahead precursor signal.

---

## Setup

Same LORO TLS scan as exp14 but the trajectory is truncated to the first 5 sub-windows (days 0-25). Sub-window 5 (days 25-30, the foreshock period) is masked out before template construction and scoring.

## Result

| Quantity | exp14 (full 6-point) | exp16 (5-point, days 0-25) |
|----------|----------------------|------------------------------|
| Macro AUC | 0.704 | **0.513** |
| Bootstrap CI95 | [0.637, 0.773] | [0.494, 0.541] |
| Permutation z | +8.01 | **+0.49** |
| Permutation p | 1.1×10⁻¹⁵ | 0.63 |
| Per-region AUCs | CA 0.73, Cas 0.81, Tur 0.61, IT 0.67 | CA 0.51, Cas 0.49, Tur 0.55, IT 0.50 |

**The signal collapses to chance.** The masked-trajectory AUCs are within ±0.05 of 0.5 across all four LORO splits.

## Interpretation

The exp14 8σ pass is **exclusively** driven by sub-window 5 (days 25–30 before target). There is no detectable Benioff-strain template signal in days 0–25. Our pre-reg's "precursor window" definition `[t-30, t)` includes the foreshock period, and that's where 100% of the discriminative signal lives.

This is the strongest possible reviewer-defense control. The honest paper headline becomes:

> Cross-regional template-correlation in the Benioff strain trajectory passes 3σ + Bonferroni at AUC=0.704 / z=+8 — but the signal lives entirely in the final 5 days of the 30-day pre-event window. Masking the final sub-window collapses macro AUC to 0.513 (z=+0.49). The result is operationally cross-regional FORESHOCK detection in days 25-30, not 25-day-ahead long-distance precursor.

## Files produced

- `masked_macro_auc.json` — machine-readable result
- `masked_template_plot.png` — 5-point template visualization (no terminal spike)
- `run.py`, `run.log`
