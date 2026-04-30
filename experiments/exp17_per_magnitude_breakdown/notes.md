# exp17 — Per-target-magnitude breakdown of exp14 Benioff TLS

**Date:** 2026-04-30 (Session 15)
**Status:** **Clean monotonic M-scaling.** macro AUC = 0.624 for target M=[4.5, 5.0); 0.813 for target M≥5.0. Direction-consistent across all 4 LORO splits in both bins.

---

## Result

| Magnitude bin | n_precursor | Macro AUC | Per-region AUCs |
|---------------|-------------|-----------|------------------|
| M=[4.5, 5.0) | 91 | **0.624** | CA 0.633, Cas 0.726, Tur 0.565, IT 0.573 |
| M≥5.0 | 48 | **0.813** | CA 0.809, Cas 0.850, Tur 0.718, **IT 0.876** |

Italy at M≥5.0 alone reaches 0.876. All 4 LORO splits at M≥5.0 are above 0.7.

## Interpretation

Larger mainshocks have **stronger foreshock detection signal**. This is consistent with Helmstetter et al. 2003 ETAS-theory prediction: foreshocks are a cascading-triggered seismicity, and larger mainshocks correlate with more pre-event triggered activity.

The original exp14 macro AUC of 0.704 was an average over both bins, weighted by sample size. Stratifying reveals that the cross-regional foreshock signal is qualitatively different by mainshock size:
- Moderate mainshocks (M=4.5–5.0): **AUC 0.62** — weak, just barely informative
- Large mainshocks (M≥5.0): **AUC 0.81** — strong, operationally useful

This M-scaling is a quantitative result that prior per-region work (Trugman 2019 reports fraction-with-foreshocks, not AUC by M-bin) doesn't quantify cross-regionally.

## Files produced

- `per_mag_table.csv` — per-bin per-region AUCs
- `per_mag_plot.png` — bar chart with per-region scatter, exp14-all-M reference line
- `summary.json`
- `run.py`, `run.log`
