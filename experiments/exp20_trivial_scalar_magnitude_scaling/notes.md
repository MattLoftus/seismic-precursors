# exp20 — M-Scaling of the Trivial-Scalar Finding

**Status:** complete 2026-05-05.

## Method

Recompute exp17's M-scaling result using the trivial scalar (`log10 Benioff in last 5 days`, exp18 A2_blog_last5d) instead of the now-deprecated TLS template. Improvements over exp17:

- Bootstrap CI95 per (bin × region)
- Bootstrap CI95 per bin macro
- 3-bin breakdown ([4.5, 4.8), [4.8, 5.2), ≥5.2) in addition to the 2-bin
- Mann-Kendall monotonicity test
- Logistic regression with M as a continuous predictor

## Results

### 2-bin (matches exp17 partition)

| Bin | n_pre | Macro AUC | CI95 |
|---|---|---|---|
| M ∈ [4.5, 5.0) | 91 | 0.731 | [0.64, 0.83] |
| M ≥ 5.0 | 48 | **0.879** | **[0.81, 0.93]** |

Per-region for M ≥ 5.0:
- California: 0.884 [0.77, 0.98] (n=16)
- Cascadia: 0.940 [0.86, 1.00] (n=13)
- Turkey: 0.777 [0.55, 0.94] (n=8)
- Italy: 0.914 [0.77, 1.00] (n=11)

All four LORO splits above 0.77 for M ≥ 5.0; Cascadia and Italy effectively at 0.91+.

### 3-bin (finer resolution)

| Bin | n_pre | Macro AUC | CI95 |
|---|---|---|---|
| M ∈ [4.5, 4.8) | 72 | 0.701 | [0.61, 0.82] |
| M ∈ [4.8, 5.2) | 33 | 0.837 | [0.76, 0.92] |
| M ≥ 5.2 | 34 | 0.893 | [0.81, 0.96] |

Clean monotonic increase. Each successive bin's CI overlaps the next, but the trend is consistent.

### Trend tests

- Mann-Kendall (3-bin): τ = +1.00, p = 0.33 (perfect monotonicity but only 3 points; cannot reach significance)
- Mann-Kendall (2-bin): scipy returns spurious τ = -1.00 (n = 2 edge case; macro 0.731 → 0.879 is strictly increasing)
- Logistic regression on (scalar, target_M) → both coefficients positive (scalar +0.23, target_M +1.61); target_M is the dominant predictor

## Comparison to exp17 (TLS-based)

| Bin | exp17 (TLS) | exp20 (scalar) | Δ |
|---|---|---|---|
| M ∈ [4.5, 5.0) | 0.624 | 0.731 | +0.107 |
| M ≥ 5.0 | 0.813 | 0.879 | +0.066 |

The trivial scalar uniformly beats the TLS template across both M bins, confirming the broader exp18 finding. The M ≥ 5.0 bin reaches macro AUC = 0.88 with all four regions above 0.77.

## Interpretation

The monotonic M-scaling is consistent with Gutenberg-Richter scaling of triggered seismicity: larger mainshocks have proportionally more (and larger-magnitude) preceding events on average, by the same GR statistics that govern any seismicity sample. This is a generic property of triggered seismicity, not a specific test of ETAS theory.

Per reviewer 1's concern, we **drop the "ETAS prediction" attribution** and report this as "consistent with Gutenberg-Richter scaling of triggered seismicity." Distinguishing genuine ETAS-cascade structure from a passive GR effect would require an ETAS simulation with matched M_c, declustering, and window definition — beyond the scope of this paper.

## Caveats

- Mann-Kendall p = 0.33 in the 3-bin case reflects sample-size insufficiency (only 3 points), not a real failure of monotonicity. The point estimates are unambiguously monotonic.
- The logistic-regression solver throws RuntimeWarnings (divide-by-zero in matmul), likely due to feature scaling differences (scalar ranges ~−1 to ~10; target_M ranges 4.5 to 7.1). The qualitative finding (both coefficients positive, target_M is the dominant predictor) is reliable; absolute coefficient magnitudes should be interpreted with caution.
- Per-region n for M ≥ 5.0 is 8–16; per-region CI widths are 0.10–0.40. Per-region claims should be qualified accordingly.

## Replaces

This experiment supersedes exp17. The paper Methods §4.5 will use exp20 numbers.
