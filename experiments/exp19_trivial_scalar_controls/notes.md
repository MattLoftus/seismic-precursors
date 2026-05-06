# exp19 — Four Controls on the exp18 Trivial-Scalar Finding

**Status:** complete 2026-05-02. **Verdict: ALL 5 CONTROLS PASS.** Path A (reframe paper) is justified.

## Motivation

exp18 found that the trivial scalar `log10 Benioff in last 5 days` reaches macro AUC=0.779 across 4 LORO regions (CA, Cas, Tur, IT), beating the elaborate exp14 TLS template apparatus (AUC=0.704). Before reframing the paper around this scalar, run the same controls the cold-read peer reviewers demanded for TLS, plus reviewer 3's dependence-aware null they specifically flagged.

## Method

For 939 windows (139 precursor + 800 Null A across 4 LORO-qualifying regions), recompute the main scalar (`log10 Benioff in last 5 days`) and four control variants:

| Control | What it tests | Pass criterion |
|---|---|---|
| C1 | Shift precursor windows back 30 days, recompute scalar | AUC < 0.60 (collapses to chance) |
| C2 | Shift BOTH precursor and null windows back 30 days (placebo) | AUC ≈ C1 (no temporal confound) |
| C3 | Mask the last 5 days; use log10 Benioff over days 0–25 | AUC < 0.60 (signal IS in last 5d) |
| C4 | Block-permutation null (circular shift of labels within region by random offset) | z > 3 (signal robust to autocorrelation) |
| C5 | Paired bootstrap of ΔAUC (main − C1) on shared resamples | CI95 excludes 0 |

## Results

| Test | Macro AUC | CI95 | z | p |
|---|---|---|---|---|
| **MAIN A2_blog_last5d** | **0.779** | **[0.70, 0.86]** | **+11.00** | **~0** |
| C1 precursor shift back 30d | 0.513 | [0.50, 0.53] | +0.51 | 0.61 |
| C2 both shifted (placebo) | 0.499 | [0.47, 0.53] | −0.04 | 0.97 |
| C3 mask last 5d (days 0–25) | 0.530 | [0.49, 0.56] | +1.09 | 0.27 |
| C4 same point estimate, block null | 0.779 | (same) | **+8.37** | ~0 |
| C5 paired ΔAUC main − C1 | +0.267 | [+0.20, +0.32] | — | ~0 |

All five controls pass. Detailed observations:

**C1 (precursor shift):** AUC = 0.513, essentially chance. The signal is fully localized to days [−5, 0] before mainshock — i.e., the immediate-foreshock period. This matches the exp16 finding (TLS mask control) and is consistent with Trugman & Ross 2019.

**C2 (placebo):** AUC = 0.499. Shifting BOTH classes by 30 days yields chance, ruling out any temporal confound (e.g., catalog quality drift, secular trends). This addresses the methods reviewer's "the window-shift control isn't a true negative control" concern: the placebo IS the proper symmetric control, and it returns chance as expected.

**C3 (foreshock mask):** AUC = 0.530. Removing the last 5 days kills the signal. Combined with C1, this confirms the signal is exclusively foreshock-period-localized. The 30-day "precursor window" is misleading — only the final 5 days carry information.

**C4 (block permutation null):** Reviewer 3 predicted that an autocorrelation-aware null would drop z from ~11 to ~3–4. Actual drop: z_iid = +11.00 → z_block = +8.37. **Reviewer 3's pessimistic prediction is empirically refuted.** Block-permutation z is still very strong (8.4σ), p effectively 0. The dependence inflation exists but is modest (~25% z reduction), not catastrophic.

**C5 (paired ΔAUC):** Formal bootstrap test of main vs C1 on shared resamples. ΔAUC = +0.267, CI95 = [+0.20, +0.32], excludes 0 with high confidence. The "foreshock-localized" claim is now formally tested, not just visually obvious.

## Implications for the Paper

**Decision: Path A — reframe paper around trivial scalar.**

The cross-regional finding is robust under all controls a peer reviewer would request:

> Across CA, Cas, Tur, IT, the log10 Benioff energy released in the final 5 days of a 30-day precursor window discriminates from random aftershock-free windows at macro AUC = 0.779 (CI95 [0.70, 0.86]; z_iid = +11.0; z_block = +8.4 under dependence-aware permutation; ΔAUC vs 30-day-shift placebo = +0.267 CI95 [+0.20, +0.32]).

This is the cross-regional restatement of Trugman & Ross 2019 with the simplest possible scalar. It's stronger than the TLS-based AUC = 0.704 and methodologically cleaner. The TLS apparatus should be moved to an "ablation: template matching does not add information" subsection rather than the headline.

## Score Trajectory

| Stage | Score |
|---|---|
| Initial Round-E2 estimate | 6.5–7.0 |
| Cold-read peer review | 5.5–6.0 |
| After exp18 (TLS dies) | 5.5–6.0 |
| **After exp19 (controls all pass)** | **5.5–6.0**, with higher confidence |

Honest score remains 5.5–6.0. The headline AUC went UP (0.70 → 0.78), the methodology got SIMPLER (no template), and the controls are now bulletproof — but the underlying phenomenon is still the well-established Trugman 2019 foreshock pervasiveness, so the cross-regional confirmation is the contribution. Target: SRL methods note or BSSA short article.

## Next Steps

1. Reframe paper.tex around the trivial-scalar finding (~30 min)
2. Apply the other peer-review fixes:
   - Demote failure modes #1, #2, #4 to a single "execution lessons" paragraph
   - Add per-bin CI95 to M-scaling table (exp17 redo with bootstrap)
   - Drop the "ETAS prediction" attribution; call it "GR-consistent triggered seismicity scaling"
   - Clarify that PRA-2 amendments were post-v1 (not binding pre-reg in the strict sense)
   - Quote z_block alongside z_iid in headline numbers
3. Re-compile + commit + push
