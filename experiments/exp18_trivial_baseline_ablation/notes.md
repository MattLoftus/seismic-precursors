# exp18 — Trivial-Baseline Ablation for the TLS Benioff Finding

**Status:** complete 2026-05-02. Result: **TLS apparatus is worse than trivial scalar.**

## Motivation

All three Round-E2 cold-read peer reviewers (domain, methods, statistics) independently flagged the same concern: the TLS Benioff trajectory scan (exp14, macro AUC=0.704, z=+8.01) may be no more than a slope detector that trivially picks up the late-window event uptick caused by foreshocks. Since Benioff strain is monotonically associated with event count in a sub-window, and the precursor template has a positive slope by construction, a Pearson correlation against it is essentially "is there a final-week event count increment?"

The reviewers requested a head-to-head comparison: **TLS template AUC vs. a trivial scalar baseline computed on the same windows, same catalogs, same LORO splits.**

## Method

For each of 939 windows (139 precursor + 800 Null A across 4 LORO-qualifying regions: California, Cascadia, Turkey, Italy), compute six trivial scalars from the same per-region catalog re-binning that exp14 used:

| Baseline | Definition |
|---|---|
| A1_n_last5d | event count in days 25–30 |
| A2_blog_last5d | log10 Benioff energy in days 25–30 |
| B1_n_total | total event count over the 30-day window |
| B2_blog_total | (intended sum log10 Benioff; bug — see Caveat) |
| C1_n_increment_last | n[5] − n[4] |
| C2_blog_increment_last | b[5] − b[4] |

For each baseline, compute LORO macro AUC via the same 4 train/test splits as exp14. Cross-region bootstrap (n=4 regions, B=1000) gives CI95. Within-region label permutation (B=1000) gives z, p.

## Headline Result

| Baseline | Macro AUC | CI95 | z | p |
|---|---|---|---|---|
| A1_n_last5d | 0.662 | [0.59, 0.74] | +6.5 | 7.1e-11 |
| **A2_blog_last5d** | **0.779** | **[0.70, 0.86]** | **+11.3** | **~0** |
| B1_n_total | 0.573 | [0.51, 0.62] | +2.9 | 4.1e-3 |
| B2_blog_total† | 0.779 | [0.70, 0.86] | +11.4 | ~0 |
| C1_n_increment_last | 0.649 | [0.60, 0.73] | +5.7 | 1.5e-8 |
| C2_blog_increment_last | 0.735 | [0.68, 0.81] | +9.1 | ~0 |
| **TLS exp14 reference** | **0.704** | [0.64, 0.77] | +8.0 | 1.1e-15 |

**Δ (TLS − best baseline) = −0.075.** TLS template correlation is **worse** than the simplest possible scalar (log10 Benioff in last 5 days).

Per-region, the trivial scalar A2 dominates TLS in 4/4 regions:
- California: 0.811 vs TLS 0.733
- Cascadia: 0.905 vs TLS 0.809
- Turkey: 0.691 vs TLS 0.608
- Italy: 0.710 vs TLS 0.665

## Interpretation

The TLS template apparatus normalizes by within-window standard deviation (Pearson correlation), discarding the magnitude information that turns out to be the dominant signal. The actual cross-regional finding is:

> **Precursor windows (the 30 days before an M ≥ 4.5 mainshock) have systematically elevated log10 Benioff energy in the final 5 days compared to randomly-sampled aftershock-free windows in the same regions, at macro AUC=0.779 across CA / Cas / Tur / IT.**

This is just the cross-regional restatement of Trugman & Ross 2019 (foreshock pervasiveness in southern California) — which is what the Round-E2 paper rewrite already framed as the contribution. But:

1. The TLS template apparatus must be **dropped from the central claim**. It does not add information; it actively loses it.
2. The trivial scalar finding (AUC=0.78) is **stronger** than the TLS result (AUC=0.70) but methodologically identical to "count foreshocks." The headline number goes UP, but the methodological novelty goes DOWN.
3. The five-failure-mode catalogue, pre-reg + amendment chain, and per-magnitude scaling work all stand independently of which scalar carries the signal.

## Caveats

† The B2 baseline was intended to be log10(sum of Benioff energies across all 6 sub-windows), but the code accessed `b_traj[:, -1]` which (because `b_traj` is per-sub-window not cumulative) equals `b_traj[:, 5]` — identical to A2. The fact that A2 and B2 produce identical AUCs is an artifact of this bug, not an independent finding. The bug does not change the headline (best trivial baseline is A2 = 0.779 either way).

The C2 baseline (log10 Benioff INCREMENT in last sub-window) at 0.735 is the cleanest "is there a late-window uptick" detector, and it also exceeds TLS. So the conclusion holds whether you measure absolute or incremental Benioff.

## Implications for the Paper

Two viable paths:

**Path A — Reframe (recommended).** Drop TLS as the headline; promote the trivial scalar as the cross-regional confirmation result. Update Methods §3.4 (TLS) → "we attempted a TLS template scan but a single-scalar baseline outperforms it; we report the scalar." Headline AUC becomes 0.78 (with same z and CI structure). The paper becomes shorter, more honest, and the Trugman 2019 framing becomes even more direct. Honest score impact: TLS-novelty subtraction balanced by clarity gain → likely stays in the 5.5–6.0 range.

**Path B — Shelve.** Concede that the only positive finding is a re-derivation of Trugman & Ross 2019 with no methodological novelty (template matching adds nothing), and that the failure-mode catalogue alone isn't enough for a standalone publication. Move on.

Decision pending user.
