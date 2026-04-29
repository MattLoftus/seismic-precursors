# Pre-Registration Amendment v2 (PRA-2)

**Project:** seismic-precursors (Cedar Loop LLC)
**Date committed:** 2026-04-28
**Author:** Matthew Loftus
**Repo:** https://github.com/MattLoftus/seismic-precursors
**This document:** `papers/pre_registration_amendment_v2.md`
**Pre-reg v1:** `papers/pre_registration.md` at SHA **`a4f1c6fa093ecea8bd00077fb8dbcab008c70eaa`**
**Surfacing experiment:** exp06 at commit **`3a896b53e2eb96fc5acf6f0a3aa0d50b39167f73`**
**This amendment SHA:** *(filled in to PLAYBOOK after this commit lands)*

---

## 0. Why this amendment exists

The pre-reg v1 protocol was applied to all 6 training regions in `experiments/exp06_cross_regional_macro/` (commit `3a896b5`). Two structural failures surfaced **before any test-region evaluation was attempted**:

| Failure | Evidence |
|---------|----------|
| **Overlap rule fails on event-dense regions.** Pre-reg §5.1's `[t'-30, t'+60]` rejection zone collides with itself when target density is high. | Japan: 4,563 declustered M≥4.5 / 25 yr → **0** precursor windows kept. Chile: 3,275 → **0**. Cascadia: 277 → 4. Turkey: 389 → 2. Only California (25) and Italy (19) yielded usable counts. |
| **ComCat is teleseismic-only for non-US regions.** Pre-reg §1.1 specifies ANSS ComCat universally; the b-vs-Mc plateau check faithfully picks the rolloff at the teleseismic Mc. | Japan Mc = **4.55**, Chile Mc = **4.45** (≈ M≥4.5 target threshold). For both, b-value features become uncomputable since N≥30 events ≥ Mc per 30-day window is never satisfied. |

Pre-reg v1 §6.5 explicitly contemplates the case where deviations are required: *"if any of the above seems necessary in retrospect, the deviation is documented and the result is downgraded from 'headline' to 'exploratory.'"* This amendment formalizes three deviations — committed BEFORE re-running Round D — to preserve the headline status of the eventual result.

Per CLAUDE.md and pre-reg v1 §11, the amendments are bound to this commit's SHA. The chain of custody is: v1 (`a4f1c6f`) → exp06 evidence (`3a896b5`) → this amendment (`<this commit>`). Any subsequent paper draft cites all three SHAs.

---

## 1. Amendment 1 — Catalog source: ISC global

**Replaces v1 §1.1.**

### What changes
- **Source becomes ISC (International Seismological Centre) bulletin** for every region, accessed via `obspy.clients.fdsn.Client('ISC').get_events(...)`. Replaces the universal ANSS ComCat in v1 §1.1.
- Implemented in `src/data.py` as `fetch_isc_catalog_bbox(...)` (added in PRA-2 implementation commit, alongside the original `fetch_comcat_catalog_bbox`).
- Caching, chunking, and downstream features identical to ComCat.

### What stays the same
- Magnitude-cutoff for catalog ingestion (M ≥ 1.5 for production, M ≥ 2.5 for preview).
- Time range 2000-01-01 to 2024-12-31.
- Per-region Mc estimation via Wiemer-Wyss + Woessner +0.20 + b-vs-Mc plateau check.

### Justification
ISC is the **authoritative global bulletin**. It merges contributions from 130+ regional network operators including JMA (Japan), CSN/IPOC (Chile), AFAD/KOERI (Turkey), INGV (Italy), SSN (Mexico), AEC (Alaska), NCSN/SCSN (California), and PNSN (Cascadia). Mc is therefore region-specific and reflects the local network density rather than the teleseismic backbone:

| Region | ComCat Mc (v1) | Expected ISC Mc (v2) |
|--------|----------------|----------------------|
| California | 2.75 | ~2.0–2.5 (NCSN/SCSN) |
| Cascadia | 2.95 | ~2.0–2.5 (PNSN) |
| **Japan** | **4.55** | ~2.0–3.0 (JMA, ~thousands of stations) |
| **Chile** | **4.45** | ~2.5–3.5 (IPOC + CSN) |
| Turkey | 3.50 | ~2.5–3.5 (AFAD + KOERI) |
| Italy | 3.10 | ~2.0–2.5 (INGV) |
| Mexico | n/a | ~2.5–3.5 (SSN) |
| Alaska | n/a | ~2.5–3.5 (AEC) |

Single-source ISC is preferred over six per-region locals (JMA + IPOC + AFAD + INGV + SSN + AEC) for three reasons:
1. **Uniform query interface** through one `obspy.clients.fdsn.Client('ISC')`, no per-region API dialect.
2. **Uniform processing chain** — magnitude-type heterogeneity is centralised at ISC's merge step, not introduced by us.
3. **Simpler chain of custody** — one source, one mode of failure, one amendment.

### Acknowledged tradeoffs
- ISC bulletin lags by ~6–18 months compared to local catalogs. For 2024 events near the end of our window, ISC may be incomplete (the ISC reviewed bulletin runs about 24 months behind real time; we will use the "ISC PROVISIONAL" + "ISC REVIEWED" combined where available, and document the per-region completeness lag).
- ISC magnitudes are still mixed types (Mw, mb, Ms, ML); pre-reg v1 acknowledges this as an unfixable noise source.
- ISC FDSN endpoint may be slower than ComCat. Acceptable for a one-off Round D pull with caching.

### Falsifiability
If ISC's Mc for Japan or Chile turns out to be ≥ 4.0 (still teleseismic), the amendment fails to address the gap. In that case we drop the affected region(s) from the macro pool with an explicit rejection note (consistent with pre-reg v1 §1.4's minimum-target-count rule, generalized).

---

## 2. Amendment 2 — Overlap rule: drop pre-event side

**Replaces v1 §5.1.**

### What changes
The forbidden zone per other target `t'` becomes **`[t', t' + 60 days]`** (60 days post-event only).

Pre-reg v1 specified `[t' - 30, t' + 60]` (90 days, including the other target's 30-day pre-event window). The pre-event side is dropped.

### What stays the same
- Precursor window definition: `[t - 30 days, t)`.
- Null window definitions in v1 §5.2 (Null A) and §5.3 (Null B). The same updated forbidden zone `[t', t' + 60]` is used in their overlap-rejection rules.

### Justification
The pre-event side `[t' - 30, t')` was over-conservative. Re-examining the failure modes:

- **What contaminates a precursor window?** Events INSIDE the window that are themselves a target's post-mainshock aftershock cascade. That contamination is captured by the post-event side `[t', t' + 60]`.
- **What does the pre-event side protect against?** Two windows coinciding on their pre-event ramps. But the pre-event period of another target is, by hypothesis, *also* a pre-event period for the precursor; if anything its "background" should look LIKE precursor activity, not like null. Over-rejecting these windows removes precursor signal, not contamination.
- **Empirical consequence:** dropping the pre-event side roughly halves the rejection rate (from observed ~75–100% in event-dense regions to a projected ~50–75%). This is a substantial recovery of precursor windows.

Concretely (rough projection from exp06 numbers):
- California: 91 → ~45 rejected, 25 → ~70 kept
- Cascadia: 273 → ~140 rejected, 4 → ~135 kept
- Japan: 4,563 → ~2,300 rejected, 0 → ~2,000+ kept
- Chile: 3,275 → ~1,650 rejected, 0 → ~1,600+ kept
- Turkey: 387 → ~200 rejected, 2 → ~190 kept
- Italy: 179 → ~90 rejected, 19 → ~110 kept

Total kept precursors goes from ~50 to ~4,000+. This is a 100× recovery. Statistical power for the macro AUC becomes adequate.

### Why this isn't post-hoc tuning
The pre-event side was justified in v1 as a defense against precursor-window-overlapping-precursor-window. Re-examining that defense, no leakage path actually requires it — the precursor window contains EVENTS in the catalog (not the target itself or other targets), and overlap of two precursor windows just means the same events contribute to two precursors' features. That's correlation between samples, not contamination. Standard cross-region bootstrap and per-region permutation z handle this correctly.

This is a CORRECTION, not a tuning. The v1 rule was over-restrictive based on a too-loose application of "what could go wrong"; PRA-2 restricts to the actual contamination path (post-mainshock aftershocks).

### Falsifiability
If empirical kept-precursor counts in the amended re-run come in dramatically below the projection (e.g., < 500 across 6 regions), the amendment is judged ineffective and we revisit before the test-region step.

---

## 3. Amendment 3 — Per-region precursor minimum N ≥ 8

**Adds a constraint not specified in v1.**

### What changes
A region is **dropped from the macro pool** if it has fewer than 8 KEPT precursor windows (after declustering and overlap rejection).

### What stays the same
- Pre-reg v1 §1.4 already specified a "minimum 8 declustered targets per region" rule. PRA-2 generalises this to KEPT precursor windows, since the overlap rejection (with either v1 or PRA-2 rule) reduces from declustered targets to actual usable windows.
- Per-region Mc, ZBZ parameters, feature compute — all unchanged.

### Justification
exp06 demonstrated that declustered-target-count is not a sufficient sample-size guarantee. Cascadia had 277 declustered targets but only 4 kept precursor windows; Turkey had 389 but only 2. With the precursor side of the AUC at N ≤ 4, the per-region z-score is dominated by sampling noise, not signal.

Threshold of 8 chosen to:
1. Match v1 §1.4's declustered-target threshold for consistency.
2. Provide enough headroom that bootstrap CI95 width on AUC is < 0.20 (rough rule: SE ≈ 1/√N_pre).
3. Keep the rule strict enough to reject genuinely under-powered regions but lenient enough that 4–5 of 6 training regions typically qualify.

### Bonferroni adjustment
The pre-reg v1 §6.4 Bonferroni correction used 36 tests (5 features × 2 nulls × 2 test regions × ... well, the v1 number was a heuristic). PRA-2 redefines it as:

   **Bonferroni α = 0.05 / (n_features × n_nulls × n_test_regions × n_qualifying_training_regions)**

Computed AFTER the dropping step is applied. If only 4 of 6 training regions qualify and we have 5 features × 2 nulls × 2 test regions, the Bonferroni denominator becomes 5 × 2 × 2 × 4 = 80; α = 6.25 × 10⁻⁴.

### Falsifiability
If fewer than 3 training regions qualify under the N ≥ 8 rule, the macro pool is too narrow to support a "cross-regional" claim regardless of feature performance. In that case, we report training-region-level results without macro framing, consistent with pre-reg v1 §6.5's deviation-downgrade policy.

---

## 4. What's unchanged from pre-reg v1

Everything not explicitly amended above remains in force:

- Time range 2000-01-01 to 2024-12-31.
- Per-region Mc estimation procedure (Wiemer-Wyss + Woessner +0.20 + plateau check).
- ZBZ 2013 declustering with `log10(η_0) = −5.0`, `b = 1.0`, `d_f = 1.6`, `p = q = 0.5`.
- Target declustering applied (M≥4.5 background-only).
- 30-day precursor windows.
- 6 training regions (California, Cascadia, Japan, Chile, Turkey, Italy) and 2 test regions (Mexico, Alaska) with the same lat/lon polygons in v1 §2.
- Per-region primary + 2 backup stations from v1 §3.
- 6-feature panel (b-value, b-drift, Benioff total + curvature, spectral slope, waveform entropy, HHT IMF1 IF) pinned to commit `b188ce3`.
- Repeating-event rate explicitly excluded from headline (deferred).
- Three null models (A, B, C) with overlap-rejection rule using the AMENDED forbidden zone.
- Macro-averaged ROC-AUC across qualifying training regions.
- Significance: bootstrap CI95-lower > 0.5 AND |effect| > 0.05 AND permutation z > 3, with Bonferroni α as redefined in §3.
- "Will NOT" list: no further Mc tuning, no station/region swaps, no feature additions, no window-length changes, no null exclusions, no event exclusions post-hoc.

---

## 5. New decision tree (replaces v1 §7)

After running the amended protocol on the test regions:

| Outcome | Result framing | Expected paper score |
|---------|---------------|----------------------|
| ≥ 1 feature passes 3σ AND Bonferroni-corrected on test regions, on the amended protocol | "Cross-regional precursor signal at AUC > 0.55, test-region 3σ" | 7.0–7.5 (or 7.5–8.0 if mechanism established) |
| ≥ 1 feature passes 3σ but not Bonferroni | "Suggestive cross-regional signal" | 6.5–7.0 |
| 0 features pass 3σ on test regions, but the v1 → PRA-2 transition + amended-protocol null is documented | "Tightened cross-regional null-corrected upper bound, with pre-registered failure-mode catalogue" | 6.5–7.0 (note: HIGHER than v1's 6.0–6.5 null framing because of the methods contribution) |
| Pipeline still failing after amendments | "Methodological contribution; result inconclusive" | 5.5–6.0 |

The "methods paper because pre-reg failed" framing is now baked into the score ceiling. Even a clean negative result under PRA-2 lands at 6.5–7.0 because the pre-reg + amendment + re-run process is itself a contribution worth documenting.

Cold-read score adjustment per `~/CLAUDE.md` still applies (−0.5 if subagent disagrees by ≥ 0.5).

---

## 6. Acknowledged caveats specific to PRA-2

1. **ISC bulletin lag.** ISC PROVISIONAL is fast but ISC REVIEWED runs ~12-24 months behind real time. For 2023–2024 events we may use PROVISIONAL data with the "review-state" flag exposed in the result. Documented as a per-event metadata field.
2. **ISC magnitude heterogeneity.** Same issue as ComCat (Mw / mb / Ms / ML mixed). Not fixed by ISC.
3. **Switching catalog source mid-protocol.** v1 used ComCat for exp01–06; PRA-2 uses ISC for the headline. The two are not directly comparable; we will report the ComCat result as a "source-sensitivity analysis" appendix in the paper, NOT the headline.
4. **Effective sample-size projection in §2 is back-of-envelope.** Empirical rejection rates may differ. If kept-precursor counts come in dramatically lower than projected, that's itself a finding to report.

---

## 7. PRA-2 commitment

**By committing this file to git, I (Matt Loftus) commit to running the amended cross-regional evaluation under the protocol specified here, with no further post-hoc tuning beyond what is explicitly listed in the amendments. Future deviations are documented and downgraded.**

Round D evaluation is now bound to **the union of pre-reg v1 (`a4f1c6f`) and PRA-2 (this commit's SHA)**. Where they conflict, PRA-2 takes precedence; where v1 is silent, v1 governs.

— end of PRA-2 —
