# exp06 — Cross-regional macro-AUC, 6 training regions (Round C end + pre-reg gap finding)

**Date:** 2026-04-28 (Session 7)
**Status:** **PASS** mechanically; **CRITICAL pre-reg gaps surfaced** that require amendment before Round D headline. The 0/10 macro result is uninterpretable as a science finding — only 4 of 6 regions yielded usable AUCs, and only 1 region (California, N_pre=25) had enough precursor count to be statistically meaningful.
**Pre-reg:** `papers/pre_registration.md` at SHA `a4f1c6fa093ecea8bd00077fb8dbcab008c70eaa`
**Bug fixed during session:** empty precursor list crashed downstream AUC code with KeyError; defensive handling now in `src/region_pipeline.py`.

---

## Goal

Apply the locked pre-reg protocol to all 6 training regions, then compute the macro-AUC pool — the first scientifically meaningful AUC test, since per-region power was always too low for AUC ≈ 0.6 effects to clear 3σ.

## Setup

Same params as exp05, applied via `src/region_pipeline.py`:

- `catalog_m_min = 2.5` (preview deviation; pre-reg = 1.5)
- `target_m_min = 4.5`
- `zbz_log_eta_threshold = -5.0`
- `null_buffer_days_before = 30`, `null_buffer_days_after = 60`
- `n_null_a = n_null_b = 200` per region

## Headline (mechanical)

**0/10 features pass 3σ. 0/10 pass Bonferroni-corrected (α = 1.39 × 10⁻³).**

This is **NOT a science finding**. It is a **pre-reg failure mode**: only 4 of 6 regions contributed any AUC, and only 1 region (California) had enough precursor windows for the b-value features to be computable.

## What actually happened, region by region

| Region | Mc | declustered targets | precursor kept | precursor rejected | null_A kept | result |
|--------|----|---------------------|----------------|---------------------|-------------|--------|
| California | 2.75 | 116 | **25** | 91 | 200 | usable |
| Cascadia | 2.95 | 277 | **4** | 273 | **23** | b NaN; benioff usable but N=4 |
| **Japan** | **4.55** | 4,563 | **0** | 4,563 | — | **CRASHED — empty feat_df, KeyError** |
| **Chile** | **4.45** | 3,275 | **0** | 3,275 | — | **CRASHED** |
| Turkey | 3.50 | 389 | **2** | 387 | 74 | b NaN; benioff usable but N=2 |
| Italy | 3.10 | 198 | **19** | 179 | 200 | usable; b may be NaN at Mc=3.10 |

**Two distinct pre-reg gaps:**

### Gap 1 — Overlap rule fails on event-dense regions

The pre-reg's "reject precursor windows that overlap any other target's [t'-30, t'+60] zone" was calibrated against California-scale event density (~5 declustered M≥4.5 per year). At Cascadia / Japan / Chile densities (10–180 per year), the 90-day forbidden zones cover the full timeline:

- **Japan: 4,563 declustered M≥4.5 events / 25 yr** = 182/yr ≈ one every 2 days. With ~200 zones × 90 days = 18,000 zone-days vs 9,131 calendar days, every window is contaminated. 0 windows kept.
- **Chile: 3,275 / 25 yr** ≈ 130/yr. Same problem. 0 windows kept.
- **Cascadia: 277 / 25 yr** ≈ 11/yr. 4 windows survive.

The macro pool is effectively reduced from a 6-region panel to a **4-region partial panel** plus 2 zero-data regions. This drops cross-region statistical power dramatically and breaks the sign test (need 6 regions to detect 6/0 direction consistency at p<0.02).

### Gap 2 — ComCat is incomplete for non-US regions

ComCat is ANSS, US-centric. International data flows through teleseismic feeds (IRIS/IDA/GSN), which have a typical Mc of 4–5. Our b-vs-Mc plateau check faithfully picked the rolloff:

- **Japan Mc = 4.55** — at the M≥4.5 target threshold!
- **Chile Mc = 4.45** — same.
- **Turkey Mc = 3.50** — partial INGV/AFAD ingest.
- **Italy Mc = 3.10** — INGV well-integrated.
- **California Mc = 2.75** — NCSN/SCSN well-integrated.
- **Cascadia Mc = 2.95** — PNSN partially integrated.

Pre-reg §1.1 says "Magnitude type heterogeneity ... we cannot fix it post-hoc," but this is more severe than that line acknowledges: for Japan and Chile, **Mc above target M means b-value features are essentially uncomputable**. The pre-reg's b-value-drift feature requires N≥30 events ≥ Mc per 30-day window; at Mc=4.55, no 30-day pre-event window has 30 events that high.

## What's required for Round D headline

This is where **pre-registration amendment** comes in. The pre-reg explicitly contemplates a "deviations log" mechanism (§6.5: "deviations are documented and downgraded"). For Round D, we need a **Pre-Registration Amendment v2 (PRA-2)** drafted next session that justifies and adopts:

1. **Catalog source per region** (replaces "ANSS ComCat" universally):
   - California, Cascadia: ComCat (ANSS-fed by NCSN/SCSN/PNSN — already complete)
   - Japan: JMA catalog directly (e.g., via `https://www.hinet.bosai.go.jp/` API)
   - Chile: IPOC/CSN catalogs directly
   - Turkey: AFAD/KOERI directly
   - Italy: INGV directly (already partially via ComCat)
   - Mexico: SSN directly (Mexican National Seismological Service)
   - Alaska: AEC directly (Alaska Earthquake Center)
   - Justification: ComCat is incomplete for non-US; mixing local catalogs in advance to a uniform Mc target.

2. **Modify the overlap rule for high-rate regions:**
   - Option A (preferred): scale forbidden zone by per-region target density. For a region with X targets/year, use a forbidden zone of min(60, 365/X) days post-event.
   - Option B: drop the pre-event 30-day forbidden contribution to the rejection (only post-event). Saves precursors at the cost of allowing some window-overlap with another target's pre-window — which, on reflection, isn't actually a leakage path because we're computing features from EVENTS in the catalog, not from the target itself. The pre-window overlap is harmless as long as the target itself isn't in our window (which it isn't, by definition: precursor window is [t-30, t)).
   - Option B is cleaner and probably the right choice. Justification: the [t'-30, t') zone of another target was over-conservative.

3. **Amend the macro-AUC computation:**
   - Require at least N_pre ≥ 8 precursor windows per region (matches the pre-reg's "minimum target count" rule but now applied to KEPT precursor windows, not declustered targets).
   - Drop regions failing this from the macro pool, document, and adjust Bonferroni accordingly.

These three amendments must be drafted and committed BEFORE running Round D, with explicit cross-references back to this exp06 finding so the chain of custody is intact.

## What the partial macro DOES tell us (suggestive only)

| Feature | Macro vs A | Macro vs B | Sign vs A | Sign vs B | Notes |
|---------|------------|------------|-----------|-----------|-------|
| b | 0.514 (1 region) | 0.429 (1 region) | 0/0 | 0/0 | only California computable |
| b-drift | 0.499 (1) | 0.550 (1) | 0/0 | 0/0 | only California |
| log10 Benioff | 0.394 (4) | 0.479 (4) | 1/3 | 2/2 | weak negative trend vs A |
| Benioff curv | 0.544 (4) | 0.513 (4) | 2/2 | 2/2 | indeterminate; California's z=−2.76 didn't replicate in Italy/Cascadia/Turkey |
| n above Mc | 0.464 (4) | 0.491 (4) | 2/2 | 2/2 | indeterminate |

**Notable non-result:** California's z=−2.76 for Benioff curvature in exp05 did **NOT replicate** at the macro level. Italy showed AUC=0.58 (precursor HIGHER, opposite California's direction). The "AMR-direction-reversed" sub-3σ pattern from exp05 was a single-region artifact, not a cross-regional signal.

That's a clean honest update: exp05's leading-candidate hypothesis is now disconfirmed.

## Score implication

The original ceiling expected one of:
- 6.0–6.5 — clean cross-regional null, tight upper bound
- 7.0–7.5 — b-drop or AMR-reversal replicates
- 7.5–8.0 — multi-feature + mechanism

Updated post-exp06:

- **6.5–7.0** — "Pre-Registered Cross-Regional Earthquake Precursor Search: Lessons from a Failed Protocol" — methodological paper documenting why a CSEP-style pre-reg failed when ported from rate forecasts to feature forecasts cross-regionally. **This is genuinely useful contribution if framed honestly.** GRL or SRL Statistical Seismology section.
- **6.0–6.5** — "Tight cross-regional null-corrected upper bound on M≥4.5 precursor detectability" if amended Round D also yields null.
- **6.5–7.5** — possible if PRA-2 amendments give a real macro-AUC > 0.55 at 3σ on California + Italy + amended-rest.

The headline finding of this session is that **pre-reg discipline catches its own failure modes** — that's the kind of meta-observation that strengthens a methods paper.

## Score now

If we shipped today: **6.0**. The work to date is methodologically sound but lacks a cross-regional result. The pre-reg amendments and Round D cross-regional re-run are required to even get to 6.5. The optimistic ceiling is now ~7.0–7.5 unless something surprises us in the amended macro pool.

## Bug fix made during this session

`src/region_pipeline.py` `run_region_pipeline` previously raised `KeyError: 'window_kind'` when `pre_rows` was empty (Japan, Chile cases). Now it gracefully returns a `RegionResult` with NaN-filled AUC tables, and `exp06`'s loop continues. The bug had no effect on the regions that did succeed.

## Files produced

- `catalog_<region>.csv` — per-region catalog cache (gitignored)
- `feature_summary.csv` — all windows, all regions concatenated
- `macro_auc_table.csv` — per-(feature, null) macro statistics
- `macro_auc_plot.png` — bar chart with cross-region bootstrap CI + per-region scatter
- `sign_test_plot.png` — direction-consistency visualization
- `summary.json` — machine-readable
- `run.py`, `run.log`

## Next session (Session 8)

1. **Draft `papers/pre_registration_amendment_v2.md`**, justifying:
   - Per-region catalog source (JMA / IPOC / AFAD / INGV / SSN / AEC for non-US)
   - Overlap rule: drop the 30-day pre-event side; keep only +60 post-event buffer
   - Per-region precursor-count minimum: ≥ 8 KEPT windows (drop region from macro otherwise)
2. Commit PRA-2 to git; record SHA. Same chain-of-custody discipline as v1.
3. Re-run exp06 under amended protocol → exp07_macro_pra2.
4. Honest comparison: did PRA-2 reveal anything new, or did it stay null?
5. Carryover: Mizrahi 2024 PDF.
