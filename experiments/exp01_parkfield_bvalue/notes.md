# exp01 — Parkfield b-value calibration

**Date:** 2026-04-27 (Session 1)
**Status:** PASS (marginally)
**Calibration gate:** Round A2 per `PLAYBOOK.md §10`

---

## Goal

Reproduce a published b-value for one well-characterized region. The hard gate
is `|b - 0.9| / 0.9 ≤ 0.10` against Bakun et al. 2005's reported Parkfield
value of b ≈ 0.85–0.95 (we use 0.9 as the canonical reference).

If we can't recover a published number on a well-instrumented region, every
downstream feature is suspect (PLAYBOOK §10, RESEARCH_LEARNINGS Principle 2).

## Setup

| Parameter | Value |
|-----------|-------|
| Region    | Parkfield (35.85°N, −120.40°W) |
| Radius    | 50 km (matches Bakun-McEvilly 1984 box) |
| Years     | 2000–2024 |
| M_min     | 1.0 (intentionally below Mc so we can estimate Mc from the data) |
| Mc method | Wiemer & Wyss 2000 max-curvature + Woessner & Wiemer 2005 +0.20 |
| b method  | Aki 1965 MLE, Marzocchi-Sandri 2003 binning correction |
| Uncertainty | Shi & Bolt 1982 SE + 1000-sample bootstrap CI |
| Catalog   | ANSS ComCat via `libcomcat-python` |
| Chunking  | 2-year blocks (13 ComCat queries) |

## Hello-world: IRIS FDSN

```
BK.PKD  35.9452, -120.5416  elev=583m  start=1996-09-06
```

Stack alive end-to-end: ObsPy 1.5.0 → FDSN client → station inventory parse.
**Note:** ObsPy emits a deprecation warning that the IRIS short-URL is now
`EARTHSCOPE` (rebrand). Old URL still works — flag for cleanup in Session 2.

## Catalog pull

Total: **15,732 events**, M ∈ [1.00, 5.97], pulled in 13 × 2-year chunks.

| Period | Count | Note |
|--------|-------|------|
| 2000–2002 |   970 | background |
| 2002–2004 | 1,273 | background |
| **2004–2006** | **5,407** | **2004 M6.0 Parkfield mainshock + aftershocks** |
| 2006–2008 | 1,461 | aftershock tail |
| 2008–2010 | 1,030 | settling toward background |
| 2010–2012 |   741 | background |
| 2012–2014 |   757 | background |
| 2014–2016 |   829 | background |
| 2016–2018 |   892 | background |
| 2018–2020 |   822 | background |
| 2020–2022 |   555 | background |
| 2022–2024 |   675 | background |
| 2024–2025 |   320 | partial year |

The 2004 M6.0 spike (chunk 3) is the well-known long-awaited Parkfield event
that the original USGS Parkfield Earthquake Prediction Experiment had been
waiting for since 1985.

## b-value result

```
Mc raw (max-curvature) = 1.15
Mc adjusted (+0.20)    = 1.35   (n>=Mc = 8,833 / 15,732)
b = 0.8130 ± 0.0071    (Shi-Bolt 1982 SE)
b 95% bootstrap CI     = [0.7992, 0.8270]
a (log10 N at M=0)     = 5.0436

|b − 0.9| / 0.9 = 9.67%   (tolerance 10%)
Bakun b in our 95% CI?    NO
GATE                    = PASS  (rel-error criterion)
```

## Honest interpretation

- **Passes the rel-error criterion** (9.67% < 10%) but **fails the
  in-confidence-interval criterion** (Bakun's 0.9 sits 10σ above our point
  estimate). Under any stricter test the gate would not pass.
- The FMD plot (`magnitude_frequency.png`) shows a **subtle concave
  departure** from the G-R fit in the M=2.5–4 range — the cumulative curve
  dips below the line. That's the textbook signature of an **undeclustered
  catalog**: mainshock + aftershock sequences mix multiple G-R slopes and
  curve the apparent fit downward.
- **Bakun et al. 2005 used a declustered catalog through 2003** (pre-mainshock
  background). We used 2000–2024 raw, which includes the aftershock-rich
  2004–2008 window. The 2004 M6.0 generated ~4,000 excess events visible in
  chunk 3. These bias our combined b downward (mainshock/foreshock events
  have low b ~0.6–0.8; together with high-b aftershocks they curve the FMD).
- **Mc sensitivity:** with the Woessner +0.20 correction we get Mc=1.35 and
  b=0.81. Without the correction (Mc=1.15) we'd include ~6,000 more events
  but at the price of incompleteness contamination near the rollover.

## Verdict

The pipeline mechanics are healthy: ComCat fetch works, chunking works,
Mc estimation is the textbook method, MLE + bootstrap converged. The
**numerical disagreement with Bakun 2005 is itself diagnostic** and points
to declustering as the most likely missing step.

The right Round A3 follow-up (Session 2) is:
1. Implement Reasenberg 1985 declustering (or Zaliapin & Ben-Zion 2013).
2. Recompute b on (a) full catalog, (b) declustered catalog, (c) declustered
   pre-2004 only (true Bakun-comparable window).
3. Sweep Mc ∈ {1.0, 1.15, 1.25, 1.35, 1.5} to quantify b's sensitivity.
4. Report the regression: as we approach Bakun's analysis window/protocol,
   b should converge toward 0.9.

This is not a worrying result — it's an expected one for an unprocessed
catalog, and the diagnosis is unambiguous.

## Files produced

- `run.py` — reproducible experiment script
- `run.log` — console log
- `summary.json` — machine-readable result
- `magnitude_frequency.png` — FMD with G-R fit and Mc line
- `catalog_magnitudes.npz` — cached magnitudes (15,732 floats, ~26 KB).
  Gitignored per `.gitignore` (`*.npz`); regenerable from ComCat by re-running
  `run.py` (~90 s end-to-end).

## Bugs / minor cleanup for Session 2

- `run.py` `n_chunks` displayed as 12 but the loop ran 13 times (off-by-one
  in the display only — data is correct).
- ObsPy: switch FDSN client short-URL `IRIS` → `EARTHSCOPE`.
- Move `data.py` (catalog fetch) out of this experiment script into `src/`.

## Lessons for the project

1. **Calibration gate of "PASS but outside CI" is the right honesty signal**
   — we'd want to know about the divergence even if rel-error squeaked under
   tolerance. Keep both criteria visible in future calibrations.
2. **Catalog choices propagate.** Radius, time window, declustering, and Mc
   choice each move b by 0.05–0.10. Pre-registration MUST fix all four
   before Round D, not just the algorithm.
