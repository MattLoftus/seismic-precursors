"""exp09_test_region_pra2: PRA-2 protocol on TEST regions (Mexico + Alaska).

The headline cross-regional result per pre-reg v1 §6.2 is the test-region AUC
on Mexico + Alaska after the training-region pool is computed. exp07 ran
training; this runs test. Together they constitute the full Round D under
PRA-2.

Pre-reg v1 SHA:  a4f1c6fa093ecea8bd00077fb8dbcab008c70eaa
PRA-2 SHA:       05a4b0f4f7d26b076fc5169c5cda493e9f343652

Two declared deviations carried forward from exp07 (preview-only):
  - catalog_m_min = 2.5 (pre-reg = 1.5)
  - Null C (colored noise) skipped (catalog-only, no waveforms)

Outputs:
    catalog_<region>.csv             — per region; gitignored
    test_region_auc_table.csv        — per (feature, null) for each test region
    test_region_plot.png             — per-feature AUC bars by test region
    summary.json                     — machine-readable
"""
from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.region_pipeline import (  # noqa: E402
    FEATURE_COLS,
    PipelineParams,
    run_region_pipeline,
)
from src.regions import TEST_REGIONS  # noqa: E402

START = dt.datetime(2000, 1, 1)
END = dt.datetime(2025, 1, 1)
EXP_DIR = Path(__file__).resolve().parent

PARAMS = PipelineParams(
    catalog_source="ISC",                  # PRA-2 amendment 1
    null_buffer_days_before=0,             # PRA-2 amendment 2
    null_buffer_days_after=60,
    min_kept_precursor_per_region=8,       # PRA-2 amendment 3
    catalog_m_min=2.5,                     # preview deviation
)


def _ts() -> str:
    return dt.datetime.now().strftime("%H:%M:%S")


def main() -> int:
    print(f"[{_ts()}] [exp09] start — PRA-2 protocol on TEST regions", flush=True)
    print(f"[{_ts()}] [exp09] regions: {[r.name for r in TEST_REGIONS]}", flush=True)
    print(f"[{_ts()}] [exp09] catalog_source = {PARAMS.catalog_source}", flush=True)
    print(f"[{_ts()}] [exp09] overlap forbidden zone = "
          f"[t', t'+{PARAMS.null_buffer_days_after}d]", flush=True)

    region_results = []
    for region in TEST_REGIONS:
        catalog_cache = EXP_DIR / f"catalog_{region.name}.csv"
        try:
            res = run_region_pipeline(
                region, START, END, PARAMS,
                catalog_cache=catalog_cache,
                log=lambda m: print(f"[{_ts()}] {m}", flush=True),
            )
        except Exception as e:
            print(f"[{_ts()}] [exp09] {region.name} FAILED: {type(e).__name__}: {e}",
                  flush=True)
            continue
        region_results.append(res)
        print(f"[{_ts()}] [exp09] {region.name} done: N_pre={res.n_precursor} "
              f"N_A={res.n_null_a} N_B={res.n_null_b}", flush=True)

    # Apply Amendment 3
    qual_results = [r for r in region_results
                    if r.n_precursor >= PARAMS.min_kept_precursor_per_region]
    dropped = [(r.region_name, r.n_precursor) for r in region_results
               if r.n_precursor < PARAMS.min_kept_precursor_per_region]
    print(f"\n[{_ts()}] [exp09] qualifying test regions "
          f"(N_pre>={PARAMS.min_kept_precursor_per_region}): "
          f"{[r.region_name for r in qual_results]}", flush=True)
    if dropped:
        print(f"[{_ts()}] [exp09] dropped: {dropped}", flush=True)

    if not qual_results:
        print(f"[{_ts()}] [exp09] WARNING: zero qualifying test regions; cross-regional "
              f"test-region claim cannot be made", flush=True)

    # Per-region AUC table
    rows = []
    for r in qual_results:
        for col in FEATURE_COLS:
            for null_kind in ("null_A", "null_B"):
                key = f"{col}__{null_kind}"
                a = r.auc_table[key]
                rows.append({
                    "region": r.region_name, "feature": col, "null": null_kind,
                    "auc": a["auc"], "ci_lo": a["ci_lo"], "ci_hi": a["ci_hi"],
                    "perm_z": a["perm_z"], "perm_p": a["perm_p"],
                    "n_p": a["n_p"], "n_q": a["n_q"],
                })
    auc_df = pd.DataFrame(rows)
    auc_df.to_csv(EXP_DIR / "test_region_auc_table.csv", index=False)
    print(f"[{_ts()}] [persist] test_region_auc_table.csv ({len(auc_df)} rows)", flush=True)

    print(f"\n[{_ts()}] [headline] test-region per-feature AUC table:", flush=True)
    print(f"{'region':<10s} {'feature':<22s} {'null':<7s} "
          f"{'AUC':>6s} {'CI95':>14s} {'z':>6s} {'p':>10s} {'N_p':>4s} {'N_q':>4s}",
          flush=True)
    for _, r in auc_df.iterrows():
        print(f"{r['region']:<10s} {r['feature']:<22s} {r['null']:<7s} "
              f"{r['auc']:>6.3f} [{r['ci_lo']:.2f},{r['ci_hi']:.2f}] "
              f"{r['perm_z']:+6.2f} {r['perm_p']:>10.2e} "
              f"{int(r['n_p']):>4d} {int(r['n_q']):>4d}", flush=True)

    # === Plot ===
    if qual_results:
        fig, axes = plt.subplots(1, 2, figsize=(13, 5), dpi=120, sharey=True)
        x = np.arange(len(FEATURE_COLS))
        width = 0.7 / max(len(qual_results), 1)
        for ax_i, null_kind in enumerate(("null_A", "null_B")):
            ax = axes[ax_i]
            for r_idx, r in enumerate(qual_results):
                aucs = [r.auc_table[f"{c}__{null_kind}"]["auc"] for c in FEATURE_COLS]
                ci_lo = [r.auc_table[f"{c}__{null_kind}"]["ci_lo"] for c in FEATURE_COLS]
                ci_hi = [r.auc_table[f"{c}__{null_kind}"]["ci_hi"] for c in FEATURE_COLS]
                err_lo = np.array(aucs) - np.array(ci_lo)
                err_hi = np.array(ci_hi) - np.array(aucs)
                positions = x + (r_idx - len(qual_results) / 2 + 0.5) * width
                ax.bar(positions, aucs, width=width,
                       color=plt.cm.tab10(r_idx % 10), alpha=0.7,
                       yerr=[err_lo, err_hi], capsize=3,
                       edgecolor="black", linewidth=0.5,
                       label=r.region_name)
            ax.axhline(0.5, color="black", ls="--", lw=0.8)
            ax.set_xticks(x)
            ax.set_xticklabels([c.replace("_", "\n") for c in FEATURE_COLS], fontsize=9)
            ax.set_ylim(0.20, 0.80)
            ax.set_ylabel("AUC")
            ax.set_title(f"vs {null_kind}")
            ax.grid(alpha=0.3, axis="y")
            if ax_i == 1:
                ax.legend(fontsize=9, loc="upper right", title="test region")
        fig.suptitle(f"exp09: PRA-2 test-region AUCs (Mexico + Alaska)", fontsize=11)
        fig.tight_layout()
        fig.savefig(EXP_DIR / "test_region_plot.png")
        plt.close(fig)
        print(f"[{_ts()}] [plot] test_region_plot.png", flush=True)

    summary = {
        "experiment": "exp09_test_region_pra2",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "pre_reg_v1_sha": "a4f1c6fa093ecea8bd00077fb8dbcab008c70eaa",
        "pra2_sha": "05a4b0f4f7d26b076fc5169c5cda493e9f343652",
        "params": {
            "catalog_source": PARAMS.catalog_source,
            "catalog_m_min": PARAMS.catalog_m_min,
            "null_buffer_days_before": PARAMS.null_buffer_days_before,
            "null_buffer_days_after": PARAMS.null_buffer_days_after,
            "min_kept_precursor_per_region": PARAMS.min_kept_precursor_per_region,
        },
        "test_regions_attempted": [r.region_name for r in region_results],
        "qualifying_test_regions": [r.region_name for r in qual_results],
        "dropped_test_regions": dropped,
        "per_region": [
            {
                "region": r.region_name,
                "mc": r.mc, "n_total": r.n_total, "n_background": r.n_background,
                "n_targets_all": r.n_targets_all,
                "n_targets_background": r.n_targets_background,
                "n_precursor": r.n_precursor, "n_precursor_rejected": r.n_precursor_rejected,
                "n_null_a": r.n_null_a, "n_null_b": r.n_null_b,
            }
            for r in region_results
        ],
        "auc_rows": rows,
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[{_ts()}] [persist] summary.json", flush=True)
    print(f"[{_ts()}] [exp09] done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
