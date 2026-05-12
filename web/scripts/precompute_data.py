"""Pre-compute static JSON for the seismic-precursors web app.

Reads from existing experiment artifacts and emits compact JSON that the
React app loads at startup. No live computation in the browser.

Outputs (all under web/public/data/):

  windows.json
    All 939 windows in the 4 qualifying regions with their 6-point
    Benioff and n_above_mc trajectories, plus metadata (region,
    window_kind, target_M, t_start, t_end).

  auc_curves.json
    For each of 6 sliding 5-day sub-window positions (days [0,5),
    [5,10), ..., [25,30)): per-region AUC of "log10 Benioff in this
    sub-window" against {precursor, Null A}. Demonstrates that signal
    lives only in the last sub-window.

  comparison.json
    Per-region AUC for TLS template scan (exp14) vs trivial scalar
    (exp18) side by side. Plus the macro values and CIs.

  m_scaling.json
    2-bin and 3-bin M-scaling tables from exp20 with per-region values.

  failure_modes.json
    Five failure-mode entries with diagnostic numbers.

  meta.json
    Top-level metadata (experiment SHAs, dates, scores).
"""
from __future__ import annotations

import datetime as dt
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
WEB_DATA = ROOT / "web" / "public" / "data"
WEB_DATA.mkdir(parents=True, exist_ok=True)

MC_PER_REGION = {
    "California": 3.50, "Cascadia": 3.50, "Turkey": 3.10, "Italy": 3.50,
}
QUAL_REGIONS = ["California", "Cascadia", "Turkey", "Italy"]
NULL_KIND = "null_A"
N_SUBWINDOWS = 6
SUBWINDOW_DAYS = 5


def _ts():
    return dt.datetime.now().strftime("%H:%M:%S")


def load_catalog(region):
    candidates = [
        ROOT / "experiments" / "exp06_cross_regional_macro" / f"catalog_{region}.csv",
        ROOT / "experiments" / "exp07_macro_pra2" / f"catalog_{region}.csv",
    ]
    for p in candidates:
        if p.is_file():
            df = pd.read_csv(p)
            df["time"] = pd.to_datetime(df["time"], utc=True, format="ISO8601")
            return df
    raise FileNotFoundError(region)


def trajectory(catalog, t_start, t_end, mc):
    """Return (n_arr, b_arr_log) per sub-window."""
    duration = (t_end - t_start).total_seconds()
    edges = [t_start + pd.Timedelta(seconds=duration * k / N_SUBWINDOWS)
             for k in range(N_SUBWINDOWS + 1)]
    n_arr = np.zeros(N_SUBWINDOWS)
    b_arr = np.zeros(N_SUBWINDOWS)
    for k in range(N_SUBWINDOWS):
        mask = ((catalog["time"] >= edges[k]) & (catalog["time"] < edges[k + 1])
                & (catalog["magnitude"] >= mc - 1e-9))
        sub = catalog.loc[mask]
        n_arr[k] = float(len(sub))
        if len(sub):
            energies = 10 ** (1.5 * sub["magnitude"].to_numpy() + 4.8)
            b_arr[k] = float(np.sum(np.sqrt(energies)))
    b_arr_log = np.log10(np.where(b_arr > 0, b_arr, 1.0))
    return n_arr, b_arr_log


def per_region_auc(scores, y):
    mask = np.isfinite(scores)
    if mask.sum() < 4 or (y[mask] == 1).sum() < 2 or (y[mask] == 0).sum() < 2:
        return float("nan")
    return float(roc_auc_score(y[mask], scores[mask]))


def main():
    print(f"[{_ts()}] [precompute] start", flush=True)

    df = pd.read_csv(ROOT / "experiments" / "exp07_macro_pra2" / "feature_summary.csv")
    df["t_start"] = pd.to_datetime(df["t_start"], utc=True, format="ISO8601")
    df["t_end"] = pd.to_datetime(df["t_end"], utc=True, format="ISO8601")
    df = df[df["region"].isin(QUAL_REGIONS) &
            df["window_kind"].isin(["precursor", NULL_KIND])].copy().reset_index(drop=True)
    df["target_M"] = df["target_M"].fillna(-1)
    print(f"[{_ts()}] [precompute] {len(df)} windows loaded", flush=True)

    catalogs = {region: load_catalog(region) for region in QUAL_REGIONS}

    # === windows.json: trajectories + metadata ===
    print(f"[{_ts()}] [precompute] computing trajectories...", flush=True)
    windows = []
    n_traj = np.zeros((len(df), N_SUBWINDOWS))
    b_traj = np.zeros((len(df), N_SUBWINDOWS))
    for i, row in df.iterrows():
        cat = catalogs[row["region"]]
        n, b = trajectory(cat, row["t_start"], row["t_end"], MC_PER_REGION[row["region"]])
        n_traj[i] = n
        b_traj[i] = b
        windows.append({
            "id": int(i),
            "region": row["region"],
            "kind": row["window_kind"],
            "target_M": float(row["target_M"]) if row["target_M"] > 0 else None,
            "t_start": row["t_start"].isoformat(),
            "t_end": row["t_end"].isoformat(),
            "n_traj": [round(x, 3) for x in n.tolist()],
            "b_traj": [round(x, 3) for x in b.tolist()],
        })
    (WEB_DATA / "windows.json").write_text(json.dumps(windows))
    print(f"[{_ts()}] [precompute] windows.json: {len(windows)} entries", flush=True)

    # === auc_curves.json: per-sub-window AUC sweep ===
    print(f"[{_ts()}] [precompute] computing AUC sweep...", flush=True)
    y = (df["window_kind"] == "precursor").astype(int).to_numpy()
    auc_curves = {"subwindow_days": SUBWINDOW_DAYS, "n_subwindows": N_SUBWINDOWS, "regions": {}}
    for region in QUAL_REGIONS:
        region_mask = (df["region"] == region).to_numpy()
        y_r = y[region_mask]
        region_data = {"sliding_blog": [], "sliding_n": [], "n_pre": int(y_r.sum()),
                       "n_null": int((y_r == 0).sum())}
        for k in range(N_SUBWINDOWS):
            scores_b = b_traj[region_mask, k]
            scores_n = n_traj[region_mask, k]
            region_data["sliding_blog"].append(round(per_region_auc(scores_b, y_r), 4))
            region_data["sliding_n"].append(round(per_region_auc(scores_n, y_r), 4))
        # Also: cumulative-total over all sub-windows
        region_data["total_n_auc"] = round(
            per_region_auc(n_traj[region_mask].sum(axis=1), y_r), 4)
        # Total log10 Benioff: log of sum of energies = log of sum of 10^b_traj_log
        # (which is the sum of within-sub-window energies)
        total_b = np.log10(np.maximum((10 ** b_traj[region_mask]).sum(axis=1), 1.0))
        region_data["total_b_auc"] = round(per_region_auc(total_b, y_r), 4)
        auc_curves["regions"][region] = region_data

    # Macro across regions (mean of per-region AUCs)
    macro_sliding_blog = []
    macro_sliding_n = []
    for k in range(N_SUBWINDOWS):
        macro_sliding_blog.append(round(
            float(np.mean([auc_curves["regions"][r]["sliding_blog"][k] for r in QUAL_REGIONS])), 4))
        macro_sliding_n.append(round(
            float(np.mean([auc_curves["regions"][r]["sliding_n"][k] for r in QUAL_REGIONS])), 4))
    auc_curves["macro_sliding_blog"] = macro_sliding_blog
    auc_curves["macro_sliding_n"] = macro_sliding_n
    auc_curves["sub_window_labels"] = [f"days {k*5}–{(k+1)*5}" for k in range(N_SUBWINDOWS)]
    (WEB_DATA / "auc_curves.json").write_text(json.dumps(auc_curves, indent=2))
    print(f"[{_ts()}] [precompute] auc_curves.json: {N_SUBWINDOWS} sub-windows × "
          f"{len(QUAL_REGIONS)} regions", flush=True)

    # === comparison.json: TLS vs trivial scalar per-region ===
    tls_table = pd.read_csv(ROOT / "experiments" / "exp14_tls_feature_scan" /
                             "tls_per_split_table.csv")
    tls_benioff = tls_table[tls_table["feature"] == "benioff_total_traj"]
    tls_per_region = {row["held_out_region"]: round(float(row["auc"]), 4)
                      for _, row in tls_benioff.iterrows()}
    tls_macro = round(float(np.mean(list(tls_per_region.values()))), 4)

    scalar_table = pd.read_csv(ROOT / "experiments" / "exp18_trivial_baseline_ablation" /
                                "baseline_per_split_table.csv")
    scalar_blog = scalar_table[scalar_table["baseline"] == "A2_blog_last5d"]
    scalar_per_region = {row["held_out_region"]: round(float(row["auc"]), 4)
                         for _, row in scalar_blog.iterrows()}
    scalar_macro = round(float(np.mean(list(scalar_per_region.values()))), 4)

    comparison = {
        "tls": {"per_region": tls_per_region, "macro": tls_macro, "ci_lo": 0.640, "ci_hi": 0.770,
                "z_iid": 8.01, "label": "TLS template correlation (exp14)"},
        "scalar": {"per_region": scalar_per_region, "macro": scalar_macro,
                   "ci_lo": 0.701, "ci_hi": 0.858, "z_iid": 11.0, "z_block": 8.37,
                   "label": "Trivial scalar: log₁₀ Benioff in last 5d (exp18)"},
        "delta": {r: round(scalar_per_region[r] - tls_per_region[r], 4)
                  for r in QUAL_REGIONS},
        "macro_delta": round(scalar_macro - tls_macro, 4),
    }
    (WEB_DATA / "comparison.json").write_text(json.dumps(comparison, indent=2))
    print(f"[{_ts()}] [precompute] comparison.json: TLS {tls_macro} vs scalar {scalar_macro}",
          flush=True)

    # === m_scaling.json: from exp20 ===
    m_macro = pd.read_csv(ROOT / "experiments" / "exp20_trivial_scalar_magnitude_scaling" /
                           "m_scaling_macro_table.csv")
    m_per_region = pd.read_csv(ROOT / "experiments" / "exp20_trivial_scalar_magnitude_scaling" /
                                "m_scaling_table.csv")
    m_data = {"bins_2": [], "bins_3": []}
    for _, row in m_macro.iterrows():
        target = m_data["bins_2"] if row["lo"] in [4.5, 5.0] and row["hi"] in [5.0, 99.0] \
            else m_data["bins_3"]
        per_region = {}
        for _, r in m_per_region[m_per_region["bin_label"] == row["bin_label"]].iterrows():
            per_region[r["region"]] = {
                "auc": round(float(r["auc"]), 4) if not pd.isna(r["auc"]) else None,
                "ci_lo": round(float(r["ci_lo"]), 4) if not pd.isna(r["ci_lo"]) else None,
                "ci_hi": round(float(r["ci_hi"]), 4) if not pd.isna(r["ci_hi"]) else None,
                "n_pre": int(r["n_pre"]),
            }
        target.append({
            "label": row["bin_label"],
            "lo": float(row["lo"]),
            "hi": float(row["hi"]),
            "n_pre": int(row["n_pre_total"]),
            "macro_auc": round(float(row["macro_auc"]), 4),
            "ci_lo": round(float(row["ci_lo"]), 4),
            "ci_hi": round(float(row["ci_hi"]), 4),
            "per_region": per_region,
        })
    (WEB_DATA / "m_scaling.json").write_text(json.dumps(m_data, indent=2))
    print(f"[{_ts()}] [precompute] m_scaling.json: 2-bin + 3-bin", flush=True)

    # === failure_modes.json: hand-written from paper ===
    failure_modes = [
        {
            "n": 1, "title": "Overlap rule × event density",
            "summary": "v1's [t'−30, t'+60] forbidden zone eliminates nearly all precursor windows in event-dense subduction zones.",
            "evidence": {
                "Japan": "0 / 4563 declustered targets retain windows",
                "Chile": "0 / 3275 retain",
                "Cascadia (v1)": "4 / 277 retain (PRA-2: 36)",
                "Turkey (v1)": "2 / 389 retain (PRA-2: 30)",
            },
            "mitigation": "PRA-2 amendment shrunk zone to [t', t'+60]; recovered Cascadia + Turkey but Japan/Chile remain inaccessible.",
        },
        {
            "n": 2, "title": "ANSS ComCat × non-US Mc",
            "summary": "Under ANSS ComCat (teleseismic for non-US regions), Japan Mc=4.55 and Chile Mc=4.45 sit at the M≥4.5 target threshold, making b-features structurally uncomputable.",
            "evidence": {
                "Japan Mc (v1)": "4.55 (at target threshold)",
                "Chile Mc (v1)": "4.45 (at target threshold)",
                "California Mc shift": "2.75 → 3.50 (ISC US-coverage gap)",
            },
            "mitigation": "PRA-2 switched to ISC global bulletin; Japan/Chile Mc dropped to 3.5.",
        },
        {
            "n": 3, "title": "b-feature × Mc-window incompatibility",
            "summary": "The Aki MLE N≥30 floor combined with 30-day windows at moderate Mc renders b-features uncomputable in 3 of 4 qualifying regions.",
            "evidence": {
                "California": "median n events ≥ Mc = 4–5 per window; 3% of windows pass N≥30",
                "Cascadia": "median 4–5 events; 3% pass",
                "Italy": "median 4–5 events; 3% pass",
                "Turkey": "median 76 events; 73% pass",
                "Cascadia b-AUC=1.00 artifact": "computed on 1 precursor × 2 null finite values",
            },
            "mitigation": "Structural — the Gulia-Wiemer 2019 b-drop signal does not transfer to short windows at moderate Mc.",
        },
        {
            "n": 4, "title": "Test-region inaccessibility",
            "summary": "Both pre-registered test regions (Mexico, Alaska) reduce to one kept precursor window each under PRA-2 — far below the N≥8 minimum. This is failure mode #1 generalized to held-out regions.",
            "evidence": {
                "Mexico": "1 / 843 declustered targets retain windows",
                "Alaska": "1 / 847 retain",
            },
            "mitigation": "None within pre-reg scope. Cross-regional test-region strengthening is structurally inaccessible without further amendments.",
        },
        {
            "n": 5, "title": "Precursor window includes foreshock period",
            "summary": "The pre-reg defines the precursor window as [t−30, t), including the immediate-foreshock period. Any signal detected by this window is foreshock detection, not 25-day-ahead long-distance precursor.",
            "evidence": {
                "Trivial scalar main": "0.779 (macro AUC, last 5d Benioff)",
                "Shifted back 30 days": "0.513 (collapses to chance)",
                "Mask last 5d (use 0–25d)": "0.530 (collapses)",
                "Placebo symmetric shift": "0.499 (no temporal confound)",
            },
            "mitigation": "Three post-hoc controls confirm the signal is exclusively foreshock-period-localized. Future pre-reg should commit explicitly on whether 'precursor window' includes the final 5–7 days.",
        },
    ]
    (WEB_DATA / "failure_modes.json").write_text(json.dumps(failure_modes, indent=2))
    print(f"[{_ts()}] [precompute] failure_modes.json: {len(failure_modes)} modes", flush=True)

    # === chain_of_custody.json: commits + sessions ===
    chain = [
        {"sha": "a4f1c6f", "label": "pre-reg v1", "date": "2026-04-27",
         "summary": "Locked protocol: catalogs, Mc, features, nulls, gates."},
        {"sha": "05a4b0f", "label": "PRA-2", "date": "2026-04-28",
         "summary": "ISC catalog, overlap rule shrink, N≥8 inclusion floor."},
        {"sha": "ea2c0a3", "label": "Session 9 audit", "date": "2026-04-29",
         "summary": "Cascadia AUC=1.00 audit → failure mode #3 identified."},
        {"sha": "c725e16", "label": "Session 15 validation", "date": "2026-04-30",
         "summary": "Mask + magnitude breakdown + novelty re-check (Trugman 2019 found)."},
        {"sha": "6e21c02", "label": "Paper v2 draft", "date": "2026-04-30",
         "summary": "Round E2 confirmation paper — TLS as centerpiece."},
        {"sha": "2671a0d", "label": "Session 17 controls", "date": "2026-05-02",
         "summary": "exp18 + exp19: TLS underperforms trivial scalar; all 5 controls pass."},
        {"sha": "b1c4299", "label": "Paper v3 rewrite", "date": "2026-05-05",
         "summary": "Trivial scalar as headline; TLS demoted to ablation; M-scaling with CIs."},
    ]
    (WEB_DATA / "chain_of_custody.json").write_text(json.dumps(chain, indent=2))
    print(f"[{_ts()}] [precompute] chain_of_custody.json: {len(chain)} commits", flush=True)

    # === meta.json ===
    meta = {
        "title": "Seismic Precursors: Cross-Regional Pre-Registered Evaluation",
        "subtitle": "Five failure modes + cross-regional confirmation of Trugman 2019 foreshock pervasiveness",
        "regions": QUAL_REGIONS,
        "mc_per_region": MC_PER_REGION,
        "n_windows": int(len(df)),
        "n_precursor": int((df["window_kind"] == "precursor").sum()),
        "n_null": int((df["window_kind"] == NULL_KIND).sum()),
        "honest_score": "5.5–6.0",
        "headline_auc": 0.779,
        "headline_z_block": 8.37,
        "headline_z_iid": 11.0,
        "ci95": [0.701, 0.858],
        "repo": "https://github.com/MattLoftus/seismic-precursors",
        "generated_at": dt.datetime.utcnow().isoformat() + "Z",
    }
    (WEB_DATA / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"[{_ts()}] [precompute] meta.json", flush=True)

    print(f"\n[{_ts()}] [precompute] all done. Files in {WEB_DATA}:", flush=True)
    for p in sorted(WEB_DATA.glob("*.json")):
        size_kb = p.stat().st_size / 1024
        print(f"[{_ts()}]   {p.name:<28s}  {size_kb:>7.1f} KB", flush=True)


if __name__ == "__main__":
    main()
