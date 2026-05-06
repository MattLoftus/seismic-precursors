"""exp19 — Four controls on the exp18 trivial-scalar finding.

exp18 found that the simplest scalar — log10 Benioff energy in the last
5 days of a precursor window — reaches macro AUC=0.779 across 4 LORO
regions, beating the elaborate exp14 TLS template (AUC=0.704).

Before reframing the paper around this scalar, run the same controls the
reviewers demanded for TLS, plus the dependence-aware null they asked
for in addition:

  C1: PRECURSOR WINDOW-SHIFT — shift precursor windows back 30 days
      (so the "last 5 days" become days [-35, -30] before mainshock).
      Recompute A2_blog_last5d. Expectation: AUC should collapse to
      ~0.5 if the signal is foreshock-localized.

  C2: SYMMETRIC PLACEBO SHIFT — shift BOTH precursor and Null A
      windows back 30 days. Expectation: macro AUC should remain ~0.5
      (this is a no-op for any temporal-confound-free baseline).

  C3: FORESHOCK-PERIOD MASK — use log10 Benioff over days 0-25 only
      (drop the last 5 days entirely). Expectation: AUC collapses to
      ~0.5 if signal is foreshock-localized.

  C4: BLOCK-PERMUTATION NULL — instead of i.i.d. label permutation,
      shuffle labels in temporal blocks of size 30 days within each
      region (or use a circular time shift). Expectation: if the
      original z=11 was inflated by autocorrelation, block-permutation
      z should be substantially smaller. Reviewer 3's prediction:
      drops from ~11 to ~3-4.

  C5: PAIRED ΔAUC TEST (main vs C1-shifted) — bootstrap the difference
      between main and shifted AUCs on the same window resamples.
      Tests reviewer 3's concern that the two AUCs were never formally
      compared.

Outputs:
  controls_table.csv       — per-control macro AUC + CI + z + p
  per_split_table.csv      — per-control × per-LORO-region detail
  controls_plot.png        — bar chart of all controls vs main result
  paired_delta_dist.png    — bootstrap distribution of ΔAUC main−shifted
  summary.json
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
from scipy.stats import norm
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
EXP_DIR = Path(__file__).resolve().parent
EXP07_CSV = ROOT / "experiments" / "exp07_macro_pra2" / "feature_summary.csv"

MC_PER_REGION = {
    "California": 3.50, "Cascadia": 3.50, "Turkey": 3.10, "Italy": 3.50,
}
QUAL_REGIONS = ["California", "Cascadia", "Turkey", "Italy"]
NULL_KIND = "null_A"
N_SUBWINDOWS = 6
SUBWINDOW_DAYS = 5

N_BOOT = 1000
N_PERM = 1000
RANDOM_SEED = 42

MAIN_AUC = 0.779   # exp18 A2_blog_last5d macro
TLS_AUC = 0.704    # exp14 reference

BLOCK_SIZE_DAYS = 30  # for block permutation


def _ts() -> str:
    return dt.datetime.now().strftime("%H:%M:%S")


def load_catalog(region: str) -> pd.DataFrame:
    candidates = [
        ROOT / "experiments" / "exp06_cross_regional_macro" / f"catalog_{region}.csv",
        ROOT / "experiments" / "exp07_macro_pra2" / f"catalog_{region}.csv",
    ]
    for p in candidates:
        if p.is_file():
            df = pd.read_csv(p)
            df["time"] = pd.to_datetime(df["time"], utc=True, format="ISO8601")
            return df
    raise FileNotFoundError(f"no catalog cache found for {region}")


def catalog_trajectory(catalog, t_start, t_end, mc, n_subwindows=N_SUBWINDOWS):
    duration = (t_end - t_start).total_seconds()
    edges = [t_start + pd.Timedelta(seconds=duration * k / n_subwindows)
             for k in range(n_subwindows + 1)]
    n_arr = np.zeros(n_subwindows)
    b_arr = np.zeros(n_subwindows)
    for k in range(n_subwindows):
        mask = ((catalog["time"] >= edges[k]) & (catalog["time"] < edges[k + 1])
                & (catalog["magnitude"] >= mc - 1e-9))
        sub = catalog.loc[mask]
        n_arr[k] = float(len(sub))
        if len(sub):
            energies = 10 ** (1.5 * sub["magnitude"].to_numpy() + 4.8)
            b_arr[k] = float(np.sum(np.sqrt(energies)))
    b_arr_log = np.log10(np.where(b_arr > 0, b_arr, 1.0))
    return n_arr, b_arr_log


def per_feature_auc(scores, y):
    mask = np.isfinite(scores)
    if mask.sum() < 4 or (y[mask] == 1).sum() < 2 or (y[mask] == 0).sum() < 2:
        return float("nan")
    return float(roc_auc_score(y[mask], scores[mask]))


def macro_with_inference(per_split, rng, perm_method="iid"):
    """Compute macro across LORO splits with bootstrap CI95 + permutation z.

    perm_method: 'iid' (label shuffle) or 'block' (30-day-block circular shift).
    """
    aucs = np.array([r["auc"] for r in per_split])
    if not np.any(np.isfinite(aucs)):
        return {"macro": float("nan"), "ci_lo": float("nan"), "ci_hi": float("nan"),
                "perm_z": float("nan"), "perm_p": float("nan"), "perm_method": perm_method,
                "per_region_aucs": {r["held_out_region"]: r["auc"] for r in per_split}}

    macro = float(np.nanmean(aucs))

    boot = np.empty(N_BOOT)
    for bi in range(N_BOOT):
        idx = rng.integers(0, len(aucs), size=len(aucs))
        boot[bi] = np.nanmean(aucs[idx])
    ci_lo, ci_hi = np.nanpercentile(boot, [2.5, 97.5])

    perm = np.empty(N_PERM)
    for p_i in range(N_PERM):
        ps = []
        for r in per_split:
            y = np.array(r["y_test"]); s = np.array(r["scores"]);
            t = np.array(r.get("t_start_unix", np.arange(len(y), dtype=float)))
            m = np.isfinite(s)
            y2 = y[m].copy(); s2 = s[m]; t2 = t[m]
            if perm_method == "iid":
                rng.shuffle(y2)
            elif perm_method == "block":
                # Circular shift labels by a random offset, after sorting
                # by time (so neighboring windows stay neighboring).
                order = np.argsort(t2)
                inv = np.argsort(order)
                y_sorted = y2[order]
                shift = int(rng.integers(0, len(y_sorted)))
                y_shifted = np.roll(y_sorted, shift)
                y2 = y_shifted[inv]
            else:
                raise ValueError(perm_method)
            if (y2 == 1).sum() < 2 or (y2 == 0).sum() < 2:
                ps.append(np.nan)
                continue
            ps.append(roc_auc_score(y2, s2))
        perm[p_i] = float(np.nanmean(ps))
    perm_mean = float(np.nanmean(perm))
    perm_std = float(np.nanstd(perm))
    z = (macro - perm_mean) / perm_std if perm_std > 0 else float("nan")
    p_two = float(2 * (1 - norm.cdf(abs(z))))
    return {
        "macro": macro, "ci_lo": float(ci_lo), "ci_hi": float(ci_hi),
        "perm_z": float(z), "perm_p": p_two, "perm_method": perm_method,
        "per_region_aucs": {r["held_out_region"]: r["auc"] for r in per_split},
    }


def compute_baseline_per_split(df_eval, baseline_extractor):
    """For a given baseline (callable: row -> scalar), run the LORO scan
    and return per_split list."""
    per_split = []
    scalars = np.array([baseline_extractor(i) for i in range(len(df_eval))])
    for held_out in QUAL_REGIONS:
        test_mask = df_eval["region"] == held_out
        test_idx = np.where(test_mask)[0]
        if len(test_idx) < 4:
            continue
        scores = scalars[test_idx]
        y_test = (df_eval.iloc[test_idx]["window_kind"] == "precursor").astype(int).to_numpy()
        t_unix = df_eval.iloc[test_idx]["t_start_unix"].to_numpy()
        auc = per_feature_auc(scores, y_test)
        per_split.append({
            "held_out_region": held_out,
            "n_test_pre": int((y_test == 1).sum()),
            "n_test_null": int((y_test == 0).sum()),
            "auc": float(auc),
            "y_test": y_test.tolist(),
            "scores": scores.tolist(),
            "t_start_unix": t_unix.tolist(),
        })
    return per_split


def main() -> int:
    print(f"[{_ts()}] [exp19] start — controls on exp18 trivial-scalar finding", flush=True)
    rng = np.random.default_rng(RANDOM_SEED)

    df = pd.read_csv(EXP07_CSV)
    df["t_start"] = pd.to_datetime(df["t_start"], utc=True, format="ISO8601")
    df["t_end"] = pd.to_datetime(df["t_end"], utc=True, format="ISO8601")
    df = df[df["region"].isin(QUAL_REGIONS) &
            df["window_kind"].isin(["precursor", NULL_KIND])].copy().reset_index(drop=True)
    df["t_start_unix"] = df["t_start"].astype("int64") // 10**9
    print(f"[{_ts()}] [exp19] {len(df)} windows loaded", flush=True)

    catalogs = {region: load_catalog(region) for region in QUAL_REGIONS}

    # ----- MAIN baseline (re-derive exp18 A2): log10 Benioff in last 5 days -----
    print(f"\n[{_ts()}] [exp19] === MAIN: A2_blog_last5d (re-derive exp18) ===", flush=True)
    n_traj_main = np.zeros((len(df), N_SUBWINDOWS))
    b_traj_main = np.zeros((len(df), N_SUBWINDOWS))
    for i, row in df.iterrows():
        cat = catalogs[row["region"]]
        n, b = catalog_trajectory(cat, row["t_start"], row["t_end"],
                                   MC_PER_REGION[row["region"]])
        n_traj_main[i] = n
        b_traj_main[i] = b
    main_per_split = compute_baseline_per_split(df, lambda i: b_traj_main[i, 5])
    main_macro = macro_with_inference(main_per_split, rng, perm_method="iid")
    main_macro["control"] = "MAIN_A2_blog_last5d"
    print(f"[{_ts()}] [headline] {main_macro['control']:<32s}  macro={main_macro['macro']:.3f}  "
          f"CI=[{main_macro['ci_lo']:.3f}, {main_macro['ci_hi']:.3f}]  "
          f"z={main_macro['perm_z']:+.2f}  p={main_macro['perm_p']:.2e}", flush=True)

    # ----- C1: shift precursor windows back 30 days -----
    print(f"\n[{_ts()}] [exp19] === C1: precursor windows shifted back 30 days ===", flush=True)
    df_c1 = df.copy()
    pre_mask = df_c1["window_kind"] == "precursor"
    df_c1.loc[pre_mask, "t_start"] = df_c1.loc[pre_mask, "t_start"] - pd.Timedelta(days=30)
    df_c1.loc[pre_mask, "t_end"] = df_c1.loc[pre_mask, "t_end"] - pd.Timedelta(days=30)
    n_traj_c1 = np.zeros((len(df_c1), N_SUBWINDOWS))
    b_traj_c1 = np.zeros((len(df_c1), N_SUBWINDOWS))
    for i, row in df_c1.iterrows():
        cat = catalogs[row["region"]]
        n, b = catalog_trajectory(cat, row["t_start"], row["t_end"],
                                   MC_PER_REGION[row["region"]])
        n_traj_c1[i] = n
        b_traj_c1[i] = b
    c1_per_split = compute_baseline_per_split(df_c1, lambda i: b_traj_c1[i, 5])
    c1_macro = macro_with_inference(c1_per_split, rng, perm_method="iid")
    c1_macro["control"] = "C1_precursor_shift_-30d"
    print(f"[{_ts()}] [headline] {c1_macro['control']:<32s}  macro={c1_macro['macro']:.3f}  "
          f"CI=[{c1_macro['ci_lo']:.3f}, {c1_macro['ci_hi']:.3f}]  "
          f"z={c1_macro['perm_z']:+.2f}  p={c1_macro['perm_p']:.2e}", flush=True)

    # ----- C2: shift BOTH precursor and null windows back 30 days (placebo) -----
    print(f"\n[{_ts()}] [exp19] === C2: BOTH classes shifted back 30 days (placebo) ===",
          flush=True)
    df_c2 = df.copy()
    df_c2["t_start"] = df_c2["t_start"] - pd.Timedelta(days=30)
    df_c2["t_end"] = df_c2["t_end"] - pd.Timedelta(days=30)
    b_traj_c2 = np.zeros((len(df_c2), N_SUBWINDOWS))
    for i, row in df_c2.iterrows():
        cat = catalogs[row["region"]]
        _, b = catalog_trajectory(cat, row["t_start"], row["t_end"],
                                   MC_PER_REGION[row["region"]])
        b_traj_c2[i] = b
    c2_per_split = compute_baseline_per_split(df_c2, lambda i: b_traj_c2[i, 5])
    c2_macro = macro_with_inference(c2_per_split, rng, perm_method="iid")
    c2_macro["control"] = "C2_both_shift_-30d_placebo"
    print(f"[{_ts()}] [headline] {c2_macro['control']:<32s}  macro={c2_macro['macro']:.3f}  "
          f"CI=[{c2_macro['ci_lo']:.3f}, {c2_macro['ci_hi']:.3f}]  "
          f"z={c2_macro['perm_z']:+.2f}  p={c2_macro['perm_p']:.2e}", flush=True)

    # ----- C3: foreshock-period mask — use log10 Benioff in days 0-25 -----
    print(f"\n[{_ts()}] [exp19] === C3: mask last 5 days, use days 0-25 only ===", flush=True)
    # Sum Benioff energy over sub-windows 0..4 (days 0-25), then take log10
    benioff_energy_first5 = (10 ** b_traj_main[:, :5]).sum(axis=1)
    blog_first5 = np.log10(np.where(benioff_energy_first5 > 1.5, benioff_energy_first5, 1.0))
    c3_per_split = compute_baseline_per_split(df, lambda i: blog_first5[i])
    c3_macro = macro_with_inference(c3_per_split, rng, perm_method="iid")
    c3_macro["control"] = "C3_mask_last5d_use_days0-25"
    print(f"[{_ts()}] [headline] {c3_macro['control']:<32s}  macro={c3_macro['macro']:.3f}  "
          f"CI=[{c3_macro['ci_lo']:.3f}, {c3_macro['ci_hi']:.3f}]  "
          f"z={c3_macro['perm_z']:+.2f}  p={c3_macro['perm_p']:.2e}", flush=True)

    # ----- C4: block-permutation null on MAIN -----
    print(f"\n[{_ts()}] [exp19] === C4: block-permutation null on MAIN ===", flush=True)
    c4_macro = macro_with_inference(main_per_split, rng, perm_method="block")
    c4_macro["control"] = "C4_main_block_permutation_null"
    c4_macro["macro"] = main_macro["macro"]  # same point estimate; just a different null
    c4_macro["ci_lo"] = main_macro["ci_lo"]; c4_macro["ci_hi"] = main_macro["ci_hi"]
    c4_macro["per_region_aucs"] = main_macro["per_region_aucs"]
    print(f"[{_ts()}] [headline] {c4_macro['control']:<32s}  macro={c4_macro['macro']:.3f}  "
          f"(same point est.)  z_block={c4_macro['perm_z']:+.2f}  p={c4_macro['perm_p']:.2e}",
          flush=True)
    print(f"[{_ts()}]   compare: z_iid={main_macro['perm_z']:+.2f} → z_block={c4_macro['perm_z']:+.2f}",
          flush=True)

    # ----- C5: paired ΔAUC test (main vs C1) -----
    print(f"\n[{_ts()}] [exp19] === C5: paired ΔAUC bootstrap (main − C1) ===", flush=True)
    deltas = np.empty(N_BOOT)
    for bi in range(N_BOOT):
        per_region_d = []
        for region in QUAL_REGIONS:
            test_idx = np.where(df["region"] == region)[0]
            if len(test_idx) < 4:
                continue
            samp = rng.choice(test_idx, size=len(test_idx), replace=True)
            y_main = (df.iloc[samp]["window_kind"] == "precursor").astype(int).to_numpy()
            s_main = b_traj_main[samp, 5]
            s_c1 = b_traj_c1[samp, 5]
            if (y_main == 1).sum() < 2 or (y_main == 0).sum() < 2:
                continue
            auc_main = roc_auc_score(y_main, s_main)
            auc_c1 = roc_auc_score(y_main, s_c1)
            per_region_d.append(auc_main - auc_c1)
        deltas[bi] = float(np.mean(per_region_d)) if per_region_d else float("nan")
    delta_mean = float(np.nanmean(deltas))
    delta_lo, delta_hi = np.nanpercentile(deltas, [2.5, 97.5])
    delta_p = float(2 * min(np.mean(deltas <= 0), np.mean(deltas >= 0)))
    print(f"[{_ts()}] [headline] paired ΔAUC (main − C1) = {delta_mean:+.3f}  "
          f"CI=[{delta_lo:+.3f}, {delta_hi:+.3f}]  p={delta_p:.2e}", flush=True)

    # ----- Persist -----
    all_controls = [main_macro, c1_macro, c2_macro, c3_macro, c4_macro]
    pd.DataFrame([{k: v for k, v in c.items() if k != "per_region_aucs"}
                  for c in all_controls]).to_csv(EXP_DIR / "controls_table.csv", index=False)
    rows = []
    for label, per_split in [("MAIN", main_per_split), ("C1", c1_per_split),
                              ("C2", c2_per_split), ("C3", c3_per_split)]:
        for r in per_split:
            rows.append({"control": label, **{k: v for k, v in r.items()
                                                if k not in ("y_test", "scores", "t_start_unix")}})
    pd.DataFrame(rows).to_csv(EXP_DIR / "per_split_table.csv", index=False)

    # ----- Plot 1: control bar chart -----
    fig, ax = plt.subplots(figsize=(11, 5), dpi=120)
    names = [c["control"] for c in all_controls] + ["TLS_exp14"]
    macros = [c["macro"] for c in all_controls] + [TLS_AUC]
    err_lo = [c["macro"] - c["ci_lo"] for c in all_controls] + [0.064]
    err_hi = [c["ci_hi"] - c["macro"] for c in all_controls] + [0.066]
    colors = ["#3a7d3a", "#7d8aa6", "#7d8aa6", "#7d8aa6", "#3a7d3a", "#cc4444"]
    x = np.arange(len(names))
    ax.bar(x, macros, color=colors, alpha=0.65, yerr=[err_lo, err_hi],
           capsize=5, edgecolor="black")
    for fi, c in enumerate(all_controls):
        for r_idx, region in enumerate(QUAL_REGIONS):
            if region not in c["per_region_aucs"]:
                continue
            ax.scatter(fi + 0.18 * (r_idx - 1.5) / 4, c["per_region_aucs"][region],
                       s=22, alpha=0.85, color=plt.cm.tab10(r_idx),
                       edgecolor="white", linewidth=0.5,
                       label=region if fi == 0 else None)
    ax.axhline(0.5, color="black", ls="--", lw=0.8, label="chance")
    ax.axhline(MAIN_AUC, color="#3a7d3a", ls=":", lw=1.2, label=f"main exp18={MAIN_AUC}")
    ax.set_xticks(x); ax.set_xticklabels(names, fontsize=8, rotation=15, ha="right")
    ax.set_ylim(0.30, 0.92)
    ax.set_ylabel("Macro AUC")
    ax.set_title("exp19: controls on the trivial-scalar finding (A2_blog_last5d)")
    ax.legend(fontsize=8, loc="lower right", ncol=2)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(EXP_DIR / "controls_plot.png")
    plt.close(fig)
    print(f"[{_ts()}] [plot] controls_plot.png", flush=True)

    # ----- Plot 2: paired ΔAUC distribution -----
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=120)
    ax.hist(deltas[np.isfinite(deltas)], bins=40, color="#3a7d3a", alpha=0.7, edgecolor="white")
    ax.axvline(0, color="black", ls="--", lw=0.8, label="ΔAUC = 0")
    ax.axvline(delta_mean, color="#cc4444", ls="-", lw=1.3, label=f"mean = {delta_mean:+.3f}")
    ax.axvline(delta_lo, color="#cc4444", ls=":", lw=1.0,
               label=f"CI95 [{delta_lo:+.3f}, {delta_hi:+.3f}]")
    ax.axvline(delta_hi, color="#cc4444", ls=":", lw=1.0)
    ax.set_xlabel("ΔAUC = AUC(main) − AUC(C1 shifted)")
    ax.set_ylabel("Bootstrap density")
    ax.set_title(f"exp19 C5: paired ΔAUC bootstrap (B={N_BOOT})")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(EXP_DIR / "paired_delta_dist.png")
    plt.close(fig)
    print(f"[{_ts()}] [plot] paired_delta_dist.png", flush=True)

    # ----- Verdict -----
    pass_c1 = c1_macro["macro"] < 0.60   # AUC collapses to ~chance
    pass_c2 = abs(c2_macro["macro"] - c1_macro["macro"]) < 0.05  # placebo ≈ shifted
    pass_c3 = c3_macro["macro"] < 0.60   # signal not in days 0-25
    pass_c4 = c4_macro["perm_z"] > 3      # signal robust to dependence-aware null
    pass_c5 = delta_lo > 0                 # main > shifted with CI excluding 0

    n_pass = sum([pass_c1, pass_c2, pass_c3, pass_c4, pass_c5])
    if n_pass == 5:
        verdict = "ALL_5_CONTROLS_PASS — trivial scalar finding is robust; reframe paper (Path A)"
    elif n_pass == 4:
        verdict = (f"4_OF_5_PASS — one weakness (failed: "
                   f"{['C1','C2','C3','C4','C5'][[pass_c1,pass_c2,pass_c3,pass_c4,pass_c5].index(False)]})"
                   "; revisit before reframe")
    elif n_pass >= 2:
        verdict = f"{n_pass}_OF_5_PASS — significant weaknesses; lean toward shelve (Path B)"
    else:
        verdict = f"{n_pass}_OF_5_PASS — trivial scalar fails; SHELVE (Path B)"

    summary = {
        "experiment": "exp19_trivial_scalar_controls",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "main_baseline_recomputed_macro": main_macro["macro"],
        "main_baseline_exp18_macro": MAIN_AUC,
        "controls": {c["control"]: {k: v for k, v in c.items()} for c in all_controls},
        "paired_delta_main_minus_c1": {
            "mean": delta_mean, "ci95_lo": float(delta_lo), "ci95_hi": float(delta_hi),
            "p_two_sided": delta_p,
        },
        "control_pass_flags": {
            "C1_precursor_shift_collapses": bool(pass_c1),
            "C2_placebo_no_effect": bool(pass_c2),
            "C3_mask_signal_in_last5d": bool(pass_c3),
            "C4_block_perm_z_above_3": bool(pass_c4),
            "C5_paired_delta_excludes_0": bool(pass_c5),
        },
        "n_controls_passed": int(n_pass),
        "verdict": verdict,
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[{_ts()}] [persist] summary.json", flush=True)

    print(f"\n[{_ts()}] [exp19] HEADLINE", flush=True)
    print(f"[{_ts()}]   main         AUC = {main_macro['macro']:.3f}  z_iid={main_macro['perm_z']:+.2f}",
          flush=True)
    print(f"[{_ts()}]   C1 shift     AUC = {c1_macro['macro']:.3f}  (collapse? {'Y' if pass_c1 else 'n'})",
          flush=True)
    print(f"[{_ts()}]   C2 placebo   AUC = {c2_macro['macro']:.3f}  (≈ C1? {'Y' if pass_c2 else 'n'})",
          flush=True)
    print(f"[{_ts()}]   C3 masked    AUC = {c3_macro['macro']:.3f}  (collapse? {'Y' if pass_c3 else 'n'})",
          flush=True)
    print(f"[{_ts()}]   C4 z_block   = {c4_macro['perm_z']:+.2f}  (>3? {'Y' if pass_c4 else 'n'})",
          flush=True)
    print(f"[{_ts()}]   C5 ΔAUC      = {delta_mean:+.3f}  CI=[{delta_lo:+.3f},{delta_hi:+.3f}]  "
          f"(>0? {'Y' if pass_c5 else 'n'})", flush=True)
    print(f"[{_ts()}]   VERDICT: {verdict}", flush=True)
    print(f"[{_ts()}] [exp19] done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
