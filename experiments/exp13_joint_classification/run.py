"""exp13 — Joint multi-feature classification with leave-one-region-out CV.

Pre-reg v1 §3.2 + §7 mentioned `scikit-learn (logistic regression, random
forest)` for cross-regional ML evaluation. exp01–12 ran per-feature AUCs;
this experiment trains a classifier over all 6 computable features
simultaneously to test whether a learned combination separates precursor
from Null A windows where no single feature does.

Approach:
    - 6 features: 3 catalog (Benioff total, Benioff curvature, n above Mc)
      + 3 waveform (spectral slope, entropy, HHT IMF1 IF). The b and b-drift
      features are excluded due to failure mode #3 (uncomputable in 3 of 4
      qualifying regions).
    - Per-region StandardScaler normalization: each region's features are
      z-scored using the union of precursor + null A within that region.
      This isolates the "precursor vs background" residual pattern from
      per-region absolute-scale heterogeneity.
    - LORO (leave-one-region-out) cross-validation: train on 3 of the 4
      qualifying regions (California, Cascadia, Turkey, Italy), test on the
      held-out region, compute ROC-AUC on the test set.
    - Two classifiers in parallel: logistic regression (class_weight='balanced')
      and random forest (200 trees, balanced class weights).

Pre-reg gates: macro AUC across the 4 LORO splits passes 3σ if
CI95-lower > 0.5 AND |effect| > 0.05 AND permutation z > 3.

Outputs:
    loro_auc_table.csv         per-classifier per-test-region AUC
    loro_roc_plot.png          ROC curves for each (classifier, test region)
    feature_importance.png     LR coefs + RF importances
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
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
EXP_DIR = Path(__file__).resolve().parent
SOURCE_CSV = ROOT / "experiments" / "exp12_macro_full_panel" / "full_panel_feature_summary.csv"

# Six features computable in essentially every window (the b features are
# excluded; failure mode #3). Five if you collapse Benioff total + curvature.
FEATURES = [
    "benioff_total_log10",
    "benioff_curv",
    "n_above_mc",
    "wf_spectral_slope",
    "wf_waveform_entropy",
    "wf_hht_imf1_if_median_hz",
]
QUAL_REGIONS = ["California", "Cascadia", "Turkey", "Italy"]
NULL_KIND = "null_A"

N_BOOT = 1000
N_PERM = 1000
RANDOM_SEED = 42


def _ts() -> str:
    return dt.datetime.now().strftime("%H:%M:%S")


def per_region_zscore(df: pd.DataFrame) -> pd.DataFrame:
    """Z-score each feature WITHIN each region using union of precursor + null A."""
    out = df.copy()
    for region in df["region"].unique():
        mask = df["region"] == region
        for col in FEATURES:
            vals = df.loc[mask, col].to_numpy(dtype=float)
            mu = np.nanmean(vals)
            sigma = np.nanstd(vals)
            if not np.isfinite(sigma) or sigma == 0:
                sigma = 1.0
            out.loc[mask, col] = (vals - mu) / sigma
    return out


def main() -> int:
    print(f"[{_ts()}] [exp13] start — joint multi-feature classification, LORO CV", flush=True)

    df = pd.read_csv(SOURCE_CSV)
    print(f"[{_ts()}] [exp13] loaded {len(df)} rows", flush=True)

    # Restrict to (qualifying region) × (precursor or null_A)
    df = df[df["region"].isin(QUAL_REGIONS) &
            df["window_kind"].isin(["precursor", NULL_KIND])].copy()
    df["y"] = (df["window_kind"] == "precursor").astype(int)
    print(f"[{_ts()}] [exp13] {len(df)} rows after region+null filter; "
          f"precursor={(df['y']==1).sum()}, null={(df['y']==0).sum()}", flush=True)

    # Drop rows with any NaN in the feature panel (after waveform pipeline some
    # rows are missing; we want a clean cross-region trainable dataset)
    n_before = len(df)
    df = df.dropna(subset=FEATURES).reset_index(drop=True)
    print(f"[{_ts()}] [exp13] dropped {n_before - len(df)} rows with NaN in feature "
          f"panel; {len(df)} clean rows remaining", flush=True)
    print(f"[{_ts()}] [exp13] per-region clean counts:", flush=True)
    for region in QUAL_REGIONS:
        sub = df[df["region"] == region]
        n_pre = int((sub["y"] == 1).sum())
        n_null = int((sub["y"] == 0).sum())
        print(f"[{_ts()}]   {region:<11s}  precursor={n_pre:>3d}  null_A={n_null:>3d}",
              flush=True)

    # Per-region z-score
    df_z = per_region_zscore(df)

    # === LORO CV ===
    rng = np.random.default_rng(RANDOM_SEED)
    loro_results = []
    classifiers = {
        "logreg": lambda: LogisticRegression(
            class_weight="balanced", max_iter=2000, random_state=RANDOM_SEED),
        "rf": lambda: RandomForestClassifier(
            n_estimators=200, max_depth=None, class_weight="balanced",
            random_state=RANDOM_SEED, n_jobs=-1),
    }

    feature_importances = {clf: {f: [] for f in FEATURES} for clf in classifiers}

    for clf_name, clf_factory in classifiers.items():
        print(f"\n[{_ts()}] [exp13] === {clf_name} LORO CV ===", flush=True)
        for held_out in QUAL_REGIONS:
            train = df_z[df_z["region"] != held_out]
            test = df_z[df_z["region"] == held_out]
            X_train = train[FEATURES].to_numpy(dtype=float)
            y_train = train["y"].to_numpy()
            X_test = test[FEATURES].to_numpy(dtype=float)
            y_test = test["y"].to_numpy()
            if (y_test == 1).sum() < 2 or (y_test == 0).sum() < 2:
                print(f"[{_ts()}] [exp13]   {held_out}: insufficient class counts", flush=True)
                continue

            clf = clf_factory()
            clf.fit(X_train, y_train)
            scores = clf.predict_proba(X_test)[:, 1]
            auc = roc_auc_score(y_test, scores)

            # Bootstrap CI on this region's AUC
            boot = np.empty(N_BOOT)
            n_test = len(y_test)
            for i in range(N_BOOT):
                idx = rng.integers(0, n_test, size=n_test)
                if (y_test[idx] == 1).sum() < 2 or (y_test[idx] == 0).sum() < 2:
                    boot[i] = np.nan
                    continue
                boot[i] = roc_auc_score(y_test[idx], scores[idx])
            ci_lo, ci_hi = np.nanpercentile(boot, [2.5, 97.5])

            # Feature importance (LR coef or RF importance)
            if clf_name == "logreg":
                imp = clf.coef_[0]
            else:
                imp = clf.feature_importances_
            for fi, fcol in enumerate(FEATURES):
                feature_importances[clf_name][fcol].append(float(imp[fi]))

            loro_results.append({
                "classifier": clf_name,
                "test_region": held_out,
                "n_test_pre": int((y_test == 1).sum()),
                "n_test_null": int((y_test == 0).sum()),
                "auc": float(auc),
                "ci_lo": float(ci_lo), "ci_hi": float(ci_hi),
                "y_test": y_test.tolist(),
                "scores": scores.tolist(),
            })
            print(f"[{_ts()}] [exp13]   test={held_out:<11s}  "
                  f"AUC={auc:.3f}  CI=[{ci_lo:.3f}, {ci_hi:.3f}]  "
                  f"n_pre={int((y_test==1).sum())}, n_null={int((y_test==0).sum())}",
                  flush=True)

    # === Macro AUC across LORO splits per classifier ===
    print(f"\n[{_ts()}] [exp13] === macro AUC + significance ===", flush=True)
    macro_table = []
    for clf_name in classifiers:
        clf_aucs = np.array([r["auc"] for r in loro_results
                              if r["classifier"] == clf_name])
        macro = float(np.mean(clf_aucs))

        # Cross-region bootstrap (resample regions)
        boot_macros = np.empty(N_BOOT)
        for i in range(N_BOOT):
            idx = rng.integers(0, len(clf_aucs), size=len(clf_aucs))
            boot_macros[i] = np.mean(clf_aucs[idx])
        ci_lo, ci_hi = np.nanpercentile(boot_macros, [2.5, 97.5])

        # Permutation z: shuffle y within each region's test set, recompute AUC
        # for that region with the trained model's scores, then macro.
        perm_macros = np.empty(N_PERM)
        for p_i in range(N_PERM):
            perm_aucs = []
            for r in loro_results:
                if r["classifier"] != clf_name:
                    continue
                y = np.array(r["y_test"])
                s = np.array(r["scores"])
                rng.shuffle(y)
                if (y == 1).sum() < 2 or (y == 0).sum() < 2:
                    perm_aucs.append(np.nan)
                    continue
                perm_aucs.append(roc_auc_score(y, s))
            perm_macros[p_i] = float(np.nanmean(perm_aucs))
        perm_mean = float(np.nanmean(perm_macros))
        perm_std = float(np.nanstd(perm_macros))
        z = (macro - perm_mean) / perm_std if perm_std > 0 else float("nan")
        from scipy.stats import norm
        p_two = 2 * (1 - norm.cdf(abs(z)))

        # Pre-reg gates
        ci_above = ci_lo > 0.5
        effect = abs(macro - 0.5)
        sigma3 = ci_above and effect > 0.05 and abs(z) > 3
        # Bonferroni: 2 classifiers; alpha = 0.05/2 = 0.025
        bonferroni_alpha = 0.05 / 2
        bonf = p_two < bonferroni_alpha

        per_region = {r["test_region"]: r["auc"] for r in loro_results
                      if r["classifier"] == clf_name}
        macro_table.append({
            "classifier": clf_name,
            "per_region_aucs": per_region,
            "macro_auc": macro, "ci_lo": float(ci_lo), "ci_hi": float(ci_hi),
            "perm_z": float(z), "perm_p": float(p_two),
            "passes_3sigma": bool(sigma3),
            "passes_bonferroni": bool(bonf),
        })
        print(f"[{_ts()}] [headline] {clf_name:<6s}  macro AUC={macro:.3f}  "
              f"CI=[{ci_lo:.3f}, {ci_hi:.3f}]  z={z:+.2f}  p={p_two:.2e}  "
              f"3σ={'Y' if sigma3 else 'n'}  Bonf={'Y' if bonf else 'n'}", flush=True)
        print(f"[{_ts()}]            per-region: {per_region}", flush=True)

    # === Plots ===
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), dpi=120)
    for ax, clf_name in zip(axes, classifiers):
        for r in loro_results:
            if r["classifier"] != clf_name:
                continue
            fpr, tpr, _ = roc_curve(np.array(r["y_test"]), np.array(r["scores"]))
            ax.plot(fpr, tpr, lw=1.5, alpha=0.85,
                    label=f"{r['test_region']} (AUC={r['auc']:.3f})")
        ax.plot([0, 1], [0, 1], color="black", ls="--", lw=0.8, alpha=0.5)
        ax.set_xlabel("FPR"); ax.set_ylabel("TPR")
        ax.set_title(f"{clf_name} LORO ROC curves")
        ax.legend(fontsize=8, loc="lower right")
        ax.grid(alpha=0.3)
    fig.suptitle("exp13: leave-one-region-out ROC, joint 6-feature classifiers",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(EXP_DIR / "loro_roc_plot.png")
    plt.close(fig)
    print(f"[{_ts()}] [plot] loro_roc_plot.png", flush=True)

    # Feature importance plot (averaged across LORO splits)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), dpi=120, sharey=True)
    for ax, clf_name in zip(axes, classifiers):
        means = np.array([np.mean(feature_importances[clf_name][f])
                          for f in FEATURES])
        stds = np.array([np.std(feature_importances[clf_name][f])
                         for f in FEATURES])
        order = np.argsort(np.abs(means))[::-1]
        x = np.arange(len(FEATURES))
        ax.barh(x, means[order], xerr=stds[order], color="#1f5fa6", alpha=0.6,
                edgecolor="black", capsize=3)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_yticks(x)
        ax.set_yticklabels([FEATURES[i] for i in order], fontsize=8)
        ax.set_title(f"{clf_name}: feature importance "
                     f"({'coef' if clf_name=='logreg' else 'gini'})")
        ax.grid(alpha=0.3, axis="x")
    fig.suptitle("Joint-classifier feature importance, mean ± std across 4 LORO splits",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(EXP_DIR / "feature_importance.png")
    plt.close(fig)
    print(f"[{_ts()}] [plot] feature_importance.png", flush=True)

    pd.DataFrame(loro_results).drop(columns=["y_test", "scores"]).to_csv(
        EXP_DIR / "loro_auc_table.csv", index=False)

    summary = {
        "experiment": "exp13_joint_classification",
        "timestamp_utc": dt.datetime.utcnow().isoformat() + "Z",
        "features": FEATURES,
        "qualifying_regions": QUAL_REGIONS,
        "null_kind": NULL_KIND,
        "macro_table": macro_table,
        "loro_per_split": [
            {k: v for k, v in r.items() if k not in ("y_test", "scores")}
            for r in loro_results
        ],
        "feature_importance_by_classifier": {
            clf: {f: float(np.mean(vals)) for f, vals in d.items()}
            for clf, d in feature_importances.items()
        },
    }
    with open(EXP_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[{_ts()}] [persist] summary.json", flush=True)

    n3 = sum(1 for r in macro_table if r["passes_3sigma"])
    nb = sum(1 for r in macro_table if r["passes_bonferroni"])
    print(f"\n[{_ts()}] [exp13] HEADLINE: {n3}/2 classifiers pass 3σ, "
          f"{nb}/2 pass Bonferroni", flush=True)
    print(f"[{_ts()}] [exp13] done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
