# seismic-precursors web

Interactive exploration tool for the seismic-precursors project. Personal-
use companion to the paper / repo — not a deployable production site.

## What's here

Seven tabs:

1. **Overview** — 3-line summary, headline stats, the arc of the project.
2. **AUC Explorer** — slide the 5-day sub-window across the 30-day pre-event
   window and watch the AUC of "log₁₀ Benioff (or event count) in this
   sub-window" against precursor / Null A. Shows the foreshock-period
   localization viscerally.
3. **Trajectory Viewer** — per-region precursor-mean and Null-A-mean
   6-point trajectories with ±1σ bands; click any window to overlay.
4. **TLS vs Scalar** — side-by-side per-region AUC; the elaborate template-
   correlation apparatus (TLS, AUC=0.704) vs the trivial single scalar
   (log₁₀ Benioff in last 5d, AUC=0.779) that beats it in 4/4 LORO splits.
5. **M-Scaling** — 2-bin and 3-bin breakdown of trivial-scalar AUC by
   target mainshock magnitude, with bootstrap CI95s per bin and per region.
6. **Failure Modes** — five execution-time failure modes catalogued with
   diagnostic evidence and mitigation.
7. **Chain of Custody** — pre-reg v1 → PRA-2 → session commits → paper v3
   timeline with GitHub commit links.

## Run locally

```sh
cd web
npm install
npm run dev
# open http://localhost:5180/
```

## Regenerate data

The app loads static JSON from `public/data/`. If experiments are updated:

```sh
/usr/bin/python3 scripts/precompute_data.py
```

This reads from `experiments/exp{07,14,18,20}/` artifacts and the per-
region catalog caches (gitignored). The output is 6 JSON files totaling
~247 KB.

## Stack

Vite 7 + React 19 + Zustand + Tailwind + Recharts. No backend.
