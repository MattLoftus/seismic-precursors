export default function Overview({ data }) {
  const { meta } = data;
  return (
    <div className="space-y-6">
      <section className="rounded-lg border border-slate-800 bg-slate-900/40 p-5">
        <h2 className="mb-3 text-base font-semibold text-slate-100">The 3-line summary</h2>
        <p className="text-sm leading-relaxed text-slate-300">
          We pre-registered a CSEP-style cross-regional precursor evaluation, ran the locked
          protocol with one documented amendment, and surfaced five execution-time failure
          modes. The simplest possible scalar — log₁₀ Benioff energy in the final 5 days of the
          30-day pre-event window — reaches macro AUC = {meta.headline_auc} across 4 LORO splits
          (California, Cascadia, Turkey, Italy), with z<sub>block</sub> = +{meta.headline_z_block}{' '}
          under a dependence-aware null. The signal is exclusively foreshock-period-localized; it
          is the cross-regional restatement of Trugman &amp; Ross 2019.
        </p>
      </section>

      <section className="grid grid-cols-2 gap-4 md:grid-cols-4">
        {[
          { label: 'Macro AUC', value: meta.headline_auc, sub: 'log₁₀ Benioff last 5d' },
          { label: 'z (block null)', value: '+' + meta.headline_z_block, sub: 'dependence-aware' },
          { label: 'CI95', value: meta.ci95.join(' – '), sub: 'cross-region bootstrap' },
          { label: 'Honest score', value: meta.honest_score, sub: 'methods note tier' },
        ].map((s) => (
          <div key={s.label} className="rounded-md border border-slate-800 bg-slate-900/40 p-4">
            <div className="text-xs uppercase tracking-wide text-slate-500">{s.label}</div>
            <div className="mt-1 font-mono text-lg text-cyan-300">{s.value}</div>
            <div className="text-xs text-slate-500">{s.sub}</div>
          </div>
        ))}
      </section>

      <section className="grid grid-cols-2 gap-4 md:grid-cols-4">
        {meta.regions.map((r) => (
          <div key={r} className="rounded-md border border-slate-800 bg-slate-900/40 p-4">
            <div className="text-xs uppercase tracking-wide text-slate-500">{r}</div>
            <div className="mt-1 font-mono text-sm text-slate-300">
              Mc = {meta.mc_per_region[r].toFixed(2)}
            </div>
          </div>
        ))}
      </section>

      <section className="rounded-lg border border-slate-800 bg-slate-900/40 p-5">
        <h2 className="mb-3 text-base font-semibold text-slate-100">The arc</h2>
        <ol className="space-y-2 text-sm text-slate-300">
          <li>
            <span className="font-mono text-cyan-300">1.</span> Pre-reg v1 + PRA-2 amendment
            commit, locking catalogs / Mc / features / nulls / gates before evaluation.
          </li>
          <li>
            <span className="font-mono text-cyan-300">2.</span> Per-window scalar (8 features × 2
            nulls) and joint LR + RF classifier — all clean nulls (0 / 16 tests pass 3σ).
          </li>
          <li>
            <span className="font-mono text-cyan-300">3.</span> TLS-style Benioff-trajectory
            template scan passes 3σ at macro AUC = 0.704 — looks like the headline.
          </li>
          <li>
            <span className="font-mono text-cyan-300">4.</span> Cold-read peer review (3
            subagents) converges on "TLS may just be a slope detector." Trivial-baseline
            ablation: log₁₀ Benioff in last 5 days hits {meta.headline_auc}, beating TLS by 0.075
            in 4/4 splits.
          </li>
          <li>
            <span className="font-mono text-cyan-300">5.</span> Four controls on the trivial
            scalar (precursor shift, placebo shift, foreshock mask, block-permutation null,
            paired ΔAUC) — all pass.
          </li>
          <li>
            <span className="font-mono text-cyan-300">6.</span> M-scaling: clean monotonic 0.70 →
            0.84 → 0.89 across three M-bins, all with bootstrap CIs.
          </li>
          <li>
            <span className="font-mono text-cyan-300">7.</span> Paper v3 rewritten with trivial
            scalar as headline, TLS demoted to ablation, peer-review fixes applied.
          </li>
        </ol>
      </section>

      <p className="text-xs text-slate-500">
        {meta.n_windows} windows total ({meta.n_precursor} precursor + {meta.n_null} Null A).
        Last regenerated: {meta.generated_at}.
      </p>
    </div>
  );
}
