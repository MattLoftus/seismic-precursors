import { useState } from 'react';
import {
  ResponsiveContainer, BarChart, Bar, XAxis, YAxis, ErrorBar, ReferenceLine,
  Tooltip, CartesianGrid,
} from 'recharts';
import { REGION_COLORS } from '../store.js';

export default function MScaling({ data }) {
  const { mScaling } = data;
  const [bins, setBins] = useState('bins_3');
  const series = mScaling[bins];

  const chartData = series.map((b) => ({
    label: b.label,
    macro: b.macro_auc,
    err: [b.macro_auc - b.ci_lo, b.ci_hi - b.macro_auc],
    n: b.n_pre,
    ...b.per_region,
  }));

  return (
    <div className="space-y-6">
      <section className="rounded-lg border border-slate-800 bg-slate-900/40 p-5">
        <h2 className="mb-1 text-base font-semibold text-slate-100">Magnitude scaling</h2>
        <p className="text-sm text-slate-400">
          AUC of the trivial scalar (log₁₀ Benioff in last 5d) stratified by target mainshock
          magnitude. Monotonic increase with target M is consistent with the Gutenberg-Richter
          scaling of triggered seismicity — bigger mainshocks have proportionally more and
          larger preceding events on average.
        </p>
      </section>

      <section className="flex items-center gap-2 rounded-lg border border-slate-800 bg-slate-900/40 p-4">
        <label className="text-xs uppercase tracking-wide text-slate-500">Bins</label>
        <button
          className={'rounded px-3 py-1 text-sm ' +
            (bins === 'bins_2' ? 'bg-cyan-500/20 text-cyan-200' : 'text-slate-400 hover:bg-slate-800')}
          onClick={() => setBins('bins_2')}
        >2-bin</button>
        <button
          className={'rounded px-3 py-1 text-sm ' +
            (bins === 'bins_3' ? 'bg-cyan-500/20 text-cyan-200' : 'text-slate-400 hover:bg-slate-800')}
          onClick={() => setBins('bins_3')}
        >3-bin (finer)</button>
      </section>

      <section className="rounded-lg border border-slate-800 bg-slate-900/40 p-4">
        <ResponsiveContainer width="100%" height={320}>
          <BarChart data={chartData}>
            <CartesianGrid stroke="#1e293b" />
            <XAxis dataKey="label" stroke="#64748b" tick={{ fontSize: 11 }} />
            <YAxis stroke="#64748b" domain={[0.3, 1.0]} tick={{ fontSize: 11 }} />
            <ReferenceLine y={0.5} stroke="#475569" strokeDasharray="3 3" />
            <Tooltip
              contentStyle={{ backgroundColor: '#0f172a', border: '1px solid #1e293b' }}
              labelStyle={{ color: '#e2e8f0' }}
            />
            <Bar dataKey="macro" fill="#06b6d4">
              <ErrorBar dataKey="err" width={4} strokeWidth={2} stroke="#a5f3fc"
                direction="y" />
            </Bar>
          </BarChart>
        </ResponsiveContainer>
      </section>

      <section className="rounded-lg border border-slate-800 bg-slate-900/40 p-4 overflow-x-auto">
        <h3 className="mb-3 text-sm font-semibold text-slate-100">Per-region per-bin AUC + CI95</h3>
        <table className="w-full font-mono text-xs">
          <thead className="text-slate-500 uppercase">
            <tr>
              <th className="text-left py-2">Bin</th>
              <th className="text-right pr-3">n</th>
              <th className="text-right pr-3">Macro AUC</th>
              {Object.keys(REGION_COLORS).map((r) => (
                <th key={r} className="text-right pr-3" style={{ color: REGION_COLORS[r] }}>{r}</th>
              ))}
            </tr>
          </thead>
          <tbody>
            {series.map((b) => (
              <tr key={b.label} className="border-t border-slate-800">
                <td className="py-2 text-slate-200">{b.label}</td>
                <td className="text-right pr-3 text-slate-400">{b.n_pre}</td>
                <td className="text-right pr-3 text-cyan-300">
                  {b.macro_auc.toFixed(3)} <span className="text-slate-500">[{b.ci_lo.toFixed(2)},{b.ci_hi.toFixed(2)}]</span>
                </td>
                {Object.keys(REGION_COLORS).map((r) => {
                  const cell = b.per_region[r];
                  return (
                    <td key={r} className="text-right pr-3" style={{ color: REGION_COLORS[r] }}>
                      {cell && cell.auc != null
                        ? `${cell.auc.toFixed(3)} (n=${cell.n_pre})`
                        : '—'}
                    </td>
                  );
                })}
              </tr>
            ))}
          </tbody>
        </table>
      </section>
    </div>
  );
}
