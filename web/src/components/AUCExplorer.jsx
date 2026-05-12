import { useMemo, useState } from 'react';
import {
  ResponsiveContainer, ComposedChart, Bar, Line, XAxis, YAxis, ReferenceLine,
  Tooltip, Legend, CartesianGrid,
} from 'recharts';
import { REGION_COLORS } from '../store.js';

export default function AUCExplorer({ data }) {
  const { aucCurves } = data;
  const [feature, setFeature] = useState('blog'); // blog | n
  const [subwindow, setSubwindow] = useState(5);   // 0..5

  const chartData = useMemo(() => {
    return aucCurves.sub_window_labels.map((label, i) => {
      const row = { label, idx: i };
      for (const r of Object.keys(aucCurves.regions)) {
        const arr = feature === 'blog'
          ? aucCurves.regions[r].sliding_blog
          : aucCurves.regions[r].sliding_n;
        row[r] = arr[i];
      }
      row.macro = feature === 'blog'
        ? aucCurves.macro_sliding_blog[i]
        : aucCurves.macro_sliding_n[i];
      return row;
    });
  }, [aucCurves, feature]);

  const selectedRow = chartData[subwindow];

  return (
    <div className="space-y-6">
      <section className="rounded-lg border border-slate-800 bg-slate-900/40 p-5">
        <h2 className="mb-1 text-base font-semibold text-slate-100">
          Sliding-sub-window AUC explorer
        </h2>
        <p className="text-sm text-slate-400">
          Each window is split into six 5-day sub-windows. For each sub-window position, this
          chart shows the AUC of "Benioff energy (or event count) in <em>this</em> sub-window
          only" against precursor vs Null A labels. Slide the sub-window earlier and watch the
          signal collapse to chance.
        </p>
      </section>

      <section className="flex flex-wrap items-center gap-4 rounded-lg border border-slate-800 bg-slate-900/40 p-4">
        <div className="flex items-center gap-2">
          <label className="text-xs uppercase tracking-wide text-slate-500">Feature</label>
          <select
            value={feature}
            onChange={(e) => setFeature(e.target.value)}
            className="rounded border border-slate-700 bg-slate-900 px-2 py-1 text-sm"
          >
            <option value="blog">log₁₀ Benioff</option>
            <option value="n">event count (n above Mc)</option>
          </select>
        </div>
        <div className="flex flex-1 items-center gap-3 min-w-[300px]">
          <label className="text-xs uppercase tracking-wide text-slate-500 whitespace-nowrap">
            Sub-window
          </label>
          <input
            type="range"
            min="0"
            max="5"
            value={subwindow}
            onChange={(e) => setSubwindow(parseInt(e.target.value))}
            className="flex-1 accent-cyan-400"
          />
          <span className="font-mono text-xs text-cyan-300 w-24 text-right">
            days {subwindow * 5}–{(subwindow + 1) * 5}
          </span>
        </div>
      </section>

      <section className="rounded-lg border border-slate-800 bg-slate-900/40 p-4">
        <ResponsiveContainer width="100%" height={360}>
          <ComposedChart data={chartData}>
            <CartesianGrid stroke="#1e293b" />
            <XAxis dataKey="label" stroke="#64748b" tick={{ fontSize: 11 }} />
            <YAxis stroke="#64748b" domain={[0.3, 1.0]} tick={{ fontSize: 11 }} />
            <ReferenceLine y={0.5} stroke="#475569" strokeDasharray="3 3" label={{
              value: 'chance', fill: '#64748b', fontSize: 10, position: 'insideTopRight',
            }} />
            <Tooltip
              contentStyle={{ backgroundColor: '#0f172a', border: '1px solid #1e293b' }}
              labelStyle={{ color: '#e2e8f0' }}
            />
            <Legend />
            {Object.entries(REGION_COLORS).map(([region, color]) => (
              <Line
                key={region}
                type="monotone"
                dataKey={region}
                stroke={color}
                strokeWidth={2}
                dot={{ r: 4 }}
                activeDot={{ r: 6 }}
              />
            ))}
            <Line
              type="monotone"
              dataKey="macro"
              stroke="#06b6d4"
              strokeWidth={3}
              strokeDasharray="6 3"
              dot={{ r: 5 }}
            />
          </ComposedChart>
        </ResponsiveContainer>
      </section>

      <section className="rounded-lg border border-slate-800 bg-slate-900/40 p-4">
        <h3 className="mb-3 text-sm font-semibold text-slate-100">
          Selected sub-window: <span className="text-cyan-300">{selectedRow.label}</span>
        </h3>
        <div className="grid grid-cols-2 gap-3 md:grid-cols-5">
          <Cell label="Macro" value={selectedRow.macro} highlight />
          {Object.keys(REGION_COLORS).map((r) => (
            <Cell key={r} label={r} value={selectedRow[r]} color={REGION_COLORS[r]} />
          ))}
        </div>
        <p className="mt-3 text-xs leading-relaxed text-slate-500">
          The signal lives <em>exclusively</em> in days 25–30 (the immediate-foreshock period).
          Days 0–25 all sit near 0.5. This is the foreshock-period-localization that justifies
          framing the result as cross-regional confirmation of Trugman 2019 rather than a
          25-day-ahead precursor.
        </p>
      </section>
    </div>
  );
}

function Cell({ label, value, color, highlight }) {
  return (
    <div
      className={
        'rounded border p-3 text-center ' +
        (highlight ? 'border-cyan-600 bg-cyan-950/30' : 'border-slate-700 bg-slate-900/40')
      }
    >
      <div className="text-xs uppercase tracking-wide" style={{ color: color || '#94a3b8' }}>
        {label}
      </div>
      <div className="mt-1 font-mono text-sm">
        {Number.isFinite(value) ? value.toFixed(3) : '—'}
      </div>
    </div>
  );
}
