import { useMemo, useState } from 'react';
import {
  ResponsiveContainer, LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend,
} from 'recharts';
import { REGION_COLORS } from '../store.js';

export default function TrajectoryViewer({ data }) {
  const { windows } = data;
  const [region, setRegion] = useState('California');
  const [feature, setFeature] = useState('b_traj'); // b_traj | n_traj
  const [windowId, setWindowId] = useState(null);

  const filtered = useMemo(
    () => windows.filter((w) => w.region === region),
    [windows, region]
  );
  const precursors = filtered.filter((w) => w.kind === 'precursor');
  const nulls = filtered.filter((w) => w.kind === 'null_A');

  const selected = windowId != null ? windows.find((w) => w.id === windowId) : null;

  // Precursor mean trajectory + 1 standard-deviation band
  const stats = useMemo(() => {
    if (!filtered.length) return null;
    const dim = 6;
    const meanPre = Array(dim).fill(0);
    const meanNull = Array(dim).fill(0);
    const sqPre = Array(dim).fill(0);
    const sqNull = Array(dim).fill(0);
    for (const w of precursors) {
      const t = w[feature];
      for (let i = 0; i < dim; i++) {
        meanPre[i] += t[i];
        sqPre[i] += t[i] * t[i];
      }
    }
    for (const w of nulls) {
      const t = w[feature];
      for (let i = 0; i < dim; i++) {
        meanNull[i] += t[i];
        sqNull[i] += t[i] * t[i];
      }
    }
    const np = precursors.length, nn = nulls.length;
    const data = [];
    for (let i = 0; i < dim; i++) {
      const mp = meanPre[i] / np;
      const mn = meanNull[i] / nn;
      const sdP = Math.sqrt(Math.max(sqPre[i] / np - mp * mp, 0));
      const sdN = Math.sqrt(Math.max(sqNull[i] / nn - mn * mn, 0));
      data.push({
        bin: `${i * 5}–${(i + 1) * 5}d`,
        precMean: mp,
        nullMean: mn,
        precHi: mp + sdP,
        precLo: mp - sdP,
        nullHi: mn + sdN,
        nullLo: mn - sdN,
        sel: selected ? selected[feature][i] : null,
      });
    }
    return data;
  }, [precursors, nulls, feature, selected]);

  return (
    <div className="space-y-6">
      <section className="rounded-lg border border-slate-800 bg-slate-900/40 p-5">
        <h2 className="mb-1 text-base font-semibold text-slate-100">Trajectory viewer</h2>
        <p className="text-sm text-slate-400">
          Each window's 6-point trajectory. The mean precursor and mean Null A curves diverge in
          the last sub-window — that's where the foreshock energy lives. Click any window from
          the list below to overlay its individual trajectory.
        </p>
      </section>

      <section className="flex flex-wrap gap-3 rounded-lg border border-slate-800 bg-slate-900/40 p-4">
        <Selector label="Region" value={region} onChange={setRegion}
          options={Object.keys(REGION_COLORS)} />
        <Selector label="Feature" value={feature} onChange={setFeature}
          options={[['b_traj', 'log₁₀ Benioff'], ['n_traj', 'n above Mc']]} />
      </section>

      <section className="rounded-lg border border-slate-800 bg-slate-900/40 p-4">
        <ResponsiveContainer width="100%" height={340}>
          <LineChart data={stats}>
            <CartesianGrid stroke="#1e293b" />
            <XAxis dataKey="bin" stroke="#64748b" tick={{ fontSize: 11 }} />
            <YAxis stroke="#64748b" tick={{ fontSize: 11 }} />
            <Tooltip
              contentStyle={{ backgroundColor: '#0f172a', border: '1px solid #1e293b' }}
              labelStyle={{ color: '#e2e8f0' }}
            />
            <Legend />
            <Line dataKey="precMean" stroke="#f87171" strokeWidth={3} name="precursor mean" />
            <Line dataKey="precHi" stroke="#7f1d1d" strokeWidth={1} strokeDasharray="3 3"
              name="precursor +1σ" dot={false} />
            <Line dataKey="precLo" stroke="#7f1d1d" strokeWidth={1} strokeDasharray="3 3"
              name="precursor −1σ" dot={false} />
            <Line dataKey="nullMean" stroke="#60a5fa" strokeWidth={3} name="null A mean" />
            <Line dataKey="nullHi" stroke="#1e3a8a" strokeWidth={1} strokeDasharray="3 3"
              name="null +1σ" dot={false} />
            <Line dataKey="nullLo" stroke="#1e3a8a" strokeWidth={1} strokeDasharray="3 3"
              name="null −1σ" dot={false} />
            {selected && (
              <Line dataKey="sel" stroke="#facc15" strokeWidth={3} dot={{ r: 5 }}
                name={`selected: ${selected.kind} #${selected.id}`} />
            )}
          </LineChart>
        </ResponsiveContainer>
      </section>

      <section className="grid grid-cols-1 gap-4 md:grid-cols-2">
        <WindowList title={`Precursor windows (${precursors.length})`} windows={precursors}
          selectedId={windowId} onSelect={setWindowId} />
        <WindowList title={`Null A windows (${nulls.length})`} windows={nulls}
          selectedId={windowId} onSelect={setWindowId} />
      </section>
    </div>
  );
}

function Selector({ label, value, onChange, options }) {
  const opts = options.map((o) => (Array.isArray(o) ? o : [o, o]));
  return (
    <div className="flex items-center gap-2">
      <label className="text-xs uppercase tracking-wide text-slate-500">{label}</label>
      <select
        value={value}
        onChange={(e) => onChange(e.target.value)}
        className="rounded border border-slate-700 bg-slate-900 px-2 py-1 text-sm"
      >
        {opts.map(([v, l]) => (
          <option key={v} value={v}>{l}</option>
        ))}
      </select>
    </div>
  );
}

function WindowList({ title, windows, selectedId, onSelect }) {
  return (
    <div className="rounded-lg border border-slate-800 bg-slate-900/40 p-4">
      <h3 className="mb-3 text-sm font-semibold text-slate-100">{title}</h3>
      <div className="h-64 overflow-y-auto space-y-1 font-mono text-xs">
        {windows.map((w) => (
          <button
            key={w.id}
            onClick={() => onSelect(w.id === selectedId ? null : w.id)}
            className={
              'block w-full text-left px-2 py-1 rounded transition-colors ' +
              (w.id === selectedId
                ? 'bg-yellow-500/20 text-yellow-200'
                : 'text-slate-400 hover:bg-slate-800 hover:text-slate-200')
            }
          >
            #{String(w.id).padStart(3, '0')} {w.t_start.slice(0, 10)}
            {w.target_M != null && (
              <span className="ml-2 text-slate-500">M{w.target_M.toFixed(1)}</span>
            )}
          </button>
        ))}
      </div>
    </div>
  );
}
