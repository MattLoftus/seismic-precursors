import {
  ResponsiveContainer, BarChart, Bar, XAxis, YAxis, ReferenceLine, Tooltip, Legend,
  CartesianGrid, Cell,
} from 'recharts';
import { REGION_COLORS } from '../store.js';

export default function TLSvsScalar({ data }) {
  const { comparison } = data;
  const chartData = Object.keys(REGION_COLORS).map((r) => ({
    region: r,
    TLS: comparison.tls.per_region[r],
    Scalar: comparison.scalar.per_region[r],
    delta: comparison.delta[r],
  }));

  return (
    <div className="space-y-6">
      <section className="rounded-lg border border-slate-800 bg-slate-900/40 p-5">
        <h2 className="mb-1 text-base font-semibold text-slate-100">
          TLS template vs trivial scalar
        </h2>
        <p className="text-sm text-slate-400">
          The TLS scan applies Pearson correlation between each window's 6-point Benioff
          trajectory and a precursor-mean template. The trivial scalar is just the last sub-
          window's log₁₀ Benioff energy. Pearson correlation normalizes by within-window
          standard deviation — discarding exactly the absolute-magnitude information that turns
          out to be the dominant signal.
        </p>
      </section>

      <section className="grid grid-cols-1 gap-4 md:grid-cols-2">
        <StatPanel
          title={comparison.tls.label}
          color="#f87171"
          macro={comparison.tls.macro}
          ci={[comparison.tls.ci_lo, comparison.tls.ci_hi]}
          z_iid={comparison.tls.z_iid}
        />
        <StatPanel
          title={comparison.scalar.label}
          color="#34d399"
          macro={comparison.scalar.macro}
          ci={[comparison.scalar.ci_lo, comparison.scalar.ci_hi]}
          z_iid={comparison.scalar.z_iid}
          z_block={comparison.scalar.z_block}
        />
      </section>

      <section className="rounded-lg border border-slate-800 bg-slate-900/40 p-4">
        <h3 className="mb-2 text-sm font-semibold text-slate-100">
          Per-region AUC: TLS vs trivial scalar (Δ in 4 of 4 LORO splits)
        </h3>
        <ResponsiveContainer width="100%" height={320}>
          <BarChart data={chartData}>
            <CartesianGrid stroke="#1e293b" />
            <XAxis dataKey="region" stroke="#64748b" />
            <YAxis stroke="#64748b" domain={[0.5, 1.0]} tick={{ fontSize: 11 }} />
            <ReferenceLine y={0.5} stroke="#475569" strokeDasharray="3 3" />
            <Tooltip
              contentStyle={{ backgroundColor: '#0f172a', border: '1px solid #1e293b' }}
              labelStyle={{ color: '#e2e8f0' }}
            />
            <Legend />
            <Bar dataKey="TLS" fill="#f87171" name="TLS template" />
            <Bar dataKey="Scalar" fill="#34d399" name="Trivial scalar" />
          </BarChart>
        </ResponsiveContainer>
      </section>

      <section className="rounded-lg border border-slate-800 bg-slate-900/40 p-4">
        <h3 className="mb-3 text-sm font-semibold text-slate-100">Δ AUC (Scalar − TLS)</h3>
        <table className="w-full font-mono text-sm">
          <thead className="text-slate-500 text-xs uppercase">
            <tr>
              <th className="text-left py-2">Region</th>
              <th className="text-right">TLS</th>
              <th className="text-right">Scalar</th>
              <th className="text-right">Δ</th>
            </tr>
          </thead>
          <tbody>
            {chartData.map((row) => (
              <tr key={row.region} className="border-t border-slate-800">
                <td className="py-2 text-slate-300">{row.region}</td>
                <td className="text-right text-red-300">{row.TLS.toFixed(3)}</td>
                <td className="text-right text-emerald-300">{row.Scalar.toFixed(3)}</td>
                <td className="text-right text-cyan-300">+{row.delta.toFixed(3)}</td>
              </tr>
            ))}
            <tr className="border-t border-slate-700 font-semibold">
              <td className="py-2 text-slate-200">Macro</td>
              <td className="text-right text-red-300">{comparison.tls.macro.toFixed(3)}</td>
              <td className="text-right text-emerald-300">{comparison.scalar.macro.toFixed(3)}</td>
              <td className="text-right text-cyan-300">+{comparison.macro_delta.toFixed(3)}</td>
            </tr>
          </tbody>
        </table>
        <p className="mt-4 text-xs leading-relaxed text-slate-500">
          The trivial scalar beats TLS in 4 of 4 LORO splits. This is now lesson #122 in the
          cross-project research-learnings doc: "before publishing any similarity-based
          methodology, run the dumbest-possible-scalar baseline on the same data."
        </p>
      </section>
    </div>
  );
}

function StatPanel({ title, color, macro, ci, z_iid, z_block }) {
  return (
    <div className="rounded-lg border p-4" style={{ borderColor: color + '50' }}>
      <div className="text-xs uppercase tracking-wide" style={{ color }}>{title}</div>
      <div className="mt-2 font-mono text-2xl text-slate-100">{macro.toFixed(3)}</div>
      <div className="mt-1 text-xs text-slate-500">
        CI95: [{ci[0].toFixed(3)}, {ci[1].toFixed(3)}]
      </div>
      <div className="mt-1 text-xs text-slate-500">
        z_iid = +{z_iid.toFixed(2)}
        {z_block != null && (
          <>
            {' · '}
            z_block = +{z_block.toFixed(2)}
          </>
        )}
      </div>
    </div>
  );
}
