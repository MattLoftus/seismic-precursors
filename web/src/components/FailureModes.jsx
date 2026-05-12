export default function FailureModes({ data }) {
  return (
    <div className="space-y-4">
      <section className="rounded-lg border border-slate-800 bg-slate-900/40 p-5">
        <h2 className="mb-1 text-base font-semibold text-slate-100">
          Five execution-time failure modes
        </h2>
        <p className="text-sm text-slate-400">
          The methodological catalogue the paper contributes. Three (#1, #2, #4) are operational
          gotchas familiar to subduction-zone researchers; their combination under our protocol
          made them load-bearing. Two (#3, #5) are the more substantive contributions.
        </p>
      </section>

      {data.failureModes.map((m) => (
        <article
          key={m.n}
          className="rounded-lg border border-slate-800 bg-slate-900/40 p-5"
        >
          <header className="flex items-baseline gap-3">
            <span className="font-mono text-2xl text-cyan-400">#{m.n}</span>
            <h3 className="text-base font-semibold text-slate-100">{m.title}</h3>
          </header>
          <p className="mt-2 text-sm leading-relaxed text-slate-300">{m.summary}</p>

          <div className="mt-4 rounded-md border border-slate-800 bg-slate-950/40 p-3">
            <div className="mb-2 text-xs uppercase tracking-wide text-slate-500">Evidence</div>
            <dl className="grid grid-cols-1 gap-x-6 gap-y-1 text-xs font-mono md:grid-cols-2">
              {Object.entries(m.evidence).map(([k, v]) => (
                <div key={k} className="flex justify-between gap-3">
                  <dt className="text-slate-400">{k}</dt>
                  <dd className="text-slate-200 text-right">{v}</dd>
                </div>
              ))}
            </dl>
          </div>

          <div className="mt-3 text-xs">
            <span className="text-slate-500 uppercase tracking-wide">Mitigation: </span>
            <span className="text-slate-300">{m.mitigation}</span>
          </div>
        </article>
      ))}
    </div>
  );
}
