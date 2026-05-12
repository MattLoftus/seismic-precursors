export default function ChainOfCustody({ data }) {
  const { chain, meta } = data;
  return (
    <div className="space-y-6">
      <section className="rounded-lg border border-slate-800 bg-slate-900/40 p-5">
        <h2 className="mb-1 text-base font-semibold text-slate-100">
          Chain of custody
        </h2>
        <p className="text-sm text-slate-400">
          Each commit is a single work session. Pre-registration and amendment both committed
          before evaluation; every subsequent session cites both SHAs. This is the
          publicly-verifiable artifact the paper points to.
        </p>
      </section>

      <ol className="relative space-y-3 border-l border-slate-700 pl-6">
        {chain.map((c, i) => (
          <li key={c.sha} className="relative">
            <span
              className="absolute -left-[34px] mt-1.5 inline-block h-3 w-3 rounded-full"
              style={{
                backgroundColor: i === chain.length - 1 ? '#06b6d4' : '#64748b',
                boxShadow: i === chain.length - 1 ? '0 0 0 4px #155e75' : 'none',
              }}
            />
            <article className="rounded-md border border-slate-800 bg-slate-900/40 p-4">
              <header className="flex flex-wrap items-baseline gap-3">
                <a
                  href={`${meta.repo}/commit/${c.sha}`}
                  className="font-mono text-sm text-cyan-300 hover:text-cyan-200"
                >
                  {c.sha}
                </a>
                <span className="text-sm font-semibold text-slate-100">{c.label}</span>
                <span className="ml-auto font-mono text-xs text-slate-500">{c.date}</span>
              </header>
              <p className="mt-1 text-sm text-slate-300">{c.summary}</p>
            </article>
          </li>
        ))}
      </ol>
    </div>
  );
}
