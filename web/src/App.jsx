import { useEffect } from 'react';
import { useStore } from './store.js';
import Overview from './components/Overview.jsx';
import AUCExplorer from './components/AUCExplorer.jsx';
import TrajectoryViewer from './components/TrajectoryViewer.jsx';
import TLSvsScalar from './components/TLSvsScalar.jsx';
import MScaling from './components/MScaling.jsx';
import FailureModes from './components/FailureModes.jsx';
import ChainOfCustody from './components/ChainOfCustody.jsx';

const TABS = [
  { id: 'overview', label: 'Overview' },
  { id: 'auc', label: 'AUC Explorer' },
  { id: 'traj', label: 'Trajectory Viewer' },
  { id: 'compare', label: 'TLS vs Scalar' },
  { id: 'mscale', label: 'M-Scaling' },
  { id: 'modes', label: 'Failure Modes' },
  { id: 'chain', label: 'Chain of Custody' },
];

export default function App() {
  const { data, loading, error, tab, setTab, loadAll } = useStore();

  useEffect(() => {
    loadAll();
  }, [loadAll]);

  if (loading) {
    return (
      <div className="flex h-screen items-center justify-center text-slate-400">
        Loading experiment data…
      </div>
    );
  }
  if (error) {
    return (
      <div className="flex h-screen items-center justify-center text-red-400">
        Load error: {error}
      </div>
    );
  }

  return (
    <div className="min-h-screen text-slate-200">
      <header className="border-b border-slate-800 bg-slate-900/50 px-6 py-4">
        <div className="mx-auto max-w-6xl">
          <h1 className="text-xl font-semibold text-slate-100">{data.meta.title}</h1>
          <p className="mt-1 text-sm text-slate-400">{data.meta.subtitle}</p>
        </div>
      </header>
      <nav className="border-b border-slate-800 bg-slate-900/30 px-6">
        <div className="mx-auto flex max-w-6xl gap-1 overflow-x-auto">
          {TABS.map((t) => (
            <button
              key={t.id}
              onClick={() => setTab(t.id)}
              className={
                'whitespace-nowrap border-b-2 px-3 py-3 text-sm font-medium transition-colors ' +
                (tab === t.id
                  ? 'border-cyan-400 text-cyan-300'
                  : 'border-transparent text-slate-400 hover:text-slate-200')
              }
            >
              {t.label}
            </button>
          ))}
        </div>
      </nav>
      <main className="mx-auto max-w-6xl px-6 py-6">
        {tab === 'overview' && <Overview data={data} />}
        {tab === 'auc' && <AUCExplorer data={data} />}
        {tab === 'traj' && <TrajectoryViewer data={data} />}
        {tab === 'compare' && <TLSvsScalar data={data} />}
        {tab === 'mscale' && <MScaling data={data} />}
        {tab === 'modes' && <FailureModes data={data} />}
        {tab === 'chain' && <ChainOfCustody data={data} />}
      </main>
      <footer className="border-t border-slate-800 bg-slate-900/30 px-6 py-4">
        <div className="mx-auto max-w-6xl text-xs text-slate-500">
          Honest score {data.meta.honest_score}. Headline macro AUC ={' '}
          <span className="text-slate-300">{data.meta.headline_auc}</span>, z_block ={' '}
          <span className="text-slate-300">+{data.meta.headline_z_block}</span>.{' '}
          <a className="text-cyan-400 hover:text-cyan-300" href={data.meta.repo}>
            GitHub repo
          </a>
        </div>
      </footer>
    </div>
  );
}
