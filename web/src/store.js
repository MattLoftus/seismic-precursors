import { create } from 'zustand';

export const useStore = create((set) => ({
  data: null,
  loading: true,
  error: null,
  tab: 'overview',
  setTab: (tab) => set({ tab }),

  loadAll: async () => {
    try {
      const [meta, windows, aucCurves, comparison, mScaling, failureModes, chain] =
        await Promise.all([
          fetch('/data/meta.json').then((r) => r.json()),
          fetch('/data/windows.json').then((r) => r.json()),
          fetch('/data/auc_curves.json').then((r) => r.json()),
          fetch('/data/comparison.json').then((r) => r.json()),
          fetch('/data/m_scaling.json').then((r) => r.json()),
          fetch('/data/failure_modes.json').then((r) => r.json()),
          fetch('/data/chain_of_custody.json').then((r) => r.json()),
        ]);
      set({
        data: { meta, windows, aucCurves, comparison, mScaling, failureModes, chain },
        loading: false,
      });
    } catch (e) {
      console.error('load failed', e);
      set({ error: String(e), loading: false });
    }
  },
}));

export const REGION_COLORS = {
  California: '#60a5fa',
  Cascadia: '#34d399',
  Turkey: '#f59e0b',
  Italy: '#f472b6',
};
