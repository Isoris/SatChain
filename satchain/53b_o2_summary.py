#!/usr/bin/env python3
"""53b_o2_summary.py  (2026-11-01, HB04 / HB14 / CA26 / CA28): the second-order (Miga) test summarised over the catalogue RUNS.
Input: 53b_topologies.tsv of one or more 53b runs (--scratch S for the cat53/{G}_{src}/53b layout, or --topologies G:src=path ...)
and the catalogue 59_catalogue_{G}.tsv (merge on run_id; 53b's array_id column IS the run_id in run mode, so family columns fill).
Procedure (now in the repo, was an undocumented scratch step): per genome, 'array' unit rows (= whole run) with order2_z defined are
TESTABLE; the reported set at min_transitions m = testable runs with var_n_trans >= m, i.e. transitions of the VARIANT-ONLY path the
test is run on (not n_trans of the full path); p = one-sided normal tail of z, z = (gain - mean) / sd over the first-order simulations
(their count is the order2_sims column, --sims when the column is absent; a normal tail of a 100-simulation z is an approximation and
is stated as such in the output); q = Benjamini-Hochberg WITHIN the reported set (the pre-2026-11-01 numbers took BH over all
z-defined rows of a genome and subset afterwards); a hit = q <= q_cut and order2_excess_bits >= excess floor, every combination of
--min-trans x --q x --excess is written so 'barely changes' claims can be checked.
Outputs (--dest): 53b_o2_runs.tsv (every testable run: p, q_bh_min<m> per set, bp, family), 53b_o2_sweep.tsv (genome, min_transitions,
q_cut, excess_floor, n_runs, n_testable, n_reported, n_hits, hits_Mb, reported_Mb, n_sims, p_model) and a console summary. Units are RUNS."""
import argparse, math, os, sys
import numpy as np, pandas as pd


def bh(p):
    """Benjamini-Hochberg q-values of a p array (monotone step-up)."""
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    adj = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]; q[o] = np.minimum(adj, 1.0)
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--scratch', help='folder with cat53/{G}_{src}/53b/53b_topologies.tsv and cat59/59_catalogue_{G}.tsv')
    ap.add_argument('--topologies', nargs='*', default=[], help='G:src=path to a 53b_topologies.tsv (alternative to --scratch)')
    ap.add_argument('--catalogue', nargs='*', default=[], help='G=path to 59_catalogue_{G}.tsv (alternative to --scratch)')
    ap.add_argument('--genomes', default='Gar,Mac'); ap.add_argument('--sources', default='existing,new')
    ap.add_argument('--min-trans', default='20,50,100,200', help='reported sets: variant-only transitions >= m')
    ap.add_argument('--q', default='0.05,0.01'); ap.add_argument('--excess', default='0,0.05,0.1')
    ap.add_argument('--sims', type=int, default=100, help='simulation count when the table has no order2_sims column (pre-2026-11-01 53b)')
    ap.add_argument('--dest', required=True)
    a = ap.parse_args(); os.makedirs(a.dest, exist_ok=True)
    G_ = a.genomes.split(','); S_ = a.sources.split(','); M_ = [int(x) for x in a.min_trans.split(',')]
    Q_ = [float(x) for x in a.q.split(',')]; X_ = [float(x) for x in a.excess.split(',')]
    tops = {tuple(k.split(':')): v for k, v in (x.split('=', 1) for x in a.topologies)}
    cats = dict(x.split('=', 1) for x in a.catalogue)
    if a.scratch:
        for G in G_:
            cats.setdefault(G, f'{a.scratch}/cat59/59_catalogue_{G}.tsv')
            for src in S_: tops.setdefault((G, src), f'{a.scratch}/cat53/{G}_{src}/53b/53b_topologies.tsv')
    rows = []
    for G in G_:
        if G not in cats: sys.exit(f'ERROR no catalogue for {G}')
        cat = pd.read_csv(cats[G], sep='\t', low_memory=False)
        fcol = 'family_v46' if 'family_v46' in cat else 'family'; ncol = 'name_v46' if 'name_v46' in cat else fcol
        fam = cat.groupby('run_id').agg(family=(fcol, 'first'), name=(ncol, 'first'), unit_bp=('unit_bp', 'median')).reset_index()
        for src in S_:
            p = tops.get((G, src))
            if not p or not os.path.exists(p): sys.exit(f'ERROR 53b_topologies.tsv for {G}:{src} missing ({p})')
            d = pd.read_csv(p, sep='\t', low_memory=False)
            d = d[d.unit == 'array'].rename(columns={'array_id': 'run_id'}).copy(); d['G'] = G; d['src'] = src; d['bp'] = d.end - d.start0
            if 'order2_sims' not in d: d['order2_sims'] = a.sims
            rows.append(d.merge(fam, on='run_id', how='left'))      # merge on run_id (CA28): family columns fill in run mode
    rows = [r for r in rows if len(r)]
    if not rows: sys.exit('ERROR no array-unit rows in any 53b_topologies.tsv')
    A = pd.concat(rows, ignore_index=True)
    for c in ('order2_z', 'order2_excess_bits', 'memory_z', 'var_n_trans', 'n_trans'): A[c] = pd.to_numeric(A[c], errors='coerce')
    A['p_one_sided'] = 0.5 * (1 + np.vectorize(math.erf)(-A.order2_z / math.sqrt(2)))   # normal tail of the simulation z (approximation)
    for m in M_: A[f'q_bh_min{m}'] = np.nan
    sweep = []
    for G, g in A.groupby('G'):
        t = g[g.order2_z.notna()]
        nsims = sorted(set(t.order2_sims.dropna().astype(int))) if len(t) else []
        print(f'{G}: {len(g)} runs, {len(t)} testable (order-2 z defined), testable bp {t.bp.sum() / 1e6:.2f} of {g.bp.sum() / 1e6:.2f} Mb; '
              f'first-order simulations per run: {nsims}; p = one-sided normal tail of z (approximation), q = BH within each reported set')
        for m in M_:
            rep = t[t.var_n_trans >= m]
            if len(rep): A.loc[rep.index, f'q_bh_min{m}'] = bh(rep.p_one_sided.values)
            q = A.loc[rep.index, f'q_bh_min{m}']
            for qc in Q_:
                for ex in X_:
                    h = rep[(q <= qc) & (rep.order2_excess_bits >= ex)]
                    sweep.append(dict(genome=G, min_transitions_variant_only=m, q_cut=qc, excess_floor_bits=ex, n_runs=len(g), n_testable=len(t), n_reported=len(rep),
                                      n_hits=len(h), hits_Mb=round(h.bp.sum() / 1e6, 3), reported_Mb=round(rep.bp.sum() / 1e6, 3),
                                      n_sims=','.join(map(str, nsims)), p_model='normal_tail_of_simulation_z'))
                    print(f'   var_trans>={m:<4} q<={qc:<5} excess>={ex:<4}: {len(h):5d} of {len(rep):5d} runs, {h.bp.sum() / 1e6:6.2f} of {rep.bp.sum() / 1e6:6.2f} Mb')
        m0 = M_[0]; rep = t[t.var_n_trans >= m0]; q = A.loc[rep.index, f'q_bh_min{m0}']; h = rep[q <= Q_[0]]
        print(f'   overlap with first-order memory (z>=3) at var_trans>={m0}, q<={Q_[0]}: order-2 hits {len(h)}, of which with memory {(h.memory_z >= 3).sum()}; memory hits {(rep.memory_z >= 3).sum()}')
        if 'family' in h and len(h):
            f = h.groupby('family').agg(n_runs=('run_id', 'size'), bp=('bp', 'sum'), z=('order2_z', 'median'), ex=('order2_excess_bits', 'median')).sort_values('bp', ascending=False).head(8)
            print('   top families by order-2 bp:'); print(f.to_string(float_format=lambda v: f'{v:.3g}'))
            print('   top contexts:', h.order2_top_context.value_counts().head(6).to_dict())
    A.to_csv(os.path.join(a.dest, '53b_o2_runs.tsv'), sep='\t', index=False)
    pd.DataFrame(sweep).to_csv(os.path.join(a.dest, '53b_o2_sweep.tsv'), sep='\t', index=False)
    print('wrote', os.path.join(a.dest, '53b_o2_runs.tsv'), os.path.join(a.dest, '53b_o2_sweep.tsv'))


if __name__ == '__main__':
    main()
