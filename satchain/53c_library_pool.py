#!/usr/bin/env python3
"""
53c_library_pool.py  (v4.5.2; one topology library for both subgenomes, on the runs of record)

Pools the satellite segments (not 'unresolved' / 'too_short', >= --min-trans transitions) of several 53b runs, coarsens their transition
counts and compositions to the library alphabet (53b.coarse), computes D = D_trans + D_comp between all pairs, clusters them
(average linkage) at cuts taken from the pooled SPLIT-HALF NULL of 53b (53b_library_null.tsv: D of one segment against its own other
half; cuts = --null-quantiles of it, 2026-11-01 HB09 - no longer quantiles of the data's nearest-neighbour distance), and describes
each type: members per subgenome, centroid transition probabilities among V1 / V2 / V3+ / U and the annotation states, median variant
structure index, median annotation fraction, member segment profiles. A type with members in >= 2 tags is 'shared', else '<tag>_only'.
2026-11-01 (CA29): the composition is the stored 53b 'composition' column (state counts), not the out-transition row sums; tags,
--min-trans and the cut quantiles are arguments. Column array_id = run_id in run mode (the unit is a run segment).
Outputs (<dest>/): 53c_library_assignments.tsv, 53c_library_types.tsv (one row per type and cut), 53c_library_null.tsv
usage: 53c_library_pool.py <dest> <tag>=<53b dir> [<tag>=<53b dir> ...] [--min-trans 40] [--null-quantiles 0.5,0.9,0.95,0.99]
"""
import argparse, collections, csv, importlib.util, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('t53b', os.path.join(HERE, '53b_topology_compare.py')); T = importlib.util.module_from_spec(spec); spec.loader.exec_module(T)
ap = argparse.ArgumentParser(); ap.add_argument('dest'); ap.add_argument('runs', nargs='+', help='tag=<53b dir>')
ap.add_argument('--min-trans', type=int, default=40); ap.add_argument('--null-quantiles', type=lambda s: [float(x) for x in s.split(',')], default=[0.5, 0.9, 0.95, 0.99])
a = ap.parse_args(); DEST = a.dest; RUNS = [x.split('=', 1) for x in a.runs]; MIN_TRANS = a.min_trans; QS = a.null_quantiles; TAGS = [t for t, _ in RUNS]
os.makedirs(DEST, exist_ok=True)
pool = []; null = []
for tag, d in RUNS:
    top = {(r['array_id'], r['unit']): r for r in csv.DictReader(open(os.path.join(d, '53b_topologies.tsv')), delimiter='\t')}
    tr = collections.defaultdict(collections.Counter)
    for r in csv.DictReader(open(os.path.join(d, '53b_transitions.tsv')), delimiter='\t'): tr[(r['array_id'], r['unit'])][(r['from'], r['to'])] += int(r['count'])
    pn = os.path.join(d, '53b_library_null.tsv')
    if not os.path.exists(pn): sys.exit(f'ERROR {pn} missing: 53b must be the 2026-11-01 version (split-half null)')
    null += [float(r['split_half_D']) for r in csv.DictReader(open(pn), delimiter='\t') if r['split_half_D'] not in ('NA', '')]
    for k, r in top.items():
        if r['unit'] == 'array' or r['segment_profile'] in ('unresolved', 'too_short') or int(r['n_trans'] or 0) < MIN_TRANS: continue
        comp = collections.Counter({x.rsplit(':', 1)[0]: int(x.rsplit(':', 1)[1]) for x in r.get('composition', '').split(';') if x})
        if not comp:
            for (x, y), v in tr[k].items(): comp[x] += v    # fallback for a table without the composition column (pre-2026-11-01)
        t2, c2 = T.coarsen(tr[k], comp)
        pool.append(dict(tag=tag, array=k[0], unit=k[1], chrom=r['chrom'], start=r['start0'], end=r['end'], profile=r['segment_profile'],
                         si=r['var_structure_index'], si_z=r['var_structure_z'], annot=r['annotation_frac'], tr=t2, comp=c2))
n = len(pool); D = [[0.0] * n for _ in range(n)]
for i in range(n):
    for j in range(i + 1, n):
        dt = T.D_trans(pool[i]['tr'], pool[i]['comp'], pool[j]['tr'], pool[j]['comp']); dt = 1.0 if dt == 'NA' else dt
        D[i][j] = D[j][i] = dt + T.jsd(pool[i]['comp'], pool[j]['comp'])
nn = sorted(min(D[i][j] for j in range(n) if j != i) for i in range(n)) if n > 1 else [0.0]
null.sort(); cuts = [round(null[min(len(null) - 1, int(q * len(null)))], 4) if null else 0.0 for q in QS]
with open(os.path.join(DEST, '53c_library_null.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n'); w.writerow(['quantile', 'cut_D', 'n_null_segments', 'nn_D_median'])
    for q, c in zip(QS, cuts): w.writerow([q, c, len(null), round(nn[len(nn) // 2], 4)])
labs = []
for cut in cuts:
    cl = [[i] for i in range(n)]
    while len(cl) > 1:
        best = None
        for x in range(len(cl)):
            for y in range(x + 1, len(cl)):
                dd = sum(D[i][j] for i in cl[x] for j in cl[y]) / (len(cl[x]) * len(cl[y]))
                if best is None or dd < best[0]: best = (dd, x, y)
        if best[0] > cut: break
        cl[best[1]] += cl[best[2]]; del cl[best[2]]
    lab = {}
    for k, c in enumerate(sorted(cl, key=len, reverse=True)): lab.update({i: f'L{k + 1:02d}' for i in c})
    labs.append(lab)
with open(os.path.join(DEST, '53c_library_assignments.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['subgenome', 'array_id', 'unit', 'chrom', 'start0', 'end', 'segment_profile', 'var_structure_index', 'annotation_frac'] + [f'type_nullq{q}' for q in QS])
    for i, p in enumerate(pool): w.writerow([p['tag'], p['array'], p['unit'], p['chrom'], p['start'], p['end'], p['profile'], p['si'], p['annot']] + [l[i] for l in labs])
KEYS = [('V1', 'V1'), ('V1', 'V2'), ('V2', 'V1'), ('V2', 'V2'), ('V1', 'U'), ('U', 'U')]
num = lambda v: float(v) if v not in ('NA', '') else None
med = lambda v: sorted(v)[len(v) // 2] if v else 'NA'
with open(os.path.join(DEST, '53c_library_types.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['null_quantile', 'cut_D', 'type', 'n_segments'] + [f'n_{t}' for t in TAGS] + ['sharing', 'bp_total'] + [f'P({b}|{a_})' for a_, b in KEYS] +
               ['V1_frac', 'V2_frac', 'annotation_states_frac', 'median_var_SI', 'profiles', 'members'])
    for q, cut, lab in zip(QS, cuts, labs):
        for ty in sorted(set(lab.values())):
            m = [pool[i] for i in lab if lab[i] == ty]
            tr = collections.Counter(); comp = collections.Counter()
            for p in m: tr.update(p['tr']); comp.update(p['comp'])
            row = collections.Counter()
            for (x, y), v in tr.items(): row[x] += v
            tot = sum(comp.values()) or 1
            g = collections.Counter(p['tag'] for p in m)
            sh = 'shared' if len(g) > 1 else f'{next(iter(g))}_only'
            si = [num(p['si']) for p in m if num(p['si']) is not None]
            w.writerow([q, cut, ty, len(m)] + [g.get(t, 0) for t in TAGS] + [sh, sum(int(p['end']) - int(p['start']) for p in m)] +
                       [round(tr[k] / row[k[0]], 3) if row[k[0]] else 'NA' for k in KEYS] +
                       [round(comp['V1'] / tot, 3), round(comp['V2'] / tot, 3), round(1 - sum(comp[s] for s in ('V1', 'V2', 'V3+', 'U', 'OF')) / tot, 3),
                        round(med(si), 3) if si else 'NA', ';'.join(f'{k}:{v}' for k, v in collections.Counter(p['profile'] for p in m).most_common()),
                        ','.join(f"{p['tag']}:{p['array']}:{p['unit']}" for p in m)])
print(n, 'segments from', TAGS, '; nearest-neighbour D median', round(nn[len(nn) // 2], 3), '; split-half null n', len(null), 'median', round(null[len(null) // 2], 3) if null else 'NA', '; cuts', cuts, '; types', [len(set(l.values())) for l in labs])
for q, lab in zip(QS, labs):
    sh = collections.Counter()
    for ty in set(lab.values()):
        g = {pool[i]['tag'] for i in lab if lab[i] == ty}; sh['shared' if len(g) > 1 else 'single_subgenome'] += 1
    print(f'  null q{q}: {dict(sh)}')
