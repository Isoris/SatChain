#!/usr/bin/env python3
"""
53c_library_pool.py  (v4.5.2; one topology library for both subgenomes)

Pools the satellite segments (not 'unresolved', >= --min-trans transitions) of several 53b runs, coarsens their transition counts
to the library alphabet (53b.coarse), computes D = D_trans + D_comp between all pairs, clusters them (average linkage) at cuts
set by quantiles of the nearest-neighbour distance (a sweep), and describes each type:
  members per subgenome, centroid transition probabilities among V1 / V2 / V3+ / U and the annotation states, median variant
  structure index, median annotation fraction, member segment profiles.
A type with members in both subgenomes is 'shared'; otherwise 'Gar_only' / 'Mac_only'.
Outputs (<dest>/): 53c_library_assignments.tsv, 53c_library_types.tsv (one row per type and cut)
usage: 53c_library_pool.py <dest> <tag>=<53b dir> [<tag>=<53b dir> ...]
"""
import collections, csv, importlib.util, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('t53b', os.path.join(HERE, '53b_topology_compare.py')); T = importlib.util.module_from_spec(spec); spec.loader.exec_module(T)
DEST = sys.argv[1]; RUNS = [x.split('=', 1) for x in sys.argv[2:]]; MIN_TRANS = 40; QS = [0.25, 0.5, 0.75, 0.9]
os.makedirs(DEST, exist_ok=True)
pool = []
for tag, d in RUNS:
    top = {(r['array_id'], r['unit']): r for r in csv.DictReader(open(os.path.join(d, '53b_topologies.tsv')), delimiter='\t')}
    tr = collections.defaultdict(collections.Counter)
    for r in csv.DictReader(open(os.path.join(d, '53b_transitions.tsv')), delimiter='\t'): tr[(r['array_id'], r['unit'])][(r['from'], r['to'])] += int(r['count'])
    for k, r in top.items():
        if r['unit'] == 'array' or r['segment_profile'] == 'unresolved' or int(r['n_trans'] or 0) < MIN_TRANS: continue
        comp = collections.Counter()
        for (x, y), v in tr[k].items(): comp[x] += v
        t2, c2 = T.coarsen(tr[k], comp)
        pool.append(dict(tag=tag, array=k[0], unit=k[1], chrom=r['chrom'], start=r['start0'], end=r['end'], profile=r['segment_profile'],
                         si=r['var_structure_index'], si_z=r['var_structure_z'], annot=r['annotation_frac'], tr=t2, comp=c2))
n = len(pool); D = [[0.0] * n for _ in range(n)]
for i in range(n):
    for j in range(i + 1, n):
        dt = T.D_trans(pool[i]['tr'], pool[i]['comp'], pool[j]['tr'], pool[j]['comp']); dt = 1.0 if dt == 'NA' else dt
        D[i][j] = D[j][i] = dt + T.jsd(pool[i]['comp'], pool[j]['comp'])
nn = sorted(min(D[i][j] for j in range(n) if j != i) for i in range(n))
cuts = [round(nn[min(n - 1, int(q * n))], 4) for q in QS]
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
    w.writerow(['subgenome', 'array_id', 'unit', 'chrom', 'start0', 'end', 'segment_profile', 'var_structure_index', 'annotation_frac'] + [f'type_q{q}' for q in QS])
    for i, p in enumerate(pool): w.writerow([p['tag'], p['array'], p['unit'], p['chrom'], p['start'], p['end'], p['profile'], p['si'], p['annot']] + [l[i] for l in labs])
KEYS = [('V1', 'V1'), ('V1', 'V2'), ('V2', 'V1'), ('V2', 'V2'), ('V1', 'U'), ('U', 'U')]
num = lambda v: float(v) if v not in ('NA', '') else None
med = lambda v: sorted(v)[len(v) // 2] if v else 'NA'
with open(os.path.join(DEST, '53c_library_types.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['quantile', 'cut_D', 'type', 'n_segments', 'n_Gar', 'n_Mac', 'sharing', 'bp_total'] + [f'P({b}|{a})' for a, b in KEYS] +
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
            w.writerow([q, cut, ty, len(m), g.get('Gar', 0), g.get('Mac', 0), sh, sum(int(p['end']) - int(p['start']) for p in m)] +
                       [round(tr[k] / row[k[0]], 3) if row[k[0]] else 'NA' for k in KEYS] +
                       [round(comp['V1'] / tot, 3), round(comp['V2'] / tot, 3), round(1 - sum(comp[s] for s in ('V1', 'V2', 'V3+', 'U', 'OF')) / tot, 3),
                        round(med(si), 3) if si else 'NA', ';'.join(f'{k}:{v}' for k, v in collections.Counter(p['profile'] for p in m).most_common()),
                        ','.join(f"{p['tag']}:{p['array']}:{p['unit']}" for p in m)])
print(n, 'segments; nearest-neighbour D median', round(nn[n // 2], 3), '; cuts', cuts, '; types', [len(set(l.values())) for l in labs])
for q, lab in zip(QS, labs):
    sh = collections.Counter()
    for ty in set(lab.values()):
        g = {pool[i]['tag'] for i in lab if lab[i] == ty}; sh['shared' if len(g) > 1 else 'single_subgenome'] += 1
    print(f'  q{q}: {dict(sh)}')
