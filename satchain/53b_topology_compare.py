#!/usr/bin/env python3
"""
53b_topology_compare.py  (v4.5.2; big topology vs subtopologies, with annotation states; read-only)

PATH. Each array becomes one ordered path of states along the chromosome:
  copy states   confident step-03 copies (identity >= 80, complete), classes merged by 53_class_merges.tsv, then named by RANK
                inside the array (V1 = most copies, V2, V3, Vx = the rest) so paths of different arrays share one state space;
                a copy of another family than the array's dominant one = OF; an uncertain copy = U
  annotation    TE:<superfamily> (step 45 insertions), 5S, rRNA (other barrnap), tRNA, SINE (AnnoSINE), miRNA, CDS, exon (lifton),
  states        GAP (detgaps). An annotation covering >= 50 % of a copy replaces it; consecutive copies replaced by the same
                feature collapse into ONE state (a TE is one step of the path, never masked). Annotation features lying
                between copies are inserted in order. Two neighbours farther apart than --max-gap-frac x period with no
                feature in between = path break (no transition across it).
UNITS. The whole array ('array' = big topology) and each step-53 segment ('segment' = subtopology).
PER UNIT. composition pi; transition counts C[i][j]; entropy of pi; entropy rate H(next | current); memory MI = H(next) - H(next | current);
  structure index S = (stay_obs - stay_exp) / (1 - stay_exp), stay_exp = sum_i r_i c_i / N^2 (r = row, c = column totals):
  +1 = pure blocks, 0 = no memory (independent draws), < 0 = alternation; z of stay_obs against --perms shuffles of the
  unit's own path (composition kept, order destroyed).
COMPARISON of two units P, Q (same state space):
  D_trans = sum_i w_i JSD(P_i., Q_i.), w_i = mean of the two stationary frequencies (rows with no data in either are skipped)
            -> 'do the same states lead to the same next states'
  D_comp  = JSD(pi_P, pi_Q) -> 'are the same states present in the same amounts'
  Both are in bits, 0 = identical, 1 = disjoint. Reported: every segment vs its whole array, vs the next segment, and all
  segment pairs across arrays (53b_segment_pairs.tsv) for clustering into topology types.
Outputs (<53 dir>/53b/): 53b_topologies.tsv, 53b_transitions.tsv (long format counts per unit), 53b_segment_pairs.tsv,
  53b_state_map.tsv (rank name -> step-03 classes per array).
NOT here: second-order (Miga) topology, i.e. transitions conditioned on the two previous states / cycles; planned next.
usage: 53b_topology_compare.py --out <pipeline out of one subgenome> --run <53 dir> --func <func annotation genome folder> [--arrays ...]
"""
import argparse, collections, csv, gzip, math, os, random

csv.field_size_limit(1 << 30)


def opn(p): return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)


def H(c):
    n = sum(c.values()); return -sum(v / n * math.log2(v / n) for v in c.values() if v) if n else 0.0


def jsd(p, q):
    keys = set(p) | set(q); sp, sq = sum(p.values()) or 1, sum(q.values()) or 1
    P = {k: p.get(k, 0) / sp for k in keys}; Q = {k: q.get(k, 0) / sq for k in keys}; M = {k: (P[k] + Q[k]) / 2 for k in keys}
    kl = lambda A: sum(A[k] * math.log2(A[k] / M[k]) for k in keys if A[k] > 0)
    return (kl(P) + kl(Q)) / 2


def te_label(s):
    s = s.split('|', 1)[-1].split(',')[0].rstrip('/')
    for key, lab in (('Gypsy', 'TE:Gypsy'), ('Copia', 'TE:Copia'), ('LTR', 'TE:LTR'), ('LINE', 'TE:LINE'), ('SINE', 'TE:SINE'), ('Helitron', 'TE:Helitron'),
                     ('DNA', 'TE:DNA'), ('TIR', 'TE:DNA'), ('MITE', 'TE:DNA'), ('fragment', 'TE:fragment')):
        if key.lower() in s.lower(): return lab
    return 'TE:other'


EDTA_SKIP = {'repeat_region', 'target_site_duplication', 'long_terminal_repeat'}


def load_annotations(a, chroms):
    feats = collections.defaultdict(list)
    def add(c, s, e, lab, fid):
        if c in chroms and e > s: feats[c].append((s, e, lab, fid))
    F = lambda *p: os.path.join(a.func, *p)
    G = a.genome_tag
    edta = F('_repeats_EDTA', f'fClaHyb_{G}.TEanno.sorted.gff3.gz')
    if os.path.exists(edta):                      # full EDTA annotation (step 45 insertions are derived from it, so not added twice)
        for l in opn(edta):
            if l.startswith('#'): continue
            f = l.split('\t', 9)
            if len(f) >= 9 and f[0] in chroms and f[2] not in EDTA_SKIP: add(f[0], int(f[3]) - 1, int(f[4]), te_label(f[2]), f'e{f[0]}:{f[3]}')
    else:
        for l in opn(os.path.join(a.out, '45_te_insertions.bed')):
            f = l.rstrip('\n').split('\t')
            if len(f) >= 4 and not l.startswith(('#', 'track')): add(f[0], int(f[1]), int(f[2]), te_label(f[3]), f'te{f[0]}:{f[1]}')
    rf = F('_ncRNA_rfam', f'fClaHyb_{G}.rfam.all_classes.sorted.gff3.gz')
    if os.path.exists(rf):                        # Rfam classes other than rRNA/tRNA/miRNA (those come from their own tools)
        for l in opn(rf):
            if l.startswith('#'): continue
            f = l.rstrip('\n').split('\t')
            if len(f) < 9 or f[0] not in chroms: continue
            nm = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x).get('Name', f[2])
            if any(k in nm for k in ('rRNA', 'tRNA', 'mir', 'MIR', 'let-7')): continue
            add(f[0], int(f[3]) - 1, int(f[4]), 'snRNA' if f[2] == 'snRNA' or nm.startswith('U') else 'ncRNA:' + f[2], f'f{f[0]}:{f[3]}')
    for nm in ('di', 'tri', 'tetra', 'penta_hexa'):   # microsatellites inside arrays
        p = F('_TR', 'GMATA_SSR', f'fClaHyb_{G}.SSR.{nm}.bed.gz')
        if os.path.exists(p):
            for l in opn(p):
                f = l.split('\t')
                if len(f) >= 3 and f[0] in chroms: add(f[0], int(f[1]), int(f[2]), 'SSR', f'ssr{f[0]}:{f[1]}')
    for path, kind in ((F('_rRNA_barrnap', f'fClaHyb_{G}.rRNA.sorted.gff3.gz'), 'rRNA'), (F('_tRNA_tRNAscan', f'fClaHyb_{G}.tRNA.sorted.gff3.gz'), 'tRNA'),
                       (F('_miRNA', f'fClaHyb_{G}.miRNA.curated.sorted.gff3.gz'), 'miRNA'), (F('_genes_lifton', f'fClaHyb_{G}.lifton.sorted.gff3.gz'), 'gene')):
        if not os.path.exists(path):
            alt = os.path.join(a.overlays, os.path.basename(path))
            if os.path.exists(alt): path = alt
            else: print('missing', path); continue
        for l in opn(path):
            if l.startswith('#'): continue
            f = l.rstrip('\n').split('\t')
            if len(f) < 9 or f[0] not in chroms: continue
            s, e = int(f[3]) - 1, int(f[4])
            if kind == 'rRNA': add(f[0], s, e, '5S' if '5S' in f[8] else 'rRNA', f'r{f[0]}:{s}')
            elif kind == 'tRNA' and f[2] in ('tRNA', 'gene', 'pseudogene'): add(f[0], s, e, 'tRNA', f't{f[0]}:{s}')
            elif kind == 'miRNA' and f[2] in ('miRNA', 'pre_miRNA', 'miRNA_primary_transcript', 'gene'): add(f[0], s, e, 'miRNA', f'm{f[0]}:{s}')
            elif kind == 'gene' and f[2] in ('CDS', 'exon'): add(f[0], s, e, 'CDS' if f[2] == 'CDS' else 'exon', f'g{f[2]}{f[0]}:{s}')
    sine = os.path.join(a.overlays, f'fClaHyb_{G}.SINE.AnnoSINE.merged.bed.gz')
    if os.path.exists(sine):
        for l in opn(sine):
            f = l.split('\t')
            if f[0] in chroms: add(f[0], int(f[1]), int(f[2]), 'SINE', f's{f[0]}:{f[1]}')
    gp = F('_genome', f'fClaHyb_{G}_LG.gaps.detgaps.bed')
    if os.path.exists(gp):
        for l in open(gp):
            f = l.split()
            if len(f) >= 3 and f[0] in chroms: add(f[0], int(f[1]), int(f[2]), 'GAP', f'gap{f[0]}:{f[1]}')
    # exon/CDS: keep CDS where both overlap the same span
    for c in feats: feats[c].sort()
    return feats


PRIO = {'GAP': 0, 'TE': 1, '5S': 2, 'rRNA': 2, 'tRNA': 3, 'miRNA': 3, 'snRNA': 3, 'ncRNA': 4, 'SINE': 4, 'CDS': 5, 'exon': 6, 'SSR': 7}
prio = lambda lab: PRIO.get(lab.split(':')[0], 9)


def ins_label(L):
    return 'INS:<1kb' if L < 1000 else 'INS:1-5kb' if L < 5000 else 'INS:>5kb'


def build_path(cs, feats, period, gapfrac):
    """cs: copies sorted (s, e, state). feats: [(s, e, lab, fid)] overlapping the array, sorted.
    A feature covering >= 50 % of a copy replaces it; other features are inserted in order; consecutive pieces of one feature, and
    consecutive overlapping features with the same label, collapse into one state. An uncovered stretch between neighbours longer than
    max(50 bp, gapfrac x period) becomes an INS node (non-satellite, unannotated DNA). An assembly gap is a path break (None)."""
    import bisect
    fst = [f[0] for f in feats]; maxlen = max((f[1] - f[0] for f in feats), default=0)
    def near(s, e):
        i0 = bisect.bisect_left(fst, s - maxlen); i1 = bisect.bisect_left(fst, e)
        return [feats[i] for i in range(i0, i1) if feats[i][1] > s]
    items = []; used = set()
    for s, e, st in cs:
        best = None
        for fs, fe, lab, fid in near(s, e):
            if min(e, fe) - max(s, fs) >= 0.5 * (e - s) and (best is None or prio(lab) < prio(best[0])): best = (lab, fid)
        if best: items.append((s, e, best[0], best[1])); used.add(best[1])
        else: items.append((s, e, st, None))
    kept = [i for i, it in enumerate(items) if it[3] is None]; kst = [items[i][0] for i in kept]; tags = collections.defaultdict(set)
    for fs, fe, lab, fid in feats:
        if fid in used: continue
        j = bisect.bisect_left(kst, fe); host = None
        for i in kept[max(0, j - 50):j]:
            s, e = items[i][0], items[i][1]
            if e > fs and min(e, fe) - max(s, fs) >= 0.5 * (fe - fs): host = i; break
        if host is not None and lab != 'GAP':      # feature inside a copy: the copy carries it (cross feature), e.g. V1.SINE, V1.SSR
            if items[host][2] not in ('U', 'OF'): tags[host].add(lab.split(':')[0] if lab in ('CDS', 'exon') else lab)
        elif lab != 'exon': items.append((fs, fe, lab, fid))
    for i, tg in tags.items():
        s, e, st, fid = items[i]; items[i] = (s, e, st + ''.join(f'.{x}' for x in sorted(tg)), fid)
    items.sort(key=lambda x: (x[0], x[1]))
    path = []; last = None; last_end = None
    for s, e, st, fid in items:
        if last is not None and ((fid is not None and fid == last[1]) or
                                 (fid is not None and last[1] is not None and st == last[0] and s < last_end)):
            last_end = max(last_end, e); continue
        if last_end is not None and s - last_end > max(50, gapfrac * max(1, period)): path.append((last_end, ins_label(s - last_end)))
        path.append((s, None if st == 'GAP' else st)); last = (st, fid); last_end = e if last_end is None else max(last_end, e)
    return path


def order2(path, rng, sims=100):
    """second-order (Miga) test on a path: H(next | previous) - H(next | previous two), in bits, against simulations of the path's own
    first-order chain (same run lengths, same transition matrix): excess and z; plus the most informative two-state context"""
    runs = []; cur = []
    for _, s in path:
        if s is None:
            if len(cur) > 2: runs.append(cur)
            cur = []
        else: cur.append(s)
    if len(cur) > 2: runs.append(cur)
    def gain(rr):
        p1 = collections.defaultdict(collections.Counter); p2 = collections.defaultdict(collections.Counter)
        for r in rr:
            for a, b in zip(r, r[1:]): p1[a][b] += 1
            for a, b, c in zip(r, r[1:], r[2:]): p2[(a, b)][c] += 1
        n1 = sum(sum(v.values()) for v in p1.values()); n2 = sum(sum(v.values()) for v in p2.values())
        if not n2: return None, p1, p2
        h1 = sum(sum(v.values()) / n1 * H(v) for v in p1.values()); h2 = sum(sum(v.values()) / n2 * H(v) for v in p2.values())
        return h1 - h2, p1, p2
    g, p1, p2 = gain(runs)
    if g is None: return dict(order2_gain_bits='NA', order2_null_bits='NA', order2_excess_bits='NA', order2_z='NA', order2_top_context='NA')
    P = {a: (list(v), [x / sum(v.values()) for x in v.values()]) for a, v in p1.items()}
    null = []
    for _ in range(sims):
        rr = []
        for r in runs:
            s = [r[0]]
            for _ in range(len(r) - 1):
                nx = P.get(s[-1])
                s.append(rng.choices(nx[0], nx[1])[0] if nx else s[-1])
            rr.append(s)
        gn, _, _ = gain(rr); null.append(gn or 0.0)
    mu = sum(null) / len(null); sd = (sum((x - mu) ** 2 for x in null) / max(1, len(null) - 1)) ** 0.5
    best = None
    for (a, b), v in p2.items():
        n = sum(v.values())
        if n < 10 or b not in P: continue
        q = dict(zip(*P[b])); kl = sum(c / n * math.log2((c / n) / max(1e-9, q.get(x, 1e-9))) for x, c in v.items())
        if best is None or kl * n > best[0]: best = (kl * n, f'{a}>{b}:{v.most_common(1)[0][0]}({v.most_common(1)[0][1] / n:.2f} vs {q.get(v.most_common(1)[0][0], 0):.2f})')
    return dict(order2_gain_bits=round(g, 4), order2_null_bits=round(mu, 4), order2_excess_bits=round(g - mu, 4),
                order2_z=round((g - mu) / sd, 2) if sd > 1e-6 else 'NA', order2_top_context=best[1] if best else 'none')



def stats(path, rng, perms):
    seq = [st for _, st in path]
    trans = collections.Counter(); comp = collections.Counter(st for st in seq if st is not None)
    runs = []; cur = []
    for st in seq:
        if st is None:
            if cur: runs.append(cur)
            cur = []
        else: cur.append(st)
    if cur: runs.append(cur)
    for r in runs:
        for x, y in zip(r, r[1:]): trans[(x, y)] += 1
    N = sum(trans.values())
    if N == 0: return dict(n_states=sum(comp.values()), n_trans=0), trans, comp
    row, col = collections.Counter(), collections.Counter()
    for (x, y), v in trans.items(): row[x] += v; col[y] += v
    stay = sum(v for (x, y), v in trans.items() if x == y) / N
    exp = sum(row[i] * col[i] for i in row) / N ** 2
    Hr = sum(row[i] / N * H(collections.Counter({y: v for (x, y), v in trans.items() if x == i})) for i in row)
    mi = H(col) - Hr
    def mi_of(tc):
        rw, cl = collections.Counter(), collections.Counter()
        for (x, y), v in tc.items(): rw[x] += v; cl[y] += v
        n_ = sum(tc.values())
        return H(cl) - sum(rw[i] / n_ * H(collections.Counter({y: v for (x, y), v in tc.items() if x == i})) for i in rw)
    flat = [s for r in runs for s in r]; nulls = []; mnull = []
    for _ in range(perms):
        rng.shuffle(flat); it = iter(flat); st_ = 0; tc = collections.Counter()
        for r in runs:
            rr = [next(it) for _ in r]; st_ += sum(1 for x, y in zip(rr, rr[1:]) if x == y)
            for x, y in zip(rr, rr[1:]): tc[(x, y)] += 1
        nulls.append(st_ / N); mnull.append(mi_of(tc))
    mu = sum(nulls) / len(nulls) if nulls else exp; sd = (sum((x - mu) ** 2 for x in nulls) / max(1, len(nulls) - 1)) ** 0.5 if nulls else 0
    return dict(n_states=sum(comp.values()), n_trans=N, n_breaks=len(runs) - 1, entropy_bits=round(H(comp), 3), entropy_rate_bits=round(Hr, 3),
                memory_MI_bits=round(mi, 3), memory_null_bits=round(sum(mnull) / len(mnull), 3) if mnull else 'NA',
                memory_excess_bits=round(mi - sum(mnull) / len(mnull), 3) if mnull else 'NA',
                memory_z=round((mi - sum(mnull) / len(mnull)) / ((sum((x - sum(mnull) / len(mnull)) ** 2 for x in mnull) / max(1, len(mnull) - 1)) ** 0.5), 2) if len(mnull) > 1 and len(set(round(x, 9) for x in mnull)) > 1 else 'NA',
                size_flag='small_n' if N < 40 else 'ok',
                stay_obs=round(stay, 3), stay_expected=round(exp, 3),
                structure_index=round((stay - exp) / (1 - exp), 3) if exp < 1 else 'NA',
                structure_z=round((stay - mu) / sd, 2) if sd > 1e-3 and exp < 0.95 else 'NA',
                annotation_frac=round(sum(v for k, v in comp.items() if not (k.split('.')[0].startswith('V') or k in ('U', 'OF'))) / max(1, sum(comp.values())), 4),
                composition=';'.join(f'{k}:{v}' for k, v in comp.most_common())), trans, comp


def coarse(s):
    """library alphabet: variants by rank with tags folded in; annotation states grouped"""
    b = s.split('.')[0]
    if b in ('V1', 'V2'): return b
    if b.startswith('V'): return 'V3+'
    if b in ('U', 'OF', 'SSR', 'SINE'): return b
    if b.startswith('TE'): return 'TE'
    if b.startswith('INS'): return 'INS'
    if b in ('CDS', 'exon'): return 'gene'
    if b in ('5S', 'rRNA', 'tRNA', 'snRNA', 'miRNA') or b.startswith('ncRNA'): return 'ncRNA'
    return 'other'


def coarsen(tr, comp):
    t2, c2 = collections.Counter(), collections.Counter()
    for (x, y), v in tr.items(): t2[(coarse(x), coarse(y))] += v
    for k, v in comp.items(): c2[coarse(k)] += v
    return t2, c2


def D_trans(t1, c1, t2, c2):
    rows = {x for x, _ in t1} | {x for x, _ in t2}; n1, n2 = sum(c1.values()) or 1, sum(c2.values()) or 1; d = 0.0; wsum = 0.0
    for i in rows:
        a = collections.Counter({y: v for (x, y), v in t1.items() if x == i}); b = collections.Counter({y: v for (x, y), v in t2.items() if x == i})
        if not a or not b: continue
        w = (c1.get(i, 0) / n1 + c2.get(i, 0) / n2) / 2; d += w * jsd(a, b); wsum += w
    return round(d / wsum, 4) if wsum else 'NA'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True); ap.add_argument('--run', required=True); ap.add_argument('--func', required=True)
    ap.add_argument('--overlays', required=True); ap.add_argument('--genome-tag', required=True)
    ap.add_argument('--arrays'); ap.add_argument('--max-gap-frac', type=float, default=0.5); ap.add_argument('--perms', type=int, default=200)
    ap.add_argument('--runs', help='catalogue runs tsv (run mode, as in 53 --runs)')
    ap.add_argument('--order2-sims', type=int, default=100, help='first-order simulations for the second-order (Miga) test')
    ap.add_argument('--min-trans', type=int, default=40, help='segments with fewer transitions are reported but not paired')
    ap.add_argument('--intrinsic-frac', type=float, default=0.3)
    ap.add_argument('--uniform-frac', type=float, default=0.95); ap.add_argument('--copy-scale', type=float, default=2.0)
    ap.add_argument('--major-p', type=float, default=0.15); ap.add_argument('--major-n', type=int, default=5)
    ap.add_argument('--lib-quantiles', type=lambda s: [float(x) for x in s.split(',')], default=[0.25, 0.5, 0.75, 0.9],
                    help='library cut = this quantile of the nearest-neighbour distance (sweep)')
    a = ap.parse_args()
    dest = os.path.join(a.run, '53b'); os.makedirs(dest, exist_ok=True)
    R = lambda f: list(csv.DictReader(open(os.path.join(a.run, f)), delimiter='\t'))
    segs = collections.defaultdict(list)
    for r in R('53_segments.tsv'): segs[r['array_id']].append(r)
    want = set(segs)
    if a.arrays: want &= set(open(a.arrays).read().split()) if os.path.exists(a.arrays) else set(a.arrays.split(','))
    if a.runs:                                   # run mode: each catalogue run is a unit; copies = the parent array's copies inside it
        arrays = {r['run_id']: dict(array_id=r['run_id'], parent=r['array_id'], chrom=r['chrom'], start_0based=r['start0'], end=r['end'], period_mode_bp=r['unit_bp'])
                  for r in csv.DictReader(open(a.runs), delimiter='\t') if r['run_id'] in want}
    else:
        arrays = {r['array_id']: r for r in csv.DictReader(open(os.path.join(a.out, '02_satellite_arrays.tsv')), delimiter='\t') if r['array_id'] in want}
    parents = {A.get('parent', k) for k, A in arrays.items()}
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x: x = parent[x]
        return x
    for r in R('53_class_merges.tsv'):
        if r['decision'] == 'merge': parent[find((r['array_id'], r['class_b']))] = find((r['array_id'], r['class_a']))
    cop = collections.defaultdict(list)
    with open(os.path.join(a.out, '03_monomer_copies.tsv'), newline='') as fh:
        rd = csv.reader(fh, delimiter='\t'); h = next(rd); ix = {k: i for i, k in enumerate(h)}
        for row in rd:
            aid = row[ix['array_id']]
            if aid not in parents: continue
            try: idt = float(row[ix['best_alignment_identity']])
            except ValueError: idt = 0
            conf = row[ix['complete_or_truncated']] == 'complete' and idt >= 80
            cop[aid].append((int(row[ix['genomic_start_0based']]), int(row[ix['genomic_end']]), row[ix['monomer_class']] if conf else 'U', row[ix['family_id']]))
    feats = load_annotations(a, {r['chrom'] for r in arrays.values()})
    rng = random.Random(53)
    tw = csv.writer(open(os.path.join(dest, '53b_topologies.tsv'), 'w', newline=''), delimiter='\t', lineterminator='\n')
    cols = ['n_states', 'n_trans', 'n_breaks', 'size_flag', 'entropy_bits', 'entropy_rate_bits', 'memory_MI_bits', 'memory_null_bits', 'memory_excess_bits', 'memory_z', 'stay_obs', 'stay_expected', 'structure_index', 'structure_z', 'annotation_frac', 'var_n_trans', 'var_memory_excess_bits', 'var_memory_z', 'var_structure_index', 'var_structure_z', 'order2_gain_bits', 'order2_null_bits', 'order2_excess_bits', 'order2_z', 'order2_top_context']
    tw.writerow(['array_id', 'chrom', 'unit', 'start0', 'end', 'segment_profile'] + cols + ['D_trans_vs_array', 'D_comp_vs_array', 'D_trans_vs_next', 'D_comp_vs_next', 'composition'])
    lw = csv.writer(open(os.path.join(dest, '53b_transitions.tsv'), 'w', newline=''), delimiter='\t', lineterminator='\n')
    lw.writerow(['array_id', 'unit', 'from', 'to', 'count'])
    mw = csv.writer(open(os.path.join(dest, '53b_state_map.tsv'), 'w', newline=''), delimiter='\t', lineterminator='\n')
    mw.writerow(['array_id', 'state', 'step03_classes', 'n_copies'])
    cw = csv.writer(open(os.path.join(dest, '53b_core_accessory_edges.tsv'), 'w', newline=''), delimiter='\t', lineterminator='\n')
    cw.writerow(['array_id', 'edge', 'role', 'n_segments_with_edge', 'n_segments_tested', 'segments', 'intrinsic_annotations'])
    pool = []
    for aid, A in arrays.items():
        chrom = A['chrom']; s0, e0 = int(A['start_0based']), int(A['end']); period = float(A['period_mode_bp'] or 0)
        pid = A.get('parent', aid); cs = sorted(c for c in cop.get(pid, []) if 'parent' not in A or s0 <= c[0] < e0)
        domfam = collections.Counter(f for *_, c, f in [(x[0], x[1], x[2], x[3]) for x in cs] if c != 'U').most_common(1)
        domfam = domfam[0][0] if domfam else None
        groups = collections.Counter(find((pid, c)) for _, _, c, f in cs if c != 'U' and f == domfam)
        rank = {g: (f'V{i + 1}' if i < 3 else 'Vx') for i, (g, _) in enumerate(groups.most_common())}
        members = collections.defaultdict(set)
        for _, _, c, f in cs:
            if c != 'U' and f == domfam: members[rank[find((pid, c))]].add(c)
        for st_, m in sorted(members.items()): mw.writerow([aid, st_, '+'.join(sorted(m)), sum(groups[g] for g, r in rank.items() if r == st_)])
        states = [(s, e, 'U' if c == 'U' else ('OF' if f != domfam else rank[find((pid, c))])) for s, e, c, f in cs]
        F = [x for x in feats.get(chrom, []) if x[1] > s0 and x[0] < e0]
        # an annotation class overlapping > --intrinsic-frac of the array's copies labels the satellite itself (e.g. AnnoSINE on a
        # 5S-related satellite): recorded as an attribute of the array, not inserted as path states
        # intrinsic annotations: a label carried by > --intrinsic-frac of the array's confident copies (carried = one feature of that
        # label covers >= 50 % of the copy) describes the satellite itself, so it is never inserted as path states. Three kinds:
        #   uniform   carried by >= --uniform-frac of the copies (e.g. EDTA calling the whole satellite LTR): array attribute only
        #   regional  median feature length > --copy-scale x period (one feature spans many copies): region attribute only, since
        #             tagging copies with it would create blocks by construction
        #   copy      copy-scale and partial (e.g. per-monomer AnnoSINE calls): becomes a cross-feature node, V1.SINE vs V1
        import bisect
        nconf = max(1, sum(1 for _, _, c, _ in cs if c != 'U')); intrinsic = {}; cross = []
        byl = collections.defaultdict(list)
        for x in F: byl[x[2]].append((x[0], x[1]))
        def carries(iv, ivs, ml, s, e):
            i0 = bisect.bisect_left(ivs, s - ml); i1 = bisect.bisect_left(ivs, e)
            return any(min(e, iv[i][1]) - max(s, iv[i][0]) >= 0.5 * min(e - s, iv[i][1] - iv[i][0]) for i in range(i0, i1))
        idx = {}
        for lab, iv in byl.items():
            if lab == 'GAP': continue
            iv.sort(); ivs = [x[0] for x in iv]; ml = max(b - a_ for a_, b in iv); idx[lab] = (iv, ivs, ml)
            hit = sum(1 for s, e, c, _ in cs if c != 'U' and carries(iv, ivs, ml, s, e)); fr = hit / nconf
            if fr > a.intrinsic_frac:
                medlen = sorted(b - a_ for a_, b in iv)[len(iv) // 2]
                kind = 'uniform' if fr >= a.uniform_frac else 'regional' if medlen > a.copy_scale * max(1, period) else 'copy'
                intrinsic[lab] = f'{fr:.2f}/{kind}'
                if kind == 'copy': cross.append(lab)
        F = [x for x in F if x[2] not in intrinsic or x[2] in cross]   # copy-scale intrinsic labels stay, to tag their copies
        path = build_path(states, F, period, a.max_gap_frac)
        units = [('array', s0, e0, 'NA')] + [(f"seg{r['segment']}", int(r['start0']), int(r['end']), r['profile']) for r in segs[aid]]
        res = []
        for name, us, ue, prof in units:
            p = [(x, st) for x, st in path if us <= x < ue]
            st, tr, comp = stats(p, rng, a.perms)
            # variant-only path: non-variant states become breaks, so order among monomer variants is measured alone
            pv = [(x, s_ if (s_ is not None and s_.split('.')[0].startswith('V')) else None) for x, s_ in p]
            sv, _, _ = stats(pv, rng, a.perms)
            for k_ in ('n_trans', 'memory_excess_bits', 'memory_z', 'structure_index', 'structure_z'): st['var_' + k_] = sv.get(k_, 'NA')
            st.update(order2(pv, rng, a.order2_sims))
            res.append((name, us, ue, prof, st, tr, comp))
            for (x, y), v in sorted(tr.items()): lw.writerow([aid, name, x, y, v])
        arr = res[0]
        major = {}
        for name, us, ue, prof, st, tr, comp in res[1:]:
            if st.get('n_trans', 0) < a.min_trans or prof in ('unresolved', 'too_short'): continue   # flanks / tiny runs are not subtopologies
            row = collections.Counter()
            for (x, y), v in tr.items(): row[x] += v
            major[name] = {(x, y) for (x, y), v in tr.items() if v >= a.major_n and v / row[x] >= a.major_p}
        if major:
            core = set.intersection(*major.values()); allm = set.union(*major.values())
            for e_ in sorted(allm):
                where = [n for n, m in major.items() if e_ in m]
                cw.writerow([aid, f'{e_[0]}>{e_[1]}', 'core' if e_ in core else 'accessory', len(where), len(major), ','.join(where),
                             ';'.join(f'{k}:{v}' for k, v in intrinsic.items()) or 'none'])
        for i, (name, us, ue, prof, st, tr, comp) in enumerate(res):
            dva = (D_trans(tr, comp, arr[5], arr[6]), round(jsd(comp, arr[6]), 4)) if i else ('NA', 'NA')
            nxt = res[i + 1] if i and i + 1 < len(res) else None
            dvn = (D_trans(tr, comp, nxt[5], nxt[6]), round(jsd(comp, nxt[6]), 4)) if nxt else ('NA', 'NA')
            tw.writerow([aid, chrom, name, us, ue, prof] + [st.get(c, 'NA') for c in cols] + list(dva) + list(dvn) + [st.get('composition', '')])
            if i and st.get('n_trans', 0) >= a.min_trans and prof not in ('unresolved', 'too_short'): pool.append((aid, name, prof) + coarsen(tr, comp))
    pw = csv.writer(open(os.path.join(dest, '53b_segment_pairs.tsv'), 'w', newline=''), delimiter='\t', lineterminator='\n')
    pw.writerow(['array_a', 'unit_a', 'profile_a', 'array_b', 'unit_b', 'profile_b', 'D_trans', 'D_comp'])
    for i in range(len(pool)):
        for j in range(i + 1, len(pool)):
            x, y = pool[i], pool[j]
            pw.writerow([x[0], x[1], x[2], y[0], y[1], y[2], D_trans(x[3], x[4], y[3], y[4]), round(jsd(x[4], y[4]), 4)])
    n = len(pool); D = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            dt = D_trans(pool[i][3], pool[i][4], pool[j][3], pool[j][4]); dt = 1.0 if dt == 'NA' else dt
            D[i][j] = D[j][i] = dt + jsd(pool[i][4], pool[j][4])
    nn = sorted(min(D[i][j] for j in range(n) if j != i) for i in range(n)) if n > 1 else [0]
    qcut = lambda q: nn[min(len(nn) - 1, int(q * len(nn)))]
    cuts = [round(qcut(q), 4) for q in a.lib_quantiles]
    lw2 = csv.writer(open(os.path.join(dest, '53b_topology_library.tsv'), 'w', newline=''), delimiter='\t', lineterminator='\n')
    lw2.writerow(['array_id', 'unit', 'segment_profile'] + [f'type_at_q{q}_D{c}' for q, c in zip(a.lib_quantiles, cuts)])
    labels = []
    for cut in cuts:
        cl = [[i] for i in range(n)]
        while len(cl) > 1:
            best = None
            for x in range(len(cl)):
                for y in range(x + 1, len(cl)):
                    d = sum(D[i][j] for i in cl[x] for j in cl[y]) / (len(cl[x]) * len(cl[y]))
                    if best is None or d < best[0]: best = (d, x, y)
            if best[0] > cut: break
            _, x, y = best; cl[x] += cl[y]; del cl[y]
        lab = {}
        for k, c in enumerate(sorted(cl, key=len, reverse=True)):
            for i in c: lab[i] = f'T{k + 1:02d}'
        labels.append(lab)
    for i, (aid, name, prof, *_ ) in enumerate(pool): lw2.writerow([aid, name, prof] + [l[i] for l in labels])
    print(f'{len(arrays)} arrays, {len(pool)} segments paired; nearest-neighbour D median {nn[len(nn)//2]:.3f}; cuts', cuts, '-> library types', [len(set(l.values())) for l in labels], '->', dest)


if __name__ == '__main__':
    main()
