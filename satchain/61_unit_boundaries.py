#!/usr/bin/env python3
"""
61_unit_boundaries.py  (SatChain v4.7 draft; per-array unit boundaries to the base pair, with a stability score)

An interrupted array is a chain of catalogue runs of one array_id, separated by interruptions (gap <= --max-gap). Every run edge that
faces an interruption is a junction. Two junction types, read in ARRAY orientation (the strand of the array's own consensus):
  OUT = satellite -> interruption (end of a run),  IN = interruption -> satellite (start of the next run).

1 Edge refinement. The edge region (2 consensus lengths of satellite + --flank bp outside) is aligned locally to the array's dominant
  class consensus repeated 4 times, in the copies' orientation. The refined edge is the last (OUT) / first (IN) satellite base of that
  alignment; phase = consensus position of that base (1..L); delta_hmm = refined - catalogue edge. Using the array's own consensus avoids
  scoring a chromosome-specific variant against a distant family consensus.
2 Junction homology. Context = --inner bp of satellite + --flank bp of interruption, in array orientation. Per array and junction type
  the medoid context (highest mean local score to the others; <= 40 sampled) is the reference. Each context is aligned to it (global,
  free end gaps). flank_identity = identity over the medoid's flank columns; offset = position of the junction's edge mapped into the
  medoid minus the medoid edge. A junction is homologous when flank_identity >= cutoff.
3 Stability score per array and setting (flank x cutoff x tol): share of homologous junctions, share of homologous junctions with
  |offset| <= tol, share with phase within +-3 of the modal phase (circular), median |delta_hmm|.
4 Refined units: from one refined IN junction to the next in array orientation (unit = run + following interruption), BED6.

Outputs (<dest>/): 61_junctions_<tag>.tsv (per junction, per flank length), 61_array_scores_<tag>.tsv (per array and setting),
61_units_<tag>.bed, 61_summary_<tag>.tsv.
usage: 61_unit_boundaries.py --catalogue 59_catalogue_<G>.tsv --copies 03_monomer_copies.tsv [...] --consensus 03_array_consensus.tsv [...]
                             --genome <G>.fa --dest <dir> --tag <G> [--min-runs 4] [--max-gap 50000] [--flank 100,150,300] [--inner 60]
                             [--cutoff 0.6,0.7,0.8] [--tol 10,25,50] [--arrays id1,id2]
"""
import argparse, collections, csv, os, random
import numpy as np
from Bio import Align
csv.field_size_limit(1 << 30)
RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def rc(s): return s.translate(RC)[::-1]


class Fasta:
    def __init__(self, path):
        self.fh = open(path, 'rb'); self.idx = {}
        for l in open(path + '.fai'):
            n, ln, off, lb, lby = l.split('\t')[:5]; self.idx[n] = (int(ln), int(off), int(lb), int(lby))
    def get(self, c, a, b):
        ln, off, lb, lby = self.idx[c]; a = max(0, a); b = min(ln, b)
        if b <= a: return ''
        s0 = off + (a // lb) * lby + a % lb; s1 = off + (b // lb) * lby + b % lb
        self.fh.seek(s0); return self.fh.read(s1 - s0).replace(b'\n', b'').replace(b'\r', b'').decode().upper()


LOC = Align.PairwiseAligner(mode='local', match_score=2, mismatch_score=-3, open_gap_score=-5, extend_gap_score=-2)
GLO = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-3, open_gap_score=-5, extend_gap_score=-2,
                            end_open_gap_score=0, end_extend_gap_score=0)


def refine(fa, chrom, edge, side, cons, flank):
    """side 'R' = genomic right edge of a run (satellite on the left), 'L' = genomic left edge. cons already in genomic orientation."""
    L = len(cons); ref = cons * 4
    if side == 'R':
        a = edge - 2 * L; w = fa.get(chrom, a, edge + flank)
    else:
        a = edge - flank; w = fa.get(chrom, a, edge + 2 * L)
    if len(w) < L: return None
    al = LOC.align(w, ref)[0]
    blocks = list(zip(al.aligned[0], al.aligned[1]))
    if sum(e - s for (s, e), _ in blocks) < 0.6 * L: return None
    if side == 'R':
        (s, e), (t0, t1) = blocks[-1]; pos = a + e; phase = (t1 - 1) % L + 1
    else:
        (s, e), (t0, t1) = blocks[0]; pos = a + s; phase = t0 % L + 1
    return pos, phase


def context(fa, chrom, edge, side, rev, inner, flank):
    """array-orientation context; returns (seq, edge index). OUT: satellite then flank; IN: flank then satellite."""
    if side == 'R': g = fa.get(chrom, edge - inner, edge + flank); gi = inner        # satellite | spacer (genomic)
    else: g = fa.get(chrom, edge - flank, edge + inner); gi = flank                  # spacer | satellite (genomic)
    if rev: g = rc(g); gi = len(g) - gi
    return g, gi


def map_edge(al, edge_i):
    """position in target (medoid) of query index edge_i, via aligned blocks; nearest block boundary if edge falls in a gap"""
    best = None
    for (qs, qe), (ts, te) in zip(*al.aligned):
        if qs <= edge_i <= qe: return ts + (edge_i - qs)
        d = min(abs(edge_i - qs), abs(edge_i - qe))
        if best is None or d < best[0]: best = (d, ts if abs(edge_i - qs) <= abs(edge_i - qe) else te)
    return best[1] if best else None


def flank_identity(al, q, t, t_lo, t_hi):
    m = n = 0
    for (qs, qe), (ts, te) in zip(*al.aligned):
        for i in range(qe - qs):
            tp = ts + i
            if t_lo <= tp < t_hi: n += 1; m += q[qs + i] == t[tp]
    span = t_hi - t_lo
    return m / span if span else 0.0, n / span if span else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--catalogue', required=True); ap.add_argument('--copies', nargs='+', required=True)
    ap.add_argument('--consensus', nargs='+', required=True); ap.add_argument('--genome', required=True)
    ap.add_argument('--dest', required=True); ap.add_argument('--tag', required=True)
    ap.add_argument('--min-runs', type=int, default=4); ap.add_argument('--max-gap', type=int, default=50000)
    ap.add_argument('--flank', default='100,150,300'); ap.add_argument('--inner', type=int, default=60)
    ap.add_argument('--cutoff', default='0.6,0.7,0.8'); ap.add_argument('--tol', default='10,25,50')
    ap.add_argument('--arrays', default=''); ap.add_argument('--seed', type=int, default=61)
    ap.add_argument('--gaps', default='', help='assembly-gap BED: chains never cross a gap (v4.7)')
    a = ap.parse_args(); random.seed(a.seed); os.makedirs(a.dest, exist_ok=True)
    FL = [int(x) for x in a.flank.split(',')]; CU = [float(x) for x in a.cutoff.split(',')]; TO = [int(x) for x in a.tol.split(',')]
    keep = set(a.arrays.split(',')) if a.arrays else None
    fa = Fasta(a.genome)
    runs = collections.defaultdict(list)
    for r in csv.DictReader(open(a.catalogue), delimiter='\t'):
        if r.get('element_class', 'satellite') != 'satellite': continue
        if keep and r['array_id'] not in keep: continue
        runs[r['array_id']].append(r)
    # arrays = chains of >= min_runs runs, one family, gaps <= max_gap, never across an assembly gap (v4.7)
    import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from gaps import Gaps
    GAPS = Gaps(a.gaps) if a.gaps else None; refused = 0
    chains = []
    for aid, rs in runs.items():
        rs.sort(key=lambda r: int(r['start0'])); cur = [rs[0]]
        for r in rs[1:]:
            across = GAPS is not None and GAPS.between(r['chrom'], int(cur[-1]['end']), int(r['start0'])); refused += across
            if r['family'] == cur[-1]['family'] and 0 < int(r['start0']) - int(cur[-1]['end']) <= a.max_gap and not across: cur.append(r)
            else:
                if len(cur) >= a.min_runs: chains.append((aid, cur))
                cur = [r]
        if len(cur) >= a.min_runs: chains.append((aid, cur))
    need = {aid for aid, _ in chains}
    print(f'{len(chains)} chains; links refused at assembly gaps: {refused if GAPS is not None else "no gap BED given"}', flush=True)
    # dominant class consensus and copy orientation per array
    cons = {}
    for p in a.consensus:
        for r in csv.DictReader(open(p), delimiter='\t'):
            if r['array_id'] in need:
                k = (r['array_id'], r['monomer_class']); cons[k] = (int(r['n_copies_in_class'] or 0), r['consensus_sequence'].upper())
    cls = collections.defaultdict(collections.Counter); rev = collections.defaultdict(collections.Counter)
    for p in a.copies:
        with open(p, newline='') as fh:
            rd = csv.reader(fh, delimiter='\t'); h = next(rd); ix = {k: i for i, k in enumerate(h)}
            for row in rd:
                aid = row[ix['array_id']]
                if aid not in need: continue
                cls[aid][row[ix['monomer_class']]] += 1; rev[aid][row[ix['reverse_complemented']]] += 1
    J = []; units = []
    for ci, (aid, rs) in enumerate(chains):
        top = cls[aid].most_common(1)
        if not top or (aid, top[0][0]) not in cons: continue
        c0 = cons[(aid, top[0][0])][1]; is_rev = rev[aid]['yes'] > rev[aid]['no']
        cg = rc(c0) if is_rev else c0          # consensus in genomic orientation of the copies
        chain_id = f'{aid}#{ci}'; chrom = rs[0]['chrom']
        for fl in FL:
            ctx = {'OUT': [], 'IN': []}
            for i, (r1, r2) in enumerate(zip(rs, rs[1:])):
                for side, edge, run in (('R', int(r1['end']), r1['run_id']), ('L', int(r2['start0']), r2['run_id'])):
                    ref = refine(fa, chrom, edge, side, cg, fl)
                    jt = ('OUT' if side == 'R' else 'IN') if not is_rev else ('IN' if side == 'R' else 'OUT')
                    if ref is None:
                        J.append(dict(chain=chain_id, chrom=chrom, run_id=run, gap_index=i, type=jt, flank=fl, hmm_edge=edge, refined=None)); continue
                    pos, ph = ref; s, ei = context(fa, chrom, pos, side, is_rev, a.inner, fl)
                    rec = dict(chain=chain_id, array_id=aid, family=rs[0]['family'], chrom=chrom, run_id=run, gap_index=i, type=jt, flank=fl,
                               hmm_edge=edge, refined=pos, delta_hmm=pos - edge, phase=ph, cons_len=len(cg), _s=s, _e=ei)
                    ctx[jt].append(rec); J.append(rec)
            for jt, L_ in ctx.items():
                if len(L_) < 3: continue
                smp = random.sample(L_, min(40, len(L_)))
                med = max(smp, key=lambda x: np.mean([LOC.score(x['_s'], y['_s']) for y in smp if y is not x]))
                t = med['_s']; f_lo, f_hi = (med['_e'], len(t)) if jt == 'OUT' else (0, med['_e'])
                for x in L_:
                    al = GLO.align(x['_s'], t)[0]; mp = map_edge(al, x['_e'])
                    x['offset'] = None if mp is None else mp - med['_e']
                    x['flank_identity'], x['flank_aligned'] = flank_identity(al, x['_s'], t, f_lo, f_hi)
                    x['medoid_run'] = med['run_id']
        # refined units (main flank = first in list): IN junction of run i -> IN junction of run i+1, array orientation
        ins = sorted([x for x in J if x.get('chain') == chain_id and x['flank'] == FL[0] and x['type'] == 'IN' and x.get('refined')], key=lambda x: x['refined'])
        for x, y in zip(ins, ins[1:]):
            units.append((chrom, min(x['refined'], y['refined']), max(x['refined'], y['refined']), f"{chain_id}|{x['run_id']}", '-' if is_rev else '+'))
        if ci % 50 == 0: print(f'{ci}/{len(chains)} arrays', flush=True)
    cols = ['chain', 'array_id', 'family', 'chrom', 'run_id', 'gap_index', 'type', 'flank', 'hmm_edge', 'refined', 'delta_hmm', 'phase', 'cons_len',
            'flank_identity', 'flank_aligned', 'offset', 'medoid_run']
    with open(os.path.join(a.dest, f'61_junctions_{a.tag}.tsv'), 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t'); w.writerow(cols)
        for x in J: w.writerow(['' if x.get(c) is None else x.get(c, '') for c in cols])
    # array scores
    by = collections.defaultdict(list)
    for x in J: by[(x['chain'], x['type'], x['flank'])].append(x)
    S = []
    for (ch, jt, fl), L_ in by.items():
        ok = [x for x in L_ if x.get('refined') is not None and x.get('flank_identity') is not None]
        if len(ok) < 3: continue
        ph = np.array([x['phase'] for x in ok]); L = ok[0]['cons_len']; mode = collections.Counter(ph).most_common(1)[0][0]
        circ = np.minimum((ph - mode) % L, (mode - ph) % L)
        for cu in CU:
            hom = [x for x in ok if x['flank_identity'] >= cu]
            row = dict(chain=ch, array_id=ok[0]['array_id'], family=ok[0]['family'], chrom=ok[0]['chrom'], type=jt, flank=fl, cutoff=cu,
                       junctions=len(L_), refined=len(ok), homologous_share=len(hom) / len(ok), phase_concordance=float(np.mean(circ <= 3)),
                       median_abs_delta_hmm=float(np.median([abs(x['delta_hmm']) for x in ok])))
            for tol in TO:
                row[f'precise_pm{tol}'] = float(np.mean([abs(x['offset']) <= tol for x in hom if x['offset'] is not None])) if hom else np.nan
            S.append(row)
    if S:
        with open(os.path.join(a.dest, f'61_array_scores_{a.tag}.tsv'), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(S[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(S)
    with open(os.path.join(a.dest, f'61_units_{a.tag}.bed'), 'w') as fh:
        for u in units: fh.write(f'{u[0]}\t{u[1]}\t{u[2]}\t{u[3]}\t0\t{u[4]}\n')
    # summary: per setting, median over arrays
    agg = collections.defaultdict(list)
    for r in S: agg[(r['type'], r['flank'], r['cutoff'])].append(r)
    with open(os.path.join(a.dest, f'61_summary_{a.tag}.tsv'), 'w') as fh:
        fh.write('type\tflank\tcutoff\tarrays\tmedian_homologous_share\t' + '\t'.join(f'median_precise_pm{t}' for t in TO) + '\tmedian_phase_concordance\n')
        for (jt, fl, cu), rr in sorted(agg.items()):
            fh.write(f'{jt}\t{fl}\t{cu}\t{len(rr)}\t{np.median([r["homologous_share"] for r in rr]):.3f}\t' +
                     '\t'.join(f'{np.nanmedian([r[f"precise_pm{t}"] for r in rr]):.3f}' for t in TO) +
                     f'\t{np.median([r["phase_concordance"] for r in rr]):.3f}\n')
    print(f'chains {len(chains)}, junctions {len(J)}, units {len(units)}')


if __name__ == '__main__':
    main()
