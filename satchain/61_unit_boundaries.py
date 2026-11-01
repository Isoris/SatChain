#!/usr/bin/env python3
"""
61_unit_boundaries.py  (SatChain v4.7 draft; per-array unit boundaries to the base pair, with a stability score)

An interrupted array is a chain of catalogue runs (chains.py, shared with 61b and 64 since 2026-11-01: one array key = array_id_v46 or
array_id without the 59 '~g' segment suffix, one family_v46, 0 <= gap <= --max-gap, never across an assembly gap, >= --min-runs runs).
Every run edge that faces an interruption is a junction. Two junction types, read in ARRAY orientation (the strand of the array's own
consensus):  OUT = satellite -> interruption (end of a run),  IN = interruption -> satellite (start of the next run).

1 Edge refinement. The edge region (2 consensus lengths of satellite + --flank bp outside) is aligned locally to the array's dominant
  class consensus repeated 4 times, in the copies' orientation. The refined edge is the last (OUT) / first (IN) satellite base of that
  alignment; phase = consensus position of that base (1..L); delta_hmm = refined - catalogue edge. Using the array's own consensus avoids
  scoring a chromosome-specific variant against a distant family consensus. For a '-' array (majority of copies reverse_complemented) the
  phase is a position on the REVERSE-COMPLEMENTED consensus, i.e. on the consensus as it reads along the genome (D14).
2 Junction homology. Context = --inner bp of satellite + --flank bp of interruption, in array orientation. Per array and junction type
  the medoid context (highest mean local score to the others; <= 40 sampled) is the reference. Each context is aligned to it (global,
  free end gaps). flank_identity = identity over the medoid's flank columns; offset = position of the junction's edge mapped into the
  medoid minus the medoid edge. Per-array null (2026-11-01, D06): the same statistic for null contexts that keep the junction's OWN
  satellite side (the anchor; without it unrelated flanks score ~0, with it ~0.5 at 100 bp) and carry (a) --null-shift spacer pieces per
  junction taken flank + inner further into the same spacer and (b) up to --null-cross flanks of the same type and length from OTHER
  chains; homologous_null = flank_identity > null p99 and >= the lowest --cutoff (floor: with 300-bp flanks the free-end-gap alignment
  drops unrelated flanks entirely, p99 ~ 0). The medoid is identical to itself by construction (is_medoid = 1) and enters no share.
  The fixed cutoffs (--cutoff) stay as the stability sweep.
3 Stability score per array and setting (flank x cutoff x tol): share of homologous junctions (cutoff and null-based), share of homologous
  junctions with |offset| <= tol, share with phase within +-3 of the modal phase (circular), median |delta_hmm|.
4 Refined units: from one refined IN junction to the next in array orientation (unit = run + following interruption), BED6 + two columns
  satellite_bp (refined IN -> refined OUT of the run) and interruption_bp (unit - satellite_bp) (2026-11-01, D11).

Outputs (<dest>/): 61_junctions_<tag>.tsv (per junction, per flank length), 61_array_scores_<tag>.tsv (per array and setting),
61_units_<tag>.bed (BED6 + satellite_bp, interruption_bp), 61_summary_<tag>.tsv, 61_skipped_chains_<tag>.tsv (chains without a
consensus / copy class, D12). --sweep writes <repo>/benchmarks/61_sweep.tsv (the per-setting summary; D09).
usage: 61_unit_boundaries.py --catalogue 59_catalogue_<G>.tsv --copies 03_monomer_copies.tsv [...] --consensus 03_array_consensus.tsv [...]
                             --genome <G>.fa --dest <dir> --tag <G> [--min-runs 2] [--max-gap 50000] [--flank 100,150,300] [--inner 60]
                             [--cutoff 0.6,0.7,0.8] [--tol 10,25,50] [--arrays id1,id2] [--gaps gaps.bed] [--null-cross 100] [--sweep]
"""
import argparse, collections, csv, os, random, sys
import numpy as np
from Bio import Align
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gaps import Gaps
from chains import build_chains, MIN_RUNS, MAX_GAP
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
    ap.add_argument('--min-runs', type=int, default=MIN_RUNS); ap.add_argument('--max-gap', type=int, default=MAX_GAP)
    ap.add_argument('--flank', default='100,150,300'); ap.add_argument('--inner', type=int, default=60)
    ap.add_argument('--cutoff', default='0.6,0.7,0.8'); ap.add_argument('--tol', default='10,25,50')
    ap.add_argument('--arrays', default=''); ap.add_argument('--seed', type=int, default=61)
    ap.add_argument('--gaps', default='', help='assembly-gap BED: chains never cross a gap (v4.7)')
    ap.add_argument('--null-cross', type=int, default=100, help='flanks from other chains per (chain, type, flank) null (D06)')
    ap.add_argument('--null-shift', type=int, default=5, help='shifted spacer flanks per junction in the null (D06)')
    ap.add_argument('--sweep', action='store_true'); ap.add_argument('--sweep-out', default='')
    a = ap.parse_args(); random.seed(a.seed); os.makedirs(a.dest, exist_ok=True)
    FL = [int(x) for x in a.flank.split(',')]; CU = [float(x) for x in a.cutoff.split(',')]; TO = [int(x) for x in a.tol.split(',')]
    keep = set(a.arrays.split(',')) if a.arrays else None
    fa = Fasta(a.genome)
    rows = list(csv.DictReader(open(a.catalogue), delimiter='\t'))
    GAPS = Gaps(a.gaps) if a.gaps else None
    chains, JX, st = build_chains(rows, max_gap=a.max_gap, min_runs=a.min_runs, gaps=GAPS, keep=keep)
    need = {c['array_id'] for c in chains}
    print(f'{len(chains)} chains ({st["junctions"]} junctions); links refused at assembly gaps: {st["links_refused_at_gaps"] if GAPS is not None else "no gap BED given"}, '
          f'at a family change: {st["links_refused_family"]}', flush=True)
    # dominant class consensus and copy orientation per array (03 tables are keyed on the v46 array id, D12)
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
    jx_by = collections.defaultdict(list)
    for j in JX: jx_by[j['chain_id']].append(j)
    J = []; units = []; skipped = []; info = {}
    # pass 1: refine every junction edge and build its context
    for ci, ch in enumerate(chains):
        aid, rs, chain_id, chrom = ch['array_id'], ch['runs'], ch['chain_id'], ch['chrom']
        top = cls[aid].most_common(1)
        if not top or (aid, top[0][0]) not in cons:
            skipped.append(dict(chain=chain_id, array_id=aid, chrom=chrom, runs=len(rs), reason='no_copy_class' if not top else 'no_consensus_for_class'))
            continue
        c0 = cons[(aid, top[0][0])][1]; is_rev = rev[aid]['yes'] > rev[aid]['no']
        cg = rc(c0) if is_rev else c0          # consensus in genomic orientation of the copies
        info[chain_id] = (is_rev, chrom, cg)
        for fl in FL:
            for jx in jx_by[chain_id]:
                r1, r2, i = jx['left'], jx['right'], jx['gap_index']
                for side, edge, run in (('R', jx['left_end'], r1['run_id']), ('L', jx['right_start0'], r2['run_id'])):
                    ref = refine(fa, chrom, edge, side, cg, fl)
                    jt = ('OUT' if side == 'R' else 'IN') if not is_rev else ('IN' if side == 'R' else 'OUT')
                    base = dict(chain=chain_id, array_id=aid, family=ch['family'], chrom=chrom, run_id=run, gap_index=i, type=jt, flank=fl,
                                hmm_edge=edge, spacer_bp=jx['spacer_bp'])
                    if ref is None: J.append(dict(base, refined=None)); continue
                    pos, ph = ref; s, ei = context(fa, chrom, pos, side, is_rev, a.inner, fl)
                    # array orientation: OUT = satellite then flank (edge index = inner), IN = flank then satellite
                    sat, flk = (s[:ei], s[ei:]) if jt == 'OUT' else (s[ei:], s[:ei])
                    rec = dict(base, refined=pos, delta_hmm=pos - edge, phase=ph, cons_len=len(cg), _s=s, _e=ei, _side=side, _sat=sat, _fl=flk)
                    # shifted-spacer null flank (D06): a spacer piece of the same length taken fl + inner further into the spacer; the null
                    # context keeps the junction's OWN satellite side (anchored alignments of random flanks score ~0.5 identity, unanchored ~0)
                    sp = jx['spacer_bp']; room = sp - (fl + a.inner) - fl; rec['_nullfl'] = []
                    for _ in range(a.null_shift if room > 0 else 0):
                        d = random.randint(fl + a.inner, sp - fl); se = pos + d if side == 'R' else pos - d
                        ns, ne = context(fa, chrom, se, side, is_rev, a.inner, fl)
                        rec['_nullfl'].append(ns[ne:] if jt == 'OUT' else ns[:ne])
                    J.append(rec)
        if ci % 50 == 0: print(f'{ci}/{len(chains)} arrays', flush=True)
    # pass 2: homology per (chain, type, flank) against the medoid, with the per-array null (shifted contexts + other chains' contexts)
    pool = collections.defaultdict(list)
    for x in J:
        if x.get('_s') is not None: pool[(x['type'], x['flank'])].append(x)
    by = collections.defaultdict(list)
    for x in J:
        if x.get('_s') is not None: by[(x['chain'], x['type'], x['flank'])].append(x)
    for (chain_id, jt, fl), L_ in by.items():
        if len(L_) < 3: continue
        smp = random.sample(L_, min(40, len(L_)))
        med = max(smp, key=lambda x: np.mean([LOC.score(x['_s'], y['_s']) for y in smp if y is not x]))
        t = med['_s']; f_lo, f_hi = (med['_e'], len(t)) if jt == 'OUT' else (0, med['_e'])
        for x in L_:
            al = GLO.align(x['_s'], t)[0]; mp = map_edge(al, x['_e'])
            x['offset'] = None if mp is None else mp - med['_e']
            x['flank_identity'], x['flank_aligned'] = flank_identity(al, x['_s'], t, f_lo, f_hi)
            x['medoid_run'] = med['run_id']; x['is_medoid'] = int(x is med)
        def compose(x, flank): return x['_sat'] + flank if jt == 'OUT' else flank + x['_sat']
        nq = [compose(x, f_) for x in L_ for f_ in x['_nullfl']]
        others = [y for y in pool[(jt, fl)] if y['chain'] != chain_id]
        for y in random.sample(others, min(a.null_cross, len(others))): nq.append(compose(random.choice(L_), y['_fl']))
        nv = [flank_identity(GLO.align(q, t)[0], q, t, f_lo, f_hi)[0] for q in nq]
        p99 = float(np.percentile(nv, 99)) if len(nv) >= 10 else None
        for x in L_:                               # the medoid is identical to itself by construction: not testable
            x['null_n'] = len(nv); x['null_p99'] = None if p99 is None else round(p99, 4)
            x['homologous_null'] = None if (p99 is None or x['is_medoid']) else int(x['flank_identity'] > p99 and x['flank_identity'] >= min(CU))
    # refined units (main flank = first in list): IN junction of run i -> IN junction of run i+1, array orientation (D11: satellite / interruption split)
    for chain_id, (is_rev, chrom, cg) in info.items():
        mine = [x for x in J if x['chain'] == chain_id and x['flank'] == FL[0] and x.get('refined')]
        ins = sorted([x for x in mine if x['type'] == 'IN'], key=lambda x: x['refined'])
        outs = {x['run_id']: x['refined'] for x in mine if x['type'] == 'OUT'}
        for x, y in zip(ins, ins[1:]):               # the unit's satellite run: the left one on '+', the right one on '-' (array orientation)
            u0, u1 = min(x['refined'], y['refined']), max(x['refined'], y['refined']); sx = y if is_rev else x
            sat = abs(outs[sx['run_id']] - sx['refined']) if sx['run_id'] in outs else ''
            units.append((chrom, u0, u1, f"{chain_id}|{sx['run_id']}", '-' if is_rev else '+', sat, (u1 - u0 - sat) if sat != '' else ''))
    cols = ['chain', 'array_id', 'family', 'chrom', 'run_id', 'gap_index', 'type', 'flank', 'hmm_edge', 'spacer_bp', 'refined', 'delta_hmm', 'phase', 'cons_len',
            'flank_identity', 'flank_aligned', 'offset', 'medoid_run', 'is_medoid', 'null_n', 'null_p99', 'homologous_null']
    with open(os.path.join(a.dest, f'61_junctions_{a.tag}.tsv'), 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t'); w.writerow(cols)
        for x in J: w.writerow(['' if x.get(c) is None else x.get(c, '') for c in cols])
    with open(os.path.join(a.dest, f'61_skipped_chains_{a.tag}.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['chain', 'array_id', 'chrom', 'runs', 'reason'], delimiter='\t'); w.writeheader(); w.writerows(skipped)
    for s_ in skipped: print('skipped chain', s_['chain'], s_['reason'])
    # array scores
    S = []
    for (ch, jt, fl), L_ in by.items():
        ok = [x for x in L_ if x.get('refined') is not None and x.get('flank_identity') is not None and not x.get('is_medoid')]
        if len(ok) < 2: continue
        ph = np.array([x['phase'] for x in ok]); L = ok[0]['cons_len']; mode = collections.Counter(ph).most_common(1)[0][0]
        circ = np.minimum((ph - mode) % L, (mode - ph) % L)
        hn = [x for x in ok if x.get('homologous_null')]
        for cu in CU:
            hom = [x for x in ok if x['flank_identity'] >= cu]
            row = dict(chain=ch, array_id=ok[0]['array_id'], family=ok[0]['family'], chrom=ok[0]['chrom'], type=jt, flank=fl, cutoff=cu,
                       junctions=len(L_), refined=len(ok), homologous_share=len(hom) / len(ok), homologous_share_null=len(hn) / len(ok),
                       null_p99=ok[0].get('null_p99', ''), phase_concordance=float(np.mean(circ <= 3)),
                       median_abs_delta_hmm=float(np.median([abs(x['delta_hmm']) for x in ok])))
            for tol in TO:
                row[f'precise_pm{tol}'] = float(np.mean([abs(x['offset']) <= tol for x in hom if x['offset'] is not None])) if hom else np.nan
            S.append(row)
    if S:
        with open(os.path.join(a.dest, f'61_array_scores_{a.tag}.tsv'), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(S[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(S)
    with open(os.path.join(a.dest, f'61_units_{a.tag}.bed'), 'w') as fh:
        for u in units: fh.write(f'{u[0]}\t{u[1]}\t{u[2]}\t{u[3]}\t0\t{u[4]}\t{u[5]}\t{u[6]}\n')
    # summary: per setting, median over arrays
    agg = collections.defaultdict(list)
    for r in S: agg[(r['type'], r['flank'], r['cutoff'])].append(r)
    SM = []
    for (jt, fl, cu), rr in sorted(agg.items()):
        SM.append(dict(type=jt, flank=fl, cutoff=cu, arrays=len(rr), median_homologous_share=round(float(np.median([r['homologous_share'] for r in rr])), 3),
                       median_homologous_share_null=round(float(np.median([r['homologous_share_null'] for r in rr])), 3),
                       **{f'median_precise_pm{t}': round(float(np.nanmedian([r[f'precise_pm{t}'] for r in rr])), 3) for t in TO},
                       median_phase_concordance=round(float(np.median([r['phase_concordance'] for r in rr])), 3)))
    def write(path, rows_):
        with open(path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows_[0].keys()) if rows_ else ['type'], delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows_)
    write(os.path.join(a.dest, f'61_summary_{a.tag}.tsv'), SM)
    if a.sweep:
        out = a.sweep_out or os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'benchmarks', '61_sweep.tsv')
        os.makedirs(os.path.dirname(out), exist_ok=True); write(out, [dict(step='61', tag=a.tag, **r) for r in SM]); print('sweep ->', out)
    print(f'chains {len(chains)}, skipped {len(skipped)}, junctions {len(J)}, units {len(units)}')


if __name__ == '__main__':
    main()
