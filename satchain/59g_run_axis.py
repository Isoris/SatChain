#!/usr/bin/env python3
"""59g_run_axis.py  (SatChain v4.7, 2026-11-09; the 59g sequence axis per RUN, for the 5S-related families)

Purpose (measurement-purpose rule): evidence for the hand decision on MacSat06 (59g: unresolved; the unit consensus carries an intact 5S
gene, the 59b ratio 0.726 is within 0.1 of the cut). 59g reads ONE unit per family (the dominant-class consensus of the largest run), so
it cannot say whether every run of a family carries the gene. This step asks that per run. It adds NO rule and changes no 59g label.
Controls: GarSat49 (5S gene + spacer, gene_array in 59g) and GarSat04 (5S-derived satellite, no intact gene).

Families: every family whose 59g best gene hit is 5S or whose 59b class mentions 5S (--families overrides), runs in both haplotypes.
Per run:
  a. run consensus = step-03 consensus of the run's (array_id without the 59 '~g' suffix, dominant_class), doubled (circular), both strands;
  b. up to --copies confident copies of the run (step-03 identity >= 80, complete; copies inside the run coordinates; seeded sample),
     each doubled (a copy can start inside the gene);
  c. every query: 13-mer screen against the 59g 5S references (59g_refs_<G>.fa, records '5S|...', both strands; >= 3 shared k-mers,
     <= 6 candidates), parasail Smith-Waterman with the 59g scoring (match 2 / mismatch -3 / gap open 5, extend 2); best hit =
     identity and reference coverage.
  d. axis with the 59g cut-offs: intact_gene (identity >= 0.90, coverage >= 0.80) / gene_derived (>= 0.70, >= 0.50) / none.
     organisation of the run consensus: gene_plus_spacer when its length <= 2.5 x the gene, else gene_in_long_unit.
Per family: runs and bp by run-consensus axis, share of sampled copies per axis (pooled), median copy identity.
Sweep (for the record): identity 0.80 / 0.90 x coverage 0.6 / 0.8 on the run-consensus axis (bp share intact).
Outputs (<dest>/): run_axis_runs.tsv, run_axis_copies.tsv, run_axis_family.tsv, run_axis_sweep.tsv, log.txt
usage: 59g_run_axis.py <final dir> <v4.5 out root> <new03 dir pattern with {G}> <dest> [--families F1,F2] [--copies 30] [--seed 59]
Evidence level L1 (sequence)."""
import os, sys, csv, bisect, random, argparse, statistics, collections
import parasail
csv.field_size_limit(1 << 30)
ap = argparse.ArgumentParser(); ap.add_argument('final'); ap.add_argument('out45'); ap.add_argument('new03'); ap.add_argument('dest')
ap.add_argument('--families', default=''); ap.add_argument('--copies', type=int, default=30); ap.add_argument('--seed', type=int, default=59)
ap.add_argument('--k', type=int, default=13); ap.add_argument('--min-kmers', type=int, default=3); ap.add_argument('--max-cand', type=int, default=6)
a = ap.parse_args(); os.makedirs(a.dest, exist_ok=True)
LOG = open(os.path.join(a.dest, 'log.txt'), 'w', encoding='utf-8')
def P(*x): s = ' '.join(str(v) for v in x); print(s, flush=True); LOG.write(s + '\n')
MAT = parasail.matrix_create('ACGT', 2, -3); RC = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(RC)[::-1]
def tsv(p): return list(csv.DictReader(open(p, encoding='utf-8'), delimiter='\t'))
def read_fa(p):
    d, n, s = [], None, []
    for l in open(p, encoding='utf-8'):
        l = l.strip()
        if l.startswith('>'):
            if n: d.append((n, ''.join(s).upper()))
            n, s = l[1:].split()[0], []
        elif l: s.append(l)
    if n: d.append((n, ''.join(s).upper()))
    return d

class Refs:
    def __init__(self, recs, k):
        seen = {}; self.k = k
        for n, s in recs:
            if s not in seen: seen[s] = n
        self.seqs = [(n, s) for s, n in seen.items()]; self.idx = collections.defaultdict(set)
        for i, (n, s) in enumerate(self.seqs):
            for j in range(len(s) - k + 1): self.idx[s[j:j + k]].add(i)
        self.gene_len = statistics.median(len(s) for _, s in self.seqs)
    def best(self, q):
        """q = doubled unit / copy. Same definitions as 59g sw(): identity = '|' columns / aligned columns, reference coverage =
        aligned reference span / reference length, hits < 20 columns ignored, best = highest score; query vs reference on both strands."""
        best = None
        for strand, qq in (('+', q), ('-', rc(q))):   # query reverse complement == reference reverse complement for a local alignment
            cnt = collections.Counter()
            for km in {qq[j:j + self.k] for j in range(len(qq) - self.k + 1)}:
                for i in self.idx.get(km, ()): cnt[i] += 1
            for i, c in cnt.most_common(a.max_cand):
                if c < a.min_kmers: break
                ref = self.seqs[i][1]; r = parasail.sw_trace_striped_16(qq, ref, 5, 2, MAT); tb = r.traceback
                cols = len(tb.query)
                if cols < 20: continue
                ident = sum(1 for x in tb.comp if x == '|') / cols; tl = len(tb.ref.replace('-', '')); rcov = tl / len(ref)
                h = (int(r.score), round(ident, 3), round(rcov, 3), self.seqs[i][0])
                if best is None or h[0] > best[0]: best = h
        return best   # (score, identity, ref_cov, ref) or None

def axis(ident, cov, idc=0.90, covc=0.80):
    if ident >= idc and cov >= covc: return 'intact_gene'
    if ident >= 0.70 and cov >= 0.50: return 'gene_derived'
    return 'none'

fams = set(x for x in a.families.split(',') if x)
if not fams:
    for G in ('Gar', 'Mac'):
        for r in tsv(os.path.join(a.final, 'unit_content_2026-11-07', f'59g_unit_content_{G}.tsv')):
            if r['gene_best_class'] == '5S' or '5S' in r['element_class_59b']: fams.add(r['family'])
P('families', sorted(fams))
rng = random.Random(a.seed); RUNS, COPIES = [], []
for G in ('Gar', 'Mac'):
    refs = Refs([(n, s) for n, s in read_fa(os.path.join(a.final, 'unit_content_2026-11-07', f'59g_refs_{G}.fa')) if n.startswith('5S|')], a.k)
    P(G, 'unique 5S references', len(refs.seqs), 'median length', refs.gene_len)
    cons = {}
    for p in (os.path.join(a.out45, f'fClaHyb_{G}_satellites', '03_array_consensus.tsv'), os.path.join(a.new03.format(G=G), '03_array_consensus.tsv')):
        for r in tsv(p): cons[(r['array_id'], r['monomer_class'])] = r['consensus_sequence'].upper()
    runs = [r for r in tsv(os.path.join(a.final, 'cat59', f'59_catalogue_{G}.tsv')) if r['family_v46'] in fams]
    P(G, 'runs of the selected families', len(runs))
    by_chrom = collections.defaultdict(list)
    for r in runs: by_chrom[r['chrom']].append((int(r['start0']), int(r['end']), r['run_id']))
    for c in by_chrom: by_chrom[c].sort()
    starts = {c: [x[0] for x in v] for c, v in by_chrom.items()}
    copies = collections.defaultdict(list)
    for p in (os.path.join(a.out45, f'fClaHyb_{G}_satellites', '03_monomer_copies.tsv'), os.path.join(a.new03.format(G=G), '03_monomer_copies.tsv')):
        with open(p, encoding='utf-8', newline='') as fh:
            h = fh.readline().rstrip('\r\n').split('\t'); I = {k: h.index(k) for k in ('chrom', 'genomic_start_0based', 'genomic_end', 'best_alignment_identity', 'complete_or_truncated', 'sequence')}
            for line in fh:
                x = line.rstrip('\r\n').split('\t'); c = x[I['chrom']]
                if c not in by_chrom: continue
                s0, e0 = int(x[I['genomic_start_0based']]), int(x[I['genomic_end']]); j = bisect.bisect_right(starts[c], s0) - 1
                if j < 0: continue
                rs, re_, rid = by_chrom[c][j]
                if s0 >= rs and e0 <= re_ and float(x[I['best_alignment_identity']] or 0) >= 80 and x[I['complete_or_truncated']] == 'complete':
                    copies[rid].append(x[I['sequence']].upper())
    for r in runs:
        aid = r['array_id'].split('~')[0]; u = cons.get((aid, r['dominant_class']), '')
        row = dict(haplotype=G, run_id=r['run_id'], family=r['family_v46'], chrom=r['chrom'], start0=r['start0'], end=r['end'], bp=r['bp'],
                   element_class=r['element_class'], unit_bp_catalogue=r['unit_bp'], consensus_len=len(u))
        b = refs.best(u + u) if u else None
        row.update(cons_ident=b[1] if b else 0.0, cons_refcov=b[2] if b else 0.0, cons_ref=b[3] if b else '',
                   cons_axis=axis(b[1], b[2]) if b else ('no_consensus' if not u else 'none'),
                   cons_organization=('gene_plus_spacer' if len(u) <= 2.5 * refs.gene_len else 'gene_in_long_unit') if (b and axis(b[1], b[2]) != 'none') else 'no_gene')
        cs = copies.get(r['run_id'], []); sample = rng.sample(cs, min(a.copies, len(cs)))
        ax = collections.Counter(); ids = []
        for i, s in enumerate(sample):
            bb = refs.best(s + s); x_ = axis(bb[1], bb[2]) if bb else 'none'; ax[x_] += 1; ids.append(bb[1] if bb else 0.0)
            COPIES.append(dict(haplotype=G, run_id=r['run_id'], family=r['family_v46'], copy_index=i, copy_len=len(s), ident=bb[1] if bb else 0.0, refcov=bb[2] if bb else 0.0, axis=x_))
        n = len(sample)
        row.update(copies_confident=len(cs), copies_sampled=n, share_intact=round(ax['intact_gene'] / n, 3) if n else '', share_derived=round(ax['gene_derived'] / n, 3) if n else '',
                   share_none=round(ax['none'] / n, 3) if n else '', median_copy_ident=round(statistics.median(ids), 3) if ids else '')
        RUNS.append(row)
    P(G, 'runs done', sum(1 for r in RUNS if r['haplotype'] == G))

def write(name, rows):
    with open(os.path.join(a.dest, name), 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
write('run_axis_runs.tsv', RUNS); write('run_axis_copies.tsv', COPIES)
FAM = []
for G, fam in sorted({(r['haplotype'], r['family']) for r in RUNS}):
    rs = [r for r in RUNS if r['haplotype'] == G and r['family'] == fam]
    tot = sum(int(r['bp']) for r in rs); cps = [c for c in COPIES if c['haplotype'] == G and c['family'] == fam]
    d = dict(haplotype=G, family=fam, runs=len(rs), bp=tot)
    for k in ('intact_gene', 'gene_derived', 'none', 'no_consensus'):
        d[f'runs_{k}'] = sum(1 for r in rs if r['cons_axis'] == k); d[f'bp_share_{k}'] = round(sum(int(r['bp']) for r in rs if r['cons_axis'] == k) / tot, 3) if tot else ''
    d['runs_gene_plus_spacer'] = sum(1 for r in rs if r['cons_organization'] == 'gene_plus_spacer')
    d['copies_sampled'] = len(cps)
    for k in ('intact_gene', 'gene_derived', 'none'): d[f'copy_share_{k}'] = round(sum(1 for c in cps if c['axis'] == k) / len(cps), 3) if cps else ''
    d['median_copy_ident'] = round(statistics.median(c['ident'] for c in cps), 3) if cps else ''
    FAM.append(d)
write('run_axis_family.tsv', FAM)
SW = []
for idc in (0.80, 0.90):
    for covc in (0.6, 0.8):
        for G, fam in sorted({(r['haplotype'], r['family']) for r in RUNS}):
            rs = [r for r in RUNS if r['haplotype'] == G and r['family'] == fam]; tot = sum(int(r['bp']) for r in rs)
            SW.append(dict(identity_cut=idc, coverage_cut=covc, haplotype=G, family=fam,
                           bp_share_intact=round(sum(int(r['bp']) for r in rs if axis(r['cons_ident'], r['cons_refcov'], idc, covc) == 'intact_gene') / tot, 3) if tot else ''))
write('run_axis_sweep.tsv', SW)
for d in FAM: P(d)
P('DONE')
