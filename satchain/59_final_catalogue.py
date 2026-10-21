#!/usr/bin/env python3
"""
59_final_catalogue.py  (v4.6; the validated satellite catalogue of one subgenome)

Runs = (a) copy-chain HMM runs of the existing arrays that passed the per-run sequence check (56, p90 set) and (b) HMM runs of the
new scan-57 arrays (decomposed by step 03; their candidate already passed the p99.9 excess and the 3-copy chain / TRF test).
Per run, from its step-03 copies:
  n_conf (complete, identity >= 80), n_weak (65-80), n_junk (< 65); mean identity of confident copies
  G1 dominant-class consensus GC within 5 points of the run's own GC; G2 mean identity >= 75 % (unit > 60 bp) or >= 85 % (<= 60 bp);
  G3 strand coherent (reverse fraction outside 0.35-0.65)  -- the step-07 rules, applied to the run instead of the scan interval
  family = majority family of confident copies (existing arrays) or the family matched by the new array's class consensus (BLAST >= 80 %
  identity over >= 50 % of the consensus; 'new:<array>' otherwise); lineage from 55d (A / B / cen_minisatellite; 'none' otherwise,
  'low_n' when < 5 copies were sampled)
  unit class (<= 6 micro, 7-60 mini, > 60 satellite), EDTA overlap fraction (NOT TE embedding: EDTA labels satellites as TE), overlap with a putative centromeric array (31; L1)
Tiers by cumulative share of the catalogue's bp, families ranked by bp: core (first 50 %), major (to 90 %), minor (to 99 %), rare.
Every row: evidence 'L1' (sequence-validated tandem array; no functional claim).
Outputs: 59_catalogue_<G>.tsv, 59_catalogue_<G>.bed (BED9 by lineage), 59_family_tiers_<G>.tsv, 59_summary_<G>.tsv
usage: 59_final_catalogue.py <G> <pipeline out> <hmm56 dir> <new03 dir> <hmm56new dir> <famnew tsv> <lin55d tsv> <func folder> <genome fasta> <dest>
"""
import sys, os, csv, collections, bisect, gzip
csv.field_size_limit(1 << 30)
G, OUT, H, N03, HN, FAMNEW, LIN, FUNC, FA, DEST = sys.argv[1:11]; os.makedirs(DEST, exist_ok=True)
P = os.path.join(OUT, f'fClaHyb_{G}_satellites')
lin = {}
for r in csv.DictReader(open(LIN), delimiter='\t'):
    if r['subgenome'] == G: lin[r['family']] = r['lineage'] if int(r['n_sampled_copies']) >= 5 else r['lineage'] + '_low_n'
famnew = {}
if os.path.exists(FAMNEW):
    for r in csv.DictReader(open(FAMNEW), delimiter='\t'): famnew[(r['array_id'], r['monomer_class'])] = r['family']


def load_copies(path, keep):
    cop = collections.defaultdict(list)
    with open(path, newline='') as fh:
        rd = csv.reader(fh, delimiter='\t'); h = next(rd); ix = {k: i for i, k in enumerate(h)}
        for row in rd:
            if row[ix['array_id']] not in keep: continue
            try: idt = float(row[ix['best_alignment_identity']])
            except ValueError: idt = 0.0
            cop[row[ix['array_id']]].append((int(row[ix['genomic_start_0based']]), int(row[ix['genomic_end']]), idt,
                                             row[ix['complete_or_truncated']] == 'complete', row[ix['strand']], row[ix['family_id']], row[ix['monomer_class']]))
    for v in cop.values(): v.sort()
    return cop


def runs(bed):
    out = []
    for l in open(bed):
        c, s, e, nm = l.split('\t')[:4]; out.append((c, int(s), int(e), nm.split('|')[0]))
    return out


R_old = runs(os.path.join(H, '56_runs_validated_q0.9.bed')); R_new = runs(os.path.join(HN, '56_runs_post0.95.bed')) if os.path.exists(os.path.join(HN, '56_runs_post0.95.bed')) else []
C_old = load_copies(os.path.join(P, '03_monomer_copies.tsv'), {r[3] for r in R_old})
C_new = load_copies(os.path.join(N03, '03_monomer_copies.tsv'), {r[3] for r in R_new}) if R_new else {}
per = {}
for path in (os.path.join(P, '02_satellite_arrays.tsv'), os.path.join(N03, '02_satellite_arrays.tsv')):
    if os.path.exists(path):
        for r in csv.DictReader(open(path), delimiter='\t'): per[r['array_id']] = float(r['period_mode_bp'] or 0)
cgc = {}
for path in (os.path.join(P, '03_array_consensus.tsv'), os.path.join(N03, '03_array_consensus.tsv')):
    if os.path.exists(path):
        for r in csv.DictReader(open(path), delimiter='\t'): cgc[(r['array_id'], r['monomer_class'])] = float(r['consensus_GC_percent'])
te = collections.defaultdict(list)
for l in gzip.open(os.path.join(FUNC, '_repeats_EDTA', f'fClaHyb_{G}.TEanno.sorted.gff3.gz'), 'rt'):
    if l.startswith('#'): continue
    f = l.split('\t', 5)
    if len(f) > 4 and f[2] not in ('repeat_region', 'target_site_duplication', 'long_terminal_repeat'): te[f[0]].append((int(f[3]) - 1, int(f[4])))
for v in te.values(): v.sort()
cen = collections.defaultdict(list)
for r in csv.DictReader(open(os.path.join(P, '31_centromere_summary.tsv')), delimiter='\t'): cen[r['chrom']].append((int(r['cen_start_0based']), int(r['cen_end'])))
def ov(d, c, s, e):
    v = d.get(c, []); i = bisect.bisect_left(v, (s - 2000000, 0)); t = 0
    for x in v[i:]:
        if x[0] >= e: break
        t += max(0, min(e, x[1]) - max(s, x[0]))
    return t
fai = {l.split()[0]: list(map(int, l.split()[1:5])) for l in open(FA + '.fai')}; fh_ = open(FA, 'rb')
def gc_of(c, s, e):
    L, off, lb, lB = fai[c]; fh_.seek(off + (s // lb) * lB + s % lb); q = fh_.read(int((e - s) * lB / lb) + 2 * lB).replace(b'\n', b'')[:e - s].upper()
    n = len(q) - q.count(b'N'); return 100 * (q.count(b'G') + q.count(b'C')) / max(1, n)
rows = []
for src, R, C in (('existing', R_old, C_old), ('new', R_new, C_new)):
    for k, (c, s, e, aid) in enumerate(R):
        cs = [x for x in C.get(aid, []) if x[0] >= s and x[1] <= e]
        conf = [x for x in cs if x[3] and x[2] >= 80]
        if not conf: continue
        u = per.get(aid, 0); ucls = 'micro' if u <= 6 else 'mini' if u <= 60 else 'satellite'
        mid = sum(x[2] for x in conf) / len(conf); rev = sum(x[4] == '-' for x in conf) / len(conf)
        domcls = collections.Counter(x[6] for x in conf).most_common(1)[0][0]
        if src == 'existing': fam = collections.Counter(x[5] for x in conf).most_common(1)[0][0]
        else: fam = famnew.get((aid, domcls), f'new:{aid}')
        rgc = gc_of(c, s, e); cg = cgc.get((aid, domcls))
        g1 = 'PASS' if cg is not None and abs(rgc - cg) <= 5 else 'FAIL'
        g2 = 'PASS' if mid >= (85 if u <= 60 else 75) else 'FAIL'
        g3 = 'PASS' if not (0.35 < rev < 0.65) else 'FAIL'
        rows.append(dict(run_id=f'{G}_{"E" if src == "existing" else "N"}{k + 1:05d}', subgenome=G, chrom=c, start0=s, end=e, bp=e - s, source=src, array_id=aid,
                         unit_bp=int(u), unit_class=ucls, n_conf=len(conf), n_weak=sum(1 for x in cs if 65 <= x[2] < 80 or (not x[3] and x[2] >= 80)),
                         n_junk=sum(1 for x in cs if x[2] < 65), mean_identity_conf=round(mid, 2), reverse_fraction=round(rev, 3), run_GC=round(rgc, 2),
                         consensus_GC=cg if cg is not None else 'NA', G1=g1, G2=g2, G3=g3, gates_passed=(g1, g2, g3).count('PASS'),
                         family=fam, lineage=lin.get(fam, 'none'), EDTA_overlap_fraction=round(ov(te, c, s, e) / (e - s), 3),   # EDTA often labels the satellite itself as TE: not 'embedding'
                         putative_centromeric_array_L1='yes' if ov(cen, c, s, e) > 0 else 'no', evidence='L1'))
# tiers by bp share
fbp = collections.Counter()
for r in rows: fbp[r['family']] += r['bp']
tot = sum(fbp.values()); cum = 0; tier = {}
for f, b in fbp.most_common():
    cum += b; frac_before = (cum - b) / tot
    tier[f] = 'core' if frac_before < 0.5 else 'major' if frac_before < 0.9 else 'minor' if frac_before < 0.99 else 'rare'
for r in rows: r['tier'] = tier[r['family']]
cols = list(rows[0].keys())
with open(os.path.join(DEST, f'59_catalogue_{G}.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
COL = {'A': '31,119,180', 'B': '214,39,40', 'cen_minisatellite': '148,103,189'}
with open(os.path.join(DEST, f'59_catalogue_{G}.bed'), 'w') as fh:
    for r in rows:
        fh.write(f"{r['chrom']}\t{r['start0']}\t{r['end']}\t{r['run_id']}|{r['family']}|{r['lineage']}|{r['unit_bp']}bp|{r['tier']}\t0\t.\t{r['start0']}\t{r['end']}\t{COL.get(r['lineage'].replace('_low_n', ''), '150,150,150')}\n")
with open(os.path.join(DEST, f'59_family_tiers_{G}.tsv'), 'w') as fh:
    fh.write('family\tbp\tshare\ttier\tlineage\n')
    for f, b in fbp.most_common(): fh.write(f'{f}\t{b}\t{b / tot:.4f}\t{tier[f]}\t{lin.get(f, "none")}\n')
S = collections.Counter(); B = collections.Counter()
for r in rows:
    for k in (('source', r['source']), ('unit_class', r['unit_class']), ('lineage', r['lineage']), ('tier', r['tier']),
              ('gates_passed', r['gates_passed']), ('putative_centromeric_array_L1', r['putative_centromeric_array_L1']),
              ):
        S[k] += 1; B[k] += r['bp']
with open(os.path.join(DEST, f'59_summary_{G}.tsv'), 'w') as fh:
    fh.write('kind\tvalue\tn_runs\tMb\n')
    for k in sorted(S, key=str): fh.write(f'{k[0]}\t{k[1]}\t{S[k]}\t{B[k] / 1e6:.2f}\n')
print(G, len(rows), 'runs,', round(tot / 1e6, 2), 'Mb'); print(open(os.path.join(DEST, f'59_summary_{G}.tsv')).read())
