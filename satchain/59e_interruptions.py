#!/usr/bin/env python3
"""59e: interruptions between neighbouring catalogue runs of the same family (gap <= --max-gap), classified without EDTA's satellite labels.
For each gap: EDTA features whose label is intrinsic to the flanking runs (59c intrinsic_labels at --intrinsic-cut) are ignored; the rest give
TE coverage. Class: TE_insertion (TE >= --class-cut of the gap), missed_TE_54 (a step-54 known-TE insertion family covers >= cut),
other:exon_CDS / other:intron / other:ncRNA / other:SSR (the covering annotation type; 2026-11-01, audit C16: 'other_annotation' used to
pool whole gene spans with exons and never loaded SSR), unannotated.
v2 (2026-10-21): the intrinsic-label filter is per feature, not per label. An EDTA feature with a label intrinsic to the flanking runs still
counts as TE when it is a whole element sitting in the gap (>= OWN of its length inside the gap and >= MINLEN bp): EDTA gives MacSat03 and a
7.5-kb LTR element interleaved with it the same label (LTR/unknown), and the v1 filter dropped the element.
v4.7 (2026-10-26): no interruption across an assembly gap (env SATCHAIN_GAPS or --gaps = gap BED); reported as refused; 2026-11-01: a
zero-length contig join exactly on a run edge is refused too (closed test, audit HB07).
2026-11-01 (audit C16): neighbours are consecutive SATELLITE runs of the same family (family_v46 when present) on a chromosome, whatever
their array_id (the old rule chained inside one v4.5 array_id only, so runs of one family in adjacent candidate intervals were never
examined and runs of different families inside one interval counted as one array). --sweep writes 59e_sweep_<G>.tsv over class-cut
0.3 / 0.5 / 0.7 x intrinsic-cut 0.2 / 0.3 / 0.5 (the 59c table carries all three) x max-gap 10 / 20 / 50 kb.
Output: <cat59 dir>/59e_interruptions_<G><suffix>.tsv and a summary.
usage: 59e_interruptions.py <G> <cat59 dir> <func folder> <ins54 dir> [--max-gap 50000] [--own 0.8] [--min-len 500] [--class-cut 0.5]
       [--intrinsic-cut 0.3] [--gaps <bed>] [--suffix s] [--sweep] [--all-classes]"""
import sys, csv, gzip, collections, bisect, os, glob, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from satpaths import func_path, func_dir; from gaps import Gaps
ap = argparse.ArgumentParser(); ap.add_argument('G'); ap.add_argument('C'); ap.add_argument('FUNC'); ap.add_argument('INS')
ap.add_argument('--max-gap', type=int, default=50000); ap.add_argument('--own', type=float, default=0.8); ap.add_argument('--min-len', type=int, default=500)
ap.add_argument('--class-cut', type=float, default=0.5); ap.add_argument('--intrinsic-cut', default='0.3'); ap.add_argument('--gaps', default=os.environ.get('SATCHAIN_GAPS', ''))
ap.add_argument('--suffix', default=''); ap.add_argument('--sweep', action='store_true'); ap.add_argument('--all-classes', action='store_true', help='chain gene arrays too (default: satellite runs only)')
a = ap.parse_args(); G, C, FUNC, INS = a.G, a.C, a.FUNC, a.INS; MAXGAP, OWN, MINLEN, SUF = a.max_gap, a.own, a.min_len, a.suffix
runs = [r for r in csv.DictReader(open(f'{C}/59_catalogue_{G}.tsv'), delimiter='\t') if a.all_classes or r['element_class'] == 'satellite']
FAM = lambda r: r.get('family_v46') or r['family']
intr = collections.defaultdict(dict)
for r in csv.DictReader(open(f'{C}/59c_te_{G}.tsv'), delimiter='\t'):
    for k, v in r.items():
        if k.startswith('intrinsic_labels_'): intr[k.split('_')[-1]][r['run_id']] = set(x for x in v.split(',') if x != 'none')
GAPS = Gaps(a.gaps) if a.gaps else None
te = collections.defaultdict(list)
for l in gzip.open(func_path(FUNC, G, '_repeats_EDTA', f'fClaHyb_{G}.TEanno.sorted.gff3.gz'), 'rt'):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] in ('repeat_region', 'target_site_duplication', 'long_terminal_repeat'): continue
    te[f[0]].append((int(f[3]) - 1, int(f[4]), dict(x.split('=', 1) for x in f[8].split(';') if '=' in x).get('classification', f[2])))
for v in te.values(): v.sort()
other = {k: collections.defaultdict(list) for k in ('exon_CDS', 'gene', 'ncRNA', 'SSR')}
for l in gzip.open(func_path(FUNC, G, '_genes_lifton', f'fClaHyb_{G}.lifton.sorted.gff3.gz'), 'rt'):
    if l.startswith('#'): continue
    f = l.split('\t')
    if len(f) > 4 and f[2] in ('exon', 'CDS'): other['exon_CDS'][f[0]].append((int(f[3]) - 1, int(f[4])))
    elif len(f) > 4 and f[2] == 'gene': other['gene'][f[0]].append((int(f[3]) - 1, int(f[4])))
for sub, name in (('_ncRNA_rfam', 'rfam.all_classes'), ('_tRNA_tRNAscan', 'tRNA'), ('_rRNA_barrnap', 'rRNA'), ('_ncRNA_all', 'ncRNA_all')):
    p = func_path(FUNC, G, sub, f'fClaHyb_{G}.{name}.sorted.gff3.gz', required=False)
    if not p: continue
    for l in gzip.open(p, 'rt'):
        if l.startswith('#'): continue
        f = l.split('\t')
        if len(f) > 4 and f[2] not in ('gene', 'region'): other['ncRNA'][f[0]].append((int(f[3]) - 1, int(f[4])))
for p in glob.glob(os.path.join(func_dir(FUNC, G, '_TR'), 'GMATA_SSR', f'fClaHyb_{G}.SSR.*.bed.gz')):
    for l in gzip.open(p, 'rt'):
        f = l.split('\t'); other['SSR'][f[0]].append((int(f[1]), int(f[2])))
for d in other.values():
    for v in d.values(): v.sort()
known = {r['family'] for r in csv.DictReader(open(f'{INS}/54c_ins_families_classified.tsv'), delimiter='\t') if r['class'] == 'known_TE'}
kte = collections.defaultdict(list)
for r in csv.DictReader(open(f'{INS}/54b_ins_members.tsv'), delimiter='\t'):
    if r['family'] in known and r['ins_id'].startswith(G + '|'):
        c, se = r['ins_id'].split('|')[2].split(':'); s, e = map(int, se.split('-')); kte[c].append((s, e))
for v in kte.values(): v.sort()


def own(x, s, e):                      # whole element inside the gap: counts as TE even with an intrinsic label
    L = x[1] - x[0]; return L >= MINLEN and (min(e, x[1]) - max(s, x[0])) >= OWN * L


def cover(iv, c, s, e, skip=None):
    v = iv.get(c, []); i = bisect.bisect_left(v, (s - 200000,)); pts = []
    for x in v[i:]:
        if x[0] >= e: break
        if x[1] > s and (skip is None or len(x) < 3 or x[2] not in skip or own(x, s, e)): pts.append((max(s, x[0]), min(e, x[1])))
    pts.sort(); t = 0; cur = None
    for a_, b in pts:
        if cur and a_ <= cur[1]: cur[1] = max(cur[1], b)
        else:
            if cur: t += cur[1] - cur[0]
            cur = [a_, b]
    return t + (cur[1] - cur[0] if cur else 0)


by = collections.defaultdict(list)
for r in runs: by[(r['chrom'], FAM(r))].append(r)
pairs = []                                                   # (a, b) consecutive runs of one family on one chromosome
for key, rs in by.items():
    rs.sort(key=lambda r: int(r['start0']))
    for x, y in zip(rs, rs[1:]): pairs.append((x, y))


def classify(a_, b, cut, icut, maxgap):
    s, e = int(a_['end']), int(b['start0']); c = a_['chrom']
    if e <= s or e - s > maxgap: return None
    if GAPS is not None and GAPS.between(c, s, e, closed=True): return 'refused'
    skip = intr[icut].get(a_['run_id'], set()) | intr[icut].get(b['run_id'], set()); L = e - s
    t = cover(te, c, s, e, skip); k = cover({c: kte.get(c, [])}, c, s, e)
    ex = cover(other['exon_CDS'], c, s, e); gn = cover(other['gene'], c, s, e); nc = cover(other['ncRNA'], c, s, e); ss = cover(other['SSR'], c, s, e); intron = max(0, gn - ex)
    oth = max((('exon_CDS', ex), ('intron', intron), ('ncRNA', nc), ('SSR', ss)), key=lambda kv: kv[1])
    cls = 'TE_insertion' if t >= cut * L else 'missed_TE_54' if k >= cut * L else f'other:{oth[0]}' if oth[1] >= cut * L else 'unannotated'
    top = collections.Counter()
    for x in te.get(c, []):
        if x[0] < e and x[1] > s and (x[2] not in skip or own(x, s, e)): top[x[2]] += min(e, x[1]) - max(s, x[0])
    return [G, FAM(a_), a_['array_id'], b['array_id'], 'same_array' if a_['array_id_v46'] == b['array_id_v46'] else 'across_arrays', c, s, e, L, a_['run_id'], b['run_id'],
            round(t / L, 3), round(k / L, 3), round(ex / L, 3), round(intron / L, 3), round(nc / L, 3), round(ss / L, 3), cls, top.most_common(1)[0][0] if top else 'none']


HDR = ['subgenome', 'family', 'array_left', 'array_right', 'array_relation', 'chrom', 'start0', 'end', 'bp', 'run_left', 'run_right', 'TE_frac', 'missed_TE54_frac', 'exon_CDS_frac',
       'intron_frac', 'ncRNA_frac', 'SSR_frac', 'class', 'top_TE_label']
icut0 = a.intrinsic_cut if a.intrinsic_cut in intr else sorted(intr)[0]
out = []; refused = 0; S = collections.Counter(); B = collections.Counter()
for x, y in pairs:
    r = classify(x, y, a.class_cut, icut0, MAXGAP)
    if r is None: continue
    if r == 'refused': refused += 1; continue
    out.append(r); S[r[17]] += 1; B[r[17]] += r[8]
with open(f'{C}/59e_interruptions_{G}{SUF}.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n'); w.writerow(HDR); w.writerows(out)
print(G, 'interruptions between consecutive same-family satellite runs (<= %d bp, own %.2f, minlen %d, class cut %.2f, intrinsic cut %s):' % (MAXGAP, OWN, MINLEN, a.class_cut, icut0),
      {k: (S[k], round(B[k] / 1e6, 3)) for k in sorted(S)}, '| refused at assembly gaps / joins:', refused if GAPS is not None else 'no gap BED given',
      '| across arrays:', sum(1 for r in out if r[4] == 'across_arrays'))
tl = collections.Counter(x[18] for x in out if x[17] == 'TE_insertion'); print('   top TE labels:', tl.most_common(6))
if a.sweep:
    with open(f'{C}/59e_sweep_{G}{SUF}.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n'); w.writerow(['class_cut', 'intrinsic_cut', 'max_gap', 'n', 'Mb', 'class', 'n_class', 'Mb_class'])
        for cut in (0.3, 0.5, 0.7):
            for icut in sorted(intr):
                for mg in (10000, 20000, 50000):
                    Sx = collections.Counter(); Bx = collections.Counter(); n = 0; t = 0
                    for x, y in pairs:
                        r = classify(x, y, cut, icut, mg)
                        if r is None or r == 'refused': continue
                        n += 1; t += r[8]; Sx[r[17]] += 1; Bx[r[17]] += r[8]
                    for k in sorted(Sx): w.writerow([cut, icut, mg, n, round(t / 1e6, 3), k, Sx[k], round(Bx[k] / 1e6, 3)])
    print('   sweep written')
