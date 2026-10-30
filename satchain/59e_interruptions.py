#!/usr/bin/env python3
"""59e: interruptions between neighbouring catalogue runs of the same array (gap <= --max-gap), classified without EDTA's satellite labels.
For each gap: EDTA features whose label is intrinsic to the flanking runs (59c intrinsic_labels at cut 0.3) are ignored; the rest give
TE coverage. Class: TE_insertion (TE >= 50 % of the gap), missed_TE_54 (a step-54 known-TE insertion family covers >= 50 %), other_annotation
(gene / ncRNA / SSR >= 50 %), unannotated. Output 59e_interruptions_<G>.tsv and a summary.
v2 (2026-10-21): the intrinsic-label filter is per feature, not per label. An EDTA feature with a label intrinsic to the flanking runs still
counts as TE when it is a whole element sitting in the gap (>= OWN of its length inside the gap and >= MINLEN bp): EDTA gives MacSat03 and a
7.5-kb LTR element interleaved with it the same label (LTR/unknown), and the v1 filter dropped the element.
v4.7 (2026-10-26): no interruption across an assembly gap (env SATCHAIN_GAPS = gap BED); reported as refused. Renamed annotation folders handled.
usage: 59e_interruptions.py <G> <cat59 dir> <func folder> <ins54 dir> [max_gap] [own_frac] [min_len] [out suffix]"""
import sys, csv, gzip, collections, bisect, os
G, C, FUNC, INS = sys.argv[1:5]; MAXGAP = int(sys.argv[5]) if len(sys.argv) > 5 else 50000
OWN = float(sys.argv[6]) if len(sys.argv) > 6 else 0.8; MINLEN = int(sys.argv[7]) if len(sys.argv) > 7 else 500; SUF = sys.argv[8] if len(sys.argv) > 8 else ''
runs = [r for r in csv.DictReader(open(f'{C}/59_catalogue_{G}.tsv'), delimiter='\t')]
intr = {}
for r in csv.DictReader(open(f'{C}/59c_te_{G}.tsv'), delimiter='\t'):
    intr[r['run_id']] = set(x for x in r['intrinsic_labels_0.3'].split(',') if x != 'none')
te = collections.defaultdict(list)
def func_path(sub, name):
    for d in (os.path.join(FUNC, sub), os.path.join(FUNC, f'fClaHyb_{G}{sub}')):
        if os.path.exists(os.path.join(d, name)): return os.path.join(d, name)
    return os.path.join(FUNC, sub, name)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from gaps import Gaps
GAPS = Gaps(os.environ['SATCHAIN_GAPS']) if os.environ.get('SATCHAIN_GAPS') else None; refused = 0
for l in gzip.open(func_path('_repeats_EDTA', f'fClaHyb_{G}.TEanno.sorted.gff3.gz'), 'rt'):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] in ('repeat_region', 'target_site_duplication', 'long_terminal_repeat'): continue
    te[f[0]].append((int(f[3]) - 1, int(f[4]), dict(x.split('=', 1) for x in f[8].split(';') if '=' in x).get('classification', f[2])))
for v in te.values(): v.sort()
other = collections.defaultdict(list)
for sub, name in (('_genes_lifton', 'lifton'), ('_ncRNA_rfam', 'rfam.all_classes'), ('_tRNA_tRNAscan', 'tRNA'), ('_rRNA_barrnap', 'rRNA')):
    p = func_path(sub, f'fClaHyb_{G}.{name}.sorted.gff3.gz')
    if not os.path.exists(p): continue
    for l in gzip.open(p, 'rt'):
        if l.startswith('#'): continue
        f = l.split('\t')
        if len(f) > 4 and f[2] in ('exon', 'CDS', 'tRNA', 'rRNA', 'snRNA', 'ncRNA', 'gene', 'sequence_feature', 'snoRNA'): other[f[0]].append((int(f[3]) - 1, int(f[4])))
for v in other.values(): v.sort()
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
    for a, b in pts:
        if cur and a <= cur[1]: cur[1] = max(cur[1], b)
        else:
            if cur: t += cur[1] - cur[0]
            cur = [a, b]
    return t + (cur[1] - cur[0] if cur else 0)
by = collections.defaultdict(list)
for r in runs: by[r['array_id']].append(r)
out = []; S = collections.Counter(); B = collections.Counter()
for aid, rs in by.items():
    rs.sort(key=lambda r: int(r['start0']))
    for a, b in zip(rs, rs[1:]):
        s, e = int(a['end']), int(b['start0']); c = a['chrom']
        if e <= s or e - s > MAXGAP: continue
        if GAPS is not None and GAPS.between(c, s, e): refused += 1; continue           # an assembly gap is not an interruption
        skip = intr.get(a['run_id'], set()) | intr.get(b['run_id'], set()); L = e - s
        t = cover(te, c, s, e, skip); k = cover({c: kte.get(c, [])}, c, s, e); o = cover(other, c, s, e)
        cls = 'TE_insertion' if t >= 0.5 * L else 'missed_TE_54' if k >= 0.5 * L else 'other_annotation' if o >= 0.5 * L else 'unannotated'
        top = collections.Counter()
        for x in te.get(c, []):
            if x[0] < e and x[1] > s and (x[2] not in skip or own(x, s, e)): top[x[2]] += min(e, x[1]) - max(s, x[0])
        out.append([G, aid, c, s, e, L, a['run_id'], b['run_id'], round(t / L, 3), round(k / L, 3), round(o / L, 3), cls, top.most_common(1)[0][0] if top else 'none'])
        S[cls] += 1; B[cls] += L
with open(f'{C}/59e_interruptions_{G}{SUF}.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['subgenome', 'array_id', 'chrom', 'start0', 'end', 'bp', 'run_left', 'run_right', 'TE_frac', 'missed_TE54_frac', 'other_annot_frac', 'class', 'top_TE_label']); w.writerows(out)
print(G, 'interruptions between runs (<= %d bp, own %.2f, minlen %d):' % (MAXGAP, OWN, MINLEN), {k: (S[k], round(B[k] / 1e6, 2)) for k in S},
      '| refused at assembly gaps:', refused if GAPS is not None else 'no gap BED given')
tl = collections.Counter(x[12] for x in out if x[11] == 'TE_insertion'); print('   top TE labels:', tl.most_common(6))
