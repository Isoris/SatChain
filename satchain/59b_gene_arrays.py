#!/usr/bin/env python3
"""59b: tandem gene arrays in the catalogue. For every family: share of its run bp overlapping each gene annotation class
(5S rRNA, other rRNA (barrnap); snRNA U1 / U2 / U4 / U5 / U6 / other (Rfam); tRNA (tRNAscan); miRNA (curated), miR-430 separately),
and the share of its runs touching that class. A family is a gene array when one class covers >= --min-bp-share of its bp
(swept 0.1 / 0.2 / 0.3) or is carried by >= half of its runs; it is then relabelled 'gene_array:<class>' (not satellite).
usage: 59b_gene_arrays.py <G> <catalogue tsv> <overlay dir (annotation_overlays/G)> <out tsv>"""
import sys, csv, gzip, collections, bisect, re
G, CAT, OVL, OUT = sys.argv[1:5]
feats = collections.defaultdict(list)
def add(c, s, e, lab): feats[lab].append((c, s, e))
for l in gzip.open(f'{OVL}/fClaHyb_{G}.rRNA.sorted.gff3.gz', 'rt'):
    if l.startswith('#'): continue
    f = l.split('\t'); add(f[0], int(f[3]) - 1, int(f[4]), '5S_rRNA' if '5S' in f[8] else 'rRNA_other')
for l in gzip.open(f'{OVL}/fClaHyb_{G}.tRNA.sorted.gff3.gz', 'rt'):
    if l.startswith('#'): continue
    f = l.split('\t')
    if f[2] in ('tRNA', 'gene', 'pseudogene'): add(f[0], int(f[3]) - 1, int(f[4]), 'tRNA')
for l in gzip.open(f'{OVL}/fClaHyb_{G}.miRNA.curated.sorted.gff3.gz', 'rt'):
    if l.startswith('#'): continue
    f = l.split('\t'); add(f[0], int(f[3]) - 1, int(f[4]), 'miR-430' if re.search(r'mir-?430', f[8], re.I) else 'miRNA_other')
for l in gzip.open(f'{OVL}/fClaHyb_{G}.rfam.all_classes.sorted.gff3.gz', 'rt'):
    if l.startswith('#'): continue
    f = l.split('\t'); nm = dict(x.split('=', 1) for x in f[8].strip().split(';') if '=' in x).get('Name', '')
    m = re.match(r'^(U[1-6])(?:atac)?$', nm) or re.match(r'^(U[1-6])_', nm)
    if f[2] == 'snRNA' or m: add(f[0], int(f[3]) - 1, int(f[4]), f'snRNA_{m.group(1)}' if m else 'snRNA_other')
idx = {}
for lab, v in feats.items():
    d = collections.defaultdict(list)
    for c, s, e in v: d[c].append((s, e))
    for c in d: d[c].sort()
    idx[lab] = d
def ov(d, c, s, e):
    v = d.get(c, []); i = bisect.bisect_left(v, (s - 100000, 0)); t = 0
    for x in v[i:]:
        if x[0] >= e: break
        t += max(0, min(e, x[1]) - max(s, x[0]))
    return t
rows = list(csv.DictReader(open(CAT), delimiter='\t'))
units = collections.defaultdict(list); fam_bp = collections.Counter(); hit_bp = collections.defaultdict(collections.Counter); hit_runs = collections.defaultdict(collections.Counter); nr = collections.Counter()
for r in rows:
    f = r['family']; c, s, e = r['chrom'], int(r['start0']), int(r['end']); fam_bp[f] += e - s; nr[f] += 1; units[f].append(int(r['unit_bp']))
    for lab, d in idx.items():
        o = ov(d, c, s, e)
        if o: hit_bp[f][lab] += o; hit_runs[f][lab] += 1
with open(OUT, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    # one gene per repeat unit gives an expected bp share of gene_length / unit; the observed / expected ratio decides
    GLEN = {'5S_rRNA': 120, 'rRNA_other': 1800, 'tRNA': 73, 'miR-430': 70, 'miRNA_other': 70, 'snRNA_U1': 164, 'snRNA_U2': 190,
            'snRNA_U4': 145, 'snRNA_U5': 116, 'snRNA_U6': 107, 'snRNA_other': 150}
    CUTS = (0.5, 0.7, 0.8)
    w.writerow(['family', 'bp', 'n_runs', 'unit_bp_median', 'top_class', 'top_bp_share', 'expected_share_one_gene_per_unit', 'observed_over_expected', 'top_run_share']
               + [f'call_at_{c}' for c in CUTS] + ['all_classes'])
    out = collections.defaultdict(lambda: collections.Counter())
    for f, b in fam_bp.most_common():
        if not hit_bp[f]: continue
        lab, o = hit_bp[f].most_common(1)[0]; bs = o / b; rs = hit_runs[f][lab] / nr[f]
        u = sorted(units[f])[len(units[f]) // 2] if units[f] else 0
        exp = min(1.0, GLEN[lab] / u) if u else 1.0; ratio = bs / exp if exp else 0
        calls = [(f'gene_array:{lab}' if ratio >= c else f'gene_related:{lab}' if ratio >= 0.2 else 'satellite_sporadic_genes') for c in CUTS]
        for c, cl in zip(CUTS, calls): out[c][cl.split(':')[0]] += b
        w.writerow([f, b, nr[f], u, lab, round(bs, 3), round(exp, 3), round(ratio, 3), round(rs, 3)] + calls + [';'.join(f'{k}:{v / b:.3f}' for k, v in hit_bp[f].most_common())])
for c in CUTS: print(G, f'ratio cut {c}:', {k: f'{v / 1e6:.2f}Mb' for k, v in sorted(out[c].items())})
