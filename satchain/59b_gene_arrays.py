#!/usr/bin/env python3
"""59b: tandem gene arrays in the catalogue (module used by 59_final_catalogue.py; also a stand-alone report).
Gene tracks (annotation_overlays/<G>): 5S rRNA and other rRNA (barrnap), tRNA (tRNAscan), miRNA (curated; miR-430 apart), and the Rfam
hits of every class merged symmetrically into the same classes (2026-11-01, audit C13: 5S RF00001 -> 5S_rRNA, tRNA RF00005 -> tRNA,
LSU/SSU -> rRNA_other, U1-U6 -> snRNA_<U>, other snRNA -> snRNA_other, mir-430 -> miR-430, other miRNA -> miRNA_other). Overlapping
features of one class are merged (one gene counted once even when two tracks report it).
Per family, for its top class: (1) bp ratio = observed bp share / expected share with one gene per unit, where the gene length is the
median merged-feature length of that class in this genome (not a constant; audit C12) and the unit is the bp-weighted median unit of the
family's runs; (2) unit share = merged features overlapping the family's runs / (family bp / unit), capped at 1 (the per-unit criterion
of audit C01/C12, insensitive to partial hits). A family is gene_array:<class> when ratio >= cut (0.7) AND unit share >= unit_cut (0.5);
gene_related:<class> when ratio >= 0.2 (either criterion short; the family STAYS element_class satellite and carries the flag:
tRNA-derived minisatellites such as MacSat1272 fall here); satellite otherwise. gene_flag = 'gene_array:<class>' or
'gene_related:<class>' (audit C25: the flag is set for the gene arrays themselves too). The sweep table gives the Mb per call at
cut 0.5/0.7/0.8 x unit_cut 0.3/0.5/0.7.
usage (stand-alone): 59b_gene_arrays.py <G> <catalogue tsv> <overlay dir (annotation_overlays/G)> <out tsv>"""
import sys, os, csv, gzip, collections, bisect, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from satpaths import overlay_path
CUTS = (0.5, 0.7, 0.8); UCUTS = (0.3, 0.5, 0.7)


def load_tracks(G, OVL):
    feats = collections.defaultdict(list)
    def add(c, s, e, lab): feats[lab].append((c, s, e))
    for l in gzip.open(overlay_path(OVL, G, f'fClaHyb_{G}.rRNA.sorted.gff3.gz'), 'rt'):
        if l.startswith('#'): continue
        f = l.split('\t'); add(f[0], int(f[3]) - 1, int(f[4]), '5S_rRNA' if '5S' in f[8] else 'rRNA_other')
    for l in gzip.open(overlay_path(OVL, G, f'fClaHyb_{G}.tRNA.sorted.gff3.gz'), 'rt'):
        if l.startswith('#'): continue
        f = l.split('\t')
        if f[2] in ('tRNA', 'gene', 'pseudogene'): add(f[0], int(f[3]) - 1, int(f[4]), 'tRNA')
    for l in gzip.open(overlay_path(OVL, G, f'fClaHyb_{G}.miRNA.curated.sorted.gff3.gz'), 'rt'):
        if l.startswith('#'): continue
        f = l.split('\t'); add(f[0], int(f[3]) - 1, int(f[4]), 'miR-430' if re.search(r'mir-?430', f[8], re.I) else 'miRNA_other')
    for l in gzip.open(overlay_path(OVL, G, f'fClaHyb_{G}.rfam.all_classes.sorted.gff3.gz'), 'rt'):
        if l.startswith('#'): continue
        f = l.split('\t'); at = dict(x.split('=', 1) for x in f[8].strip().split(';') if '=' in x); nm = at.get('Name', ''); acc = at.get('rfam_acc', '') + at.get('Dbxref', '') + f[8]
        m = re.match(r'^(U[1-6])(?:atac)?$', nm) or re.match(r'^(U[1-6])_', nm)
        if f[2] == 'snRNA' or m: add(f[0], int(f[3]) - 1, int(f[4]), f'snRNA_{m.group(1)}' if m else 'snRNA_other')
        elif re.search(r'5S_rRNA|RF00001', acc): add(f[0], int(f[3]) - 1, int(f[4]), '5S_rRNA')
        elif re.search(r'\btRNA\b|RF00005', acc): add(f[0], int(f[3]) - 1, int(f[4]), 'tRNA')
        elif re.search(r'LSU_rRNA|SSU_rRNA|5_8S_rRNA|RF0000[12]\d|RF02541|RF01960|RF00002', acc): add(f[0], int(f[3]) - 1, int(f[4]), 'rRNA_other')
        elif re.search(r'mir-?430', nm, re.I): add(f[0], int(f[3]) - 1, int(f[4]), 'miR-430')
        elif re.search(r'^mir-|^MIR', nm, re.I): add(f[0], int(f[3]) - 1, int(f[4]), 'miRNA_other')
    idx = {}; glen = {}
    for lab, v in feats.items():
        d = collections.defaultdict(list)
        for c, s, e in v: d[c].append((s, e))
        lens = []
        for c in d:
            d[c].sort(); mg = []
            for x in d[c]:
                if mg and x[0] <= mg[-1][1]: mg[-1] = (mg[-1][0], max(mg[-1][1], x[1]))
                else: mg.append(x)
            d[c] = mg; lens += [e - s for s, e in mg]
        idx[lab] = d; lens.sort(); glen[lab] = lens[len(lens) // 2] if lens else 100
    return idx, glen


def ov(d, c, s, e):
    """(overlap bp, n features overlapping) of [s, e) with the merged features of chromosome c"""
    v = d.get(c, []); i = bisect.bisect_left(v, (s - 100000, 0)); t = 0; n = 0
    for x in v[i:]:
        if x[0] >= e: break
        o = max(0, min(e, x[1]) - max(s, x[0]))
        if o: t += o; n += 1
    return t, n


def classify(rows, G, OVL, cut=0.7, unit_cut=0.5):
    """rows: catalogue dicts with family, chrom, start0, end, unit_bp. Returns ({family: call dict}, sweep rows)."""
    idx, glen = load_tracks(G, OVL)
    fam_bp = collections.Counter(); nr = collections.Counter(); units = collections.defaultdict(list)
    hit_bp = collections.defaultdict(collections.Counter); hit_n = collections.defaultdict(collections.Counter); hit_runs = collections.defaultdict(collections.Counter)
    for r in rows:
        f = r['family']; c, s, e = r['chrom'], int(r['start0']), int(r['end']); fam_bp[f] += e - s; nr[f] += 1; units[f].append((int(r['unit_bp']), e - s))
        for lab, d in idx.items():
            o, n = ov(d, c, s, e)
            if o: hit_bp[f][lab] += o; hit_n[f][lab] += n; hit_runs[f][lab] += 1
    calls = {}; sweep = collections.defaultdict(collections.Counter)
    for f, b in fam_bp.most_common():
        if not hit_bp[f]: continue
        lab, o = hit_bp[f].most_common(1)[0]; bs = o / b; rs = hit_runs[f][lab] / nr[f]
        us = sorted(units[f]); acc = 0; u = us[-1][0]
        for uu, w in us:                                   # bp-weighted median unit
            acc += w
            if acc >= b / 2: u = uu; break
        exp = min(1.0, glen[lab] / u) if u else 1.0; ratio = bs / exp if exp else 0.0
        n_units = b / u if u else 0; ushare = min(1.0, hit_n[f][lab] / n_units) if n_units else 0.0
        def call(c_, uc):
            return f'gene_array:{lab}' if ratio >= c_ and ushare >= uc else f'gene_related:{lab}' if ratio >= 0.2 else 'satellite_sporadic_genes'
        for c_ in CUTS:
            for uc in UCUTS: sweep[(c_, uc)][call(c_, uc).split(':')[0]] += b
        cl = call(cut, unit_cut)
        calls[f] = dict(family=f, bp=b, n_runs=nr[f], unit_bp_median_bpw=u, top_class=lab, top_bp_share=round(bs, 3), gene_len_median=glen[lab], expected_share_one_gene_per_unit=round(exp, 3),
                        observed_over_expected=round(ratio, 3), hits_per_unit=round(hit_n[f][lab] / n_units, 3) if n_units else 0.0, unit_share=round(ushare, 3), top_run_share=round(rs, 3),
                        element_class=cl if cl.startswith('gene_array') else 'satellite', gene_flag=cl if cl != 'satellite_sporadic_genes' else 'none',   # gene_related families stay satellite, flagged (C25)
                        **{f'call_at_{c_}_u{uc}': call(c_, uc) for c_ in CUTS for uc in UCUTS}, all_classes=';'.join(f'{k}:{v / b:.3f}' for k, v in hit_bp[f].most_common()))
    sw = [dict(ratio_cut=c_, unit_cut=uc, **{k: round(v / 1e6, 3) for k, v in sweep[(c_, uc)].items()}) for c_ in CUTS for uc in UCUTS]
    return calls, sw


def write(calls, sweep, out):
    if calls:
        cols = list(next(iter(calls.values())).keys())
        with open(out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(calls.values())
    with open(out.replace('.tsv', '_sweep.tsv'), 'w', newline='') as fh:
        keys = ['ratio_cut', 'unit_cut', 'gene_array', 'gene_related', 'satellite_sporadic_genes']
        w = csv.DictWriter(fh, fieldnames=keys, delimiter='\t', lineterminator='\n', restval=0.0); w.writeheader(); w.writerows(sweep)


if __name__ == '__main__':
    G, CAT, OVL, OUT = sys.argv[1:5]
    rows = list(csv.DictReader(open(CAT), delimiter='\t')); calls, sweep = classify(rows, G, OVL); write(calls, sweep, OUT)
    for s in sweep: print(G, s)
    print(G, 'gene arrays:', [(c['family'], c['top_class'], c['observed_over_expected'], c['unit_share']) for c in calls.values() if c['element_class'].startswith('gene_array')])
