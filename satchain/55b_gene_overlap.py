#!/usr/bin/env python3
"""55b (2026-11-01): independent evidence for the clusters of 55c, from annotation coordinates instead of a BLAST of the sampled copies
(the old 55b was a comment in the workflow and its output was not on disk, audit C09). For every sampled copy: the fraction of its bp
covered by each gene class of the merged tracks of 59b (5S_rRNA, rRNA_other, tRNA, snRNA_*, miR-430, miRNA_other). Writes
<dest>/55_sample_gene_overlap.tsv (copy_id, frac_5S, top_class, top_frac). The 5S fraction replaces the old 'score_5S' (0-100 scale
pid x coverage): frac_5S = 1 means the whole monomer is annotated 5S.
usage: 55b_gene_overlap.py <lin55 dir> <G>=<annotation_overlays/G> ..."""
import sys, os, csv, importlib.util
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location('ga', os.path.join(os.path.dirname(os.path.abspath(__file__)), '59b_gene_arrays.py')); ga = importlib.util.module_from_spec(spec); spec.loader.exec_module(ga)
D = sys.argv[1]; tracks = {}
for spec_ in sys.argv[2:]:
    G, OVL = spec_.split('=', 1); tracks[G] = ga.load_tracks(G, OVL)[0]
with open(os.path.join(D, '55_sample_gene_overlap.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n'); w.writerow(['copy_id', 'frac_5S', 'top_class', 'top_frac']); n = 0; n5 = 0
    for m in csv.DictReader(open(os.path.join(D, '55_sample.tsv')), delimiter='\t'):
        idx = tracks[m['subgenome']]; c, s, e = m['chrom'], int(m['start0']), int(m['end']); L = max(1, e - s)
        fr = {lab: ga.ov(d, c, s, e)[0] / L for lab, d in idx.items()}; top = max(fr.items(), key=lambda kv: kv[1]) if fr else ('none', 0.0)
        w.writerow([m['copy_id'], round(fr.get('5S_rRNA', 0.0), 3), top[0] if top[1] > 0 else 'none', round(top[1], 3)]); n += 1; n5 += fr.get('5S_rRNA', 0.0) >= 0.3
print(n, 'copies scored;', n5, 'with >= 30 % of the monomer annotated 5S')
