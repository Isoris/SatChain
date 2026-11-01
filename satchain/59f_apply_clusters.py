#!/usr/bin/env python3
"""59f (2026-11-01): write the sequence clusters of 55d into an existing catalogue (seq_cluster column, BED colours, family tiers table)
without rerunning 59. Breaks the old circular order 59 -> 53 -> 55 -> 59: 59 runs first with '-' as cluster table, 55a-d sample from
its runs, 59f adds the result. Gene arrays keep 'na'.
usage: 59f_apply_clusters.py <G> <cat59 dir> <55d_family_cluster.tsv>"""
import sys, csv, collections
G, C, T = sys.argv[1:4]
clus = {r['family']: r['seq_cluster'] for r in csv.DictReader(open(T), delimiter='\t') if r['subgenome'] == G}
rows = list(csv.DictReader(open(f'{C}/59_catalogue_{G}.tsv'), delimiter='\t')); n = 0
for r in rows:
    new = clus.get(r['family'], 'none') if r['element_class'] == 'satellite' else 'na'
    n += new != r['seq_cluster']; r['seq_cluster'] = new
with open(f'{C}/59_catalogue_{G}.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
PAL = ['31,119,180', '214,39,40', '148,103,189', '44,160,44', '255,127,14', '140,86,75', '227,119,194', '23,190,207']
names = sorted({r['seq_cluster'] for r in rows if r['seq_cluster'] not in ('none', 'na', 'unassigned')}); COL = {c: PAL[i % len(PAL)] for i, c in enumerate(names)}
with open(f'{C}/59_catalogue_{G}.bed', 'w') as fh:
    for r in rows:
        fh.write(f"{r['chrom']}\t{r['start0']}\t{r['end']}\t{r['run_id']}|{r['family']}|{r['element_class']}|{r['seq_cluster']}|{r['unit_bp']}bp|{r['tier']}\t0\t.\t{r['start0']}\t{r['end']}\t{COL.get(r['seq_cluster'], '150,150,150')}\n")
tiers = list(csv.DictReader(open(f'{C}/59_family_tiers_{G}.tsv'), delimiter='\t'))
for t in tiers: t['seq_cluster'] = clus.get(t['family'], 'none') if t['element_class'] == 'satellite' else 'na'
with open(f'{C}/59_family_tiers_{G}.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(tiers[0].keys()), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(tiers)
S = collections.Counter(); B = collections.Counter()
for r in rows: S[r['seq_cluster']] += 1; B[r['seq_cluster']] += int(r['bp'])
print(G, n, 'runs changed;', {k: (S[k], round(B[k] / 1e6, 3)) for k in sorted(S)})
