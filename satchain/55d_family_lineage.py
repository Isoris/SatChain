#!/usr/bin/env python3
"""55d: family -> lineage table from the copy-level clustering (55c, k = 6, K = 3): each family gets the cluster holding the majority of
its sampled copies (A = cluster with the old A copies, B = with the old B copies, centromere-associated minisatellite = the third), with the
majority share and n copies; the OLD label is listed for comparison. usage: 55d_family_lineage.py <lin55 dir>"""
import csv, collections, sys
D = sys.argv[1]
rows = list(csv.DictReader(open(f'{D}/55_cluster_assignments.tsv'), delimiter='\t'))
cnt = collections.defaultdict(collections.Counter)
for r in rows: cnt[r['k6_K3']][r['old_lineage']] += 1
name = {c: ('B' if v['B'] > 100 else 'A' if v['A'] > 100 else 'cen_minisatellite') for c, v in cnt.items()}
fam = collections.defaultdict(collections.Counter); old = {}
for r in rows:
    k = (r['subgenome'], r['family']); fam[k][name[r['k6_K3']]] += 1; old[k] = r['old_lineage']
with open(f'{D}/55d_family_lineage.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n'); w.writerow(['subgenome', 'family', 'lineage', 'majority_share', 'n_sampled_copies', 'old_lineage', 'change'])
    ch = collections.Counter()
    for (g, f), c in sorted(fam.items()):
        lin, n = c.most_common(1)[0]; tot = sum(c.values()); o = old[(g, f)]
        change = 'same' if o == lin else ('newly_assigned' if o == 'none' else f'{o}->{lin}')
        ch[change] += 1; w.writerow([g, f, lin, round(n / tot, 3), tot, o, change])
print(dict(ch))
