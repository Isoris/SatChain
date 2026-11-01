#!/usr/bin/env python3
"""55d (2026-11-01; replaces 55d_family_lineage.py): family -> sequence cluster from the copy-level clustering of 55c at the chosen K
(k = 6; 55_K_choice.tsv). A family gets the cluster that holds the majority of its sampled copies when that majority is >= --min-share
(0.7) and >= --min-copies (5) copies were sampled; otherwise 'unassigned' (audit C05: no plurality calls, no '_low_n' suffix).
Cluster names are content-based (C<K>.<rank>) with a descriptor column (haplotype split, median length, 5S fraction) from 55_cluster_summary.tsv.
Output: 55d_family_cluster.tsv (subgenome, family, seq_cluster, majority_share, n_sampled_copies, n_runs_sampled, cluster_descriptor)
usage: 55d_family_cluster.py <lin55 dir> [--min-share 0.7] [--min-copies 5] [--K <override>]"""
import csv, collections, sys, argparse
ap = argparse.ArgumentParser(); ap.add_argument('D'); ap.add_argument('--min-share', type=float, default=0.7); ap.add_argument('--min-copies', type=int, default=5); ap.add_argument('--K', type=int, default=0)
a = ap.parse_args(); D = a.D
rows = list(csv.DictReader(open(f'{D}/55_cluster_assignments.tsv'), delimiter='\t'))
K = a.K or int(rows[0]['chosen_K']); col = f'k6_K{K}'
desc = {}
for r in csv.DictReader(open(f'{D}/55_cluster_summary.tsv'), delimiter='\t'):
    if int(r['k']) == 6 and int(r['K']) == K: desc[r['cluster']] = f"n={r['n']};Gar={r['Gar']};Mac={r['Mac']};len={r['median_len']};5S={r['median_frac_5S']};top={r['top_families']}"
fam = collections.defaultdict(collections.Counter); runs = collections.defaultdict(set)
for r in rows:
    k = (r['subgenome'], r['family']); fam[k][f"C{K}.{int(r[col]) + 1}"] += 1; runs[k].add(r['run_id'])
with open(f'{D}/55d_family_cluster.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n'); w.writerow(['subgenome', 'family', 'seq_cluster', 'majority_share', 'n_sampled_copies', 'n_runs_sampled', 'cluster_descriptor'])
    S = collections.Counter()
    for (g, f), c in sorted(fam.items()):
        cl, n = c.most_common(1)[0]; tot = sum(c.values()); share = n / tot
        out = cl if share >= a.min_share and tot >= a.min_copies else 'unassigned'; S[out] += 1
        w.writerow([g, f, out, round(share, 3), tot, len(runs[(g, f)]), desc.get(cl, '') if out != 'unassigned' else ''])
print('K =', K, dict(S))
