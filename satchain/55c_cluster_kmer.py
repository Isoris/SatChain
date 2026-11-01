#!/usr/bin/env python3
"""55c (2026-11-01; replaces 55c_lineage_cluster.py): unsupervised SEQUENCE CLUSTERS of satellite monomer copies. What this step finds is
a partition of monomers by k-mer composition; it is not a phylogeny and the clusters are not evolutionary lineages (the old name 'lineage'
is dropped everywhere). Origin statements (5S-derived, RTE-interleaved) are made per family from homology evidence, not from here.
Features: canonical k-mer frequencies per copy (orientation- and rotation-free), k = 5, 6, 7; Hellinger transform, centring, SVD to 20
dimensions; k-means (K = 2..8, 10 restarts, seeded). Reported per (k, K): silhouette (2,000-copy subsample), ARI vs the k = 6 partition
(stability across k), ARI between two independent restart seeds (stability across restarts), and the same for each haplotype alone
(the pooled labels restricted to the haplotype vs a clustering of that haplotype only). K is chosen by a stated rule
(55_K_choice.tsv): the largest K whose minimum cross-k ARI and restart ARI are both >= --ari (0.8) and whose silhouette is within 10 %
of the best silhouette over K; when no K >= 3 passes, K = 2 (always reported) and the choice is marked 'fallback'.
Clusters are named C<K>.<rank by size> and described by content (n, haplotype split, median monomer length, median 5S fraction, share of
copies with >= 30 % 5S, top families), never by old labels (audit C18). Independent evidence comes from 55b (55_sample_gene_overlap.tsv).
Outputs: 55_cluster_assignments.tsv (copy x (k, K)), 55_cluster_summary.tsv, 55_K_choice.tsv, 55_haplotype_stability.tsv.
usage: 55c_cluster_kmer.py <lin55 dir> [--kmax 8] [--ari 0.8]"""
import csv, collections, sys, os, argparse
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('D'); ap.add_argument('--kmax', type=int, default=8); ap.add_argument('--ari', type=float, default=0.8); a = ap.parse_args(); D = a.D
meta = list(csv.DictReader(open(f'{D}/55_sample.tsv'), delimiter='\t'))
seqs = {}; name = None
for l in open(f'{D}/55_sample.fa'):
    if l.startswith('>'): name = l[1:].strip()
    else: seqs[name] = l.strip()
f5 = {}; p = f'{D}/55_sample_gene_overlap.tsv'
if os.path.exists(p):
    for r in csv.DictReader(open(p), delimiter='\t'): f5[r['copy_id']] = float(r['frac_5S'])
ids = [m['copy_id'] for m in meta]; COMP = str.maketrans('ACGT', 'TGCA')


def feats(k):
    idx = {}; rows = []
    for cid in ids:
        s = seqs[cid]; c = collections.Counter()
        for i in range(len(s) - k + 1):
            w = s[i:i + k]
            if 'N' in w: continue
            r = w.translate(COMP)[::-1]; c[min(w, r)] += 1
        rows.append(c)
        for w in c: idx.setdefault(w, len(idx))
    X = np.zeros((len(ids), len(idx)), dtype=np.float32)
    for i, c in enumerate(rows):
        t = sum(c.values()) or 1
        for w, v in c.items(): X[i, idx[w]] = v / t
    X = np.sqrt(X); X -= X.mean(0)
    from sklearn.decomposition import TruncatedSVD                 # randomized 20-dim SVD: the full SVD of 30k x 8k (k = 7) ran out of memory
    return TruncatedSVD(n_components=20, random_state=55, algorithm='randomized', n_iter=7).fit_transform(X).astype(np.float64)


def kmeans(Z, K, rng, restarts=10):
    best = None
    for r in range(restarts):
        C = Z[rng.choice(len(Z), K, replace=False)]
        for _ in range(100):
            d = ((Z[:, None, :] - C[None]) ** 2).sum(2); lab = d.argmin(1)
            C2 = np.array([Z[lab == j].mean(0) if (lab == j).any() else C[j] for j in range(K)])
            if np.allclose(C2, C): break
            C = C2
        inertia = d[np.arange(len(Z)), lab].sum()
        if best is None or inertia < best[0]: best = (inertia, lab)
    lab = best[1]; order = {j: r for r, (j, _) in enumerate(collections.Counter(lab.tolist()).most_common())}      # rank by size -> stable names
    return np.array([order[j] for j in lab])


def silhouette(Z, lab, rng, n=2000):
    i = rng.choice(len(Z), min(n, len(Z)), replace=False); Zs, ls = Z[i], lab[i]
    D_ = np.sqrt(((Zs[:, None] - Zs[None]) ** 2).sum(2)); s = []
    for q in range(len(Zs)):
        same = ls == ls[q]; same[q] = False
        if not same.any() or len(set(ls)) < 2: continue
        ai = D_[q, same].mean(); bi = min(D_[q, ls == k].mean() for k in set(ls) if k != ls[q])
        s.append((bi - ai) / max(ai, bi))
    return float(np.mean(s)) if s else float('nan')


def ari(x, y):
    n = len(x); comb = lambda v: v * (v - 1) / 2
    nij = collections.Counter(zip(x, y)); A = collections.Counter(x); B = collections.Counter(y)
    s = sum(comb(v) for v in nij.values()); sa = sum(comb(v) for v in A.values()); sb = sum(comb(v) for v in B.values())
    e = sa * sb / comb(n); m = (sa + sb) / 2
    return (s - e) / (m - e) if m != e else 1.0


KS = (5, 6, 7); Krange = range(2, a.kmax + 1); sub = np.array([m['subgenome'] for m in meta]); L = np.array([int(m['length']) for m in meta]); fam = np.array([m['family'] for m in meta])
F5 = np.array([f5.get(c, 0.0) for c in ids])
Zk = {k: feats(k) for k in KS}; labs = {}; labs2 = {}; summ = []; hap = []
for k, Z in Zk.items():
    for K in Krange:
        labs[(k, K)] = kmeans(Z, K, np.random.default_rng(55 + k)); labs2[(k, K)] = kmeans(Z, K, np.random.default_rng(1000 + 55 + k))
        summ.append(dict(k=k, K=K, silhouette=round(silhouette(Z, labs[(k, K)], np.random.default_rng(7)), 3), ARI_restarts=round(ari(labs[(k, K)].tolist(), labs2[(k, K)].tolist()), 3),
                         sizes=','.join(str(int((labs[(k, K)] == j).sum())) for j in range(K))))
        for g in sorted(set(sub)):                      # per-haplotype partition vs the pooled labels restricted to that haplotype
            m = sub == g
            if m.sum() < 2 * K: continue
            lg = kmeans(Z[m], K, np.random.default_rng(55 + k)); hap.append(dict(k=k, K=K, haplotype=g, n=int(m.sum()), ARI_pooled_vs_own=round(ari(labs[(k, K)][m].tolist(), lg.tolist()), 3),
                                                                                 n_clusters_present_in_pooled=len(set(labs[(k, K)][m].tolist()))))
for r in summ: r['ARI_vs_k6'] = round(ari(labs[(r['k'], r['K'])].tolist(), labs[(6, r['K'])].tolist()), 3)
# K choice rule
byK = collections.defaultdict(list)
for r in summ: byK[r['K']].append(r)
best_sil = max(np.nanmean([r['silhouette'] for r in byK[K]]) for K in Krange); choice = []
for K in Krange:
    rs = byK[K]; min_ari = min(r['ARI_vs_k6'] for r in rs); min_rst = min(r['ARI_restarts'] for r in rs); sil = float(np.nanmean([r['silhouette'] for r in rs]))
    ok = min_ari >= a.ari and min_rst >= a.ari and sil >= 0.9 * best_sil
    choice.append(dict(K=K, min_ARI_across_k=round(min_ari, 3), min_ARI_restarts=round(min_rst, 3), mean_silhouette=round(sil, 3), passes_rule=ok))
passing = [c['K'] for c in choice if c['passes_rule'] and c['K'] >= 3]; Kstar = max(passing) if passing else 2
for c in choice: c['chosen'] = 'yes' if c['K'] == Kstar else 'no'; c['note'] = ('fallback_K2' if not passing else 'rule') if c['K'] == Kstar else ''
with open(f'{D}/55_K_choice.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(choice[0].keys()), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(choice)
    fh.write(f'# rule: largest K with min ARI across k >= {a.ari}, min ARI between restart seeds >= {a.ari}, mean silhouette >= 0.9 x best ({best_sil:.3f}); fallback K = 2. Chosen K = {Kstar} at k = 6.\n')
with open(f'{D}/55_cluster_summary.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['k', 'K', 'cluster', 'silhouette', 'ARI_vs_k6', 'ARI_restarts', 'n', 'Gar', 'Mac', 'median_len', 'median_frac_5S', 'frac_copies_5S_ge0.3', 'top_families'])
    for r in summ:
        lab = labs[(r['k'], r['K'])]
        for j in range(r['K']):
            m = lab == j; tf = collections.Counter(fam[m].tolist()).most_common(3)
            w.writerow([r['k'], r['K'], f'C{r["K"]}.{j + 1}', r['silhouette'], r['ARI_vs_k6'], r['ARI_restarts'], int(m.sum()), int((sub[m] == 'Gar').sum()), int((sub[m] == 'Mac').sum()),
                        int(np.median(L[m])), round(float(np.median(F5[m])), 3), round(float((F5[m] >= 0.3).mean()), 3), ';'.join(f'{f}:{n}' for f, n in tf)])
with open(f'{D}/55_haplotype_stability.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(hap[0].keys()), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(hap)
with open(f'{D}/55_cluster_assignments.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n'); keys = sorted(labs)
    w.writerow(['copy_id', 'subgenome', 'run_id', 'array_id', 'family', 'length', 'frac_5S', 'chosen_K', 'cluster_k6_Kchosen'] + [f'k{k}_K{K}' for k, K in keys])
    for i, m in enumerate(meta): w.writerow([m['copy_id'], m['subgenome'], m['run_id'], m['array_id'], m['family'], m['length'], round(F5[i], 3), Kstar, f'C{Kstar}.{labs[(6, Kstar)][i] + 1}'] + [int(labs[x][i]) for x in keys])
for c in choice: print(c)
print('chosen K', Kstar)
