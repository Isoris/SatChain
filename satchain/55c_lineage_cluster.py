#!/usr/bin/env python3
"""55c: unsupervised lineage structure of centromeric monomer copies (no family or old-lineage labels used as input).
Features: canonical k-mer frequencies of each copy (orientation-free, rotation-free), k swept (5, 6, 7); centred log-ratio-like
scaling, SVD to 20 dimensions; k-means (K = 2..6, 10 restarts, seeded); silhouette on a 2,000-copy subsample; stability = adjusted
Rand index between k values at the same K. Then compares the clusters with independent evidence (5S similarity, monomer length,
subgenome) and with the OLD family-level lineage labels.
Outputs: 55_cluster_assignments.tsv (copy x K x k), 55_cluster_summary.tsv, 55_segment_lineage.tsv (per segment composition at the
chosen K, k=6). usage: 55c_lineage_cluster.py <lin55 dir>"""
import csv, collections, math, sys
import numpy as np
D = sys.argv[1]; rng = np.random.default_rng(55)
meta = list(csv.DictReader(open(f'{D}/55_sample.tsv'), delimiter='\t'))
seqs = {}; name = None
for l in open(f'{D}/55_sample.fa'):
    if l.startswith('>'): name = l[1:].strip()
    else: seqs[name] = l.strip()
s5 = {}
for l in open(f'{D}/55_sample_vs_5S.tsv'):
    q, s, pid, ln, *_ , ql, sl = l.rstrip('\n').split('\t')
    s5[q] = max(s5.get(q, 0), float(pid) * int(ln) / int(ql))
ids = [m['copy_id'] for m in meta]
COMP = str.maketrans('ACGT', 'TGCA')


def feats(k):
    idx = {}; X = np.zeros((len(ids), 0))
    rows = []
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
    X = np.sqrt(X)                                    # Hellinger transform: composition distance without zero problems
    X -= X.mean(0)
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    return U[:, :20] * S[:20]


def kmeans(Z, K, restarts=10):
    best = None
    for r in range(restarts):
        C = Z[rng.choice(len(Z), K, replace=False)]
        for _ in range(100):
            d = ((Z[:, None, :] - C[None]) ** 2).sum(2); lab = d.argmin(1)
            C2 = np.array([Z[lab == j].mean(0) if (lab == j).any() else C[j] for j in range(K)])
            if np.allclose(C2, C): break
            C = C2
        inertia = d[np.arange(len(Z)), lab].sum()
        if best is None or inertia < best[0]: best = (inertia, lab, d)
    return best[1], best[2]


def silhouette(Z, lab, n=2000):
    i = rng.choice(len(Z), min(n, len(Z)), replace=False); Zs, ls = Z[i], lab[i]
    D_ = np.sqrt(((Zs[:, None] - Zs[None]) ** 2).sum(2)); s = []
    for a in range(len(Zs)):
        same = ls == ls[a]; same[a] = False
        if not same.any(): continue
        ai = D_[a, same].mean(); bi = min(D_[a, ls == k].mean() for k in set(ls) if k != ls[a])
        s.append((bi - ai) / max(ai, bi))
    return float(np.mean(s))


def ari(x, y):
    from collections import Counter
    n = len(x); comb = lambda v: v * (v - 1) / 2
    nij = Counter(zip(x, y)); a = Counter(x); b = Counter(y)
    s = sum(comb(v) for v in nij.values()); sa = sum(comb(v) for v in a.values()); sb = sum(comb(v) for v in b.values())
    e = sa * sb / comb(n); m = (sa + sb) / 2
    return (s - e) / (m - e) if m != e else 1.0


labs = {}; summ = []
Zk = {k: feats(k) for k in (5, 6, 7)}
for k, Z in Zk.items():
    for K in range(2, 7):
        lab, d = kmeans(Z, K); labs[(k, K)] = lab
        sil = silhouette(Z, lab)
        summ.append(dict(k=k, K=K, silhouette=round(sil, 3), sizes=','.join(str(int((lab == j).sum())) for j in range(K))))
for K in range(2, 7):
    for r in summ:
        if r['K'] == K: r['ARI_vs_k6'] = round(ari(list(labs[(r['k'], K)]), list(labs[(6, K)])), 3)
# describe clusters at each K (k = 6) with independent evidence
old = np.array([m['old_lineage'] for m in meta]); sub = np.array([m['subgenome'] for m in meta])
L = np.array([int(m['length']) for m in meta]); F5 = np.array([s5.get(c, 0.0) for c in ids])
with open(f'{D}/55_cluster_summary.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['k', 'K', 'silhouette', 'ARI_vs_k6', 'cluster', 'n', 'Gar', 'Mac', 'old_A', 'old_B', 'old_none', 'median_len', 'median_5S_score', 'frac_5S_score_ge40'])
    for r in summ:
        lab = labs[(r['k'], r['K'])]
        for j in range(r['K']):
            m = lab == j
            w.writerow([r['k'], r['K'], r['silhouette'], r['ARI_vs_k6'], j, int(m.sum()), int((sub[m] == 'Gar').sum()), int((sub[m] == 'Mac').sum()),
                        int((old[m] == 'A').sum()), int((old[m] == 'B').sum()), int((old[m] == 'none').sum()), int(np.median(L[m])),
                        round(float(np.median(F5[m])), 1), round(float((F5[m] >= 40).mean()), 3)])
with open(f'{D}/55_cluster_assignments.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n'); keys = sorted(labs)
    w.writerow(['copy_id', 'subgenome', 'array_id', 'segment', 'family', 'old_lineage', 'length', 'score_5S'] + [f'k{k}_K{K}' for k, K in keys])
    for i, m in enumerate(meta): w.writerow([m['copy_id'], m['subgenome'], m['array_id'], m['segment'], m['family'], m['old_lineage'], m['length'], round(F5[i], 1)] + [int(labs[x][i]) for x in keys])
for r in summ: print(r)
