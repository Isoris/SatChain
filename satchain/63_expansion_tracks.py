#!/usr/bin/env python3
"""
63_expansion_tracks.py  (SatChain v4.7 draft; sequence-only tracks of recent expansion and unit haplotypes along one array)

Units (--mode):
  copy : every confident copy is a unit (monomeric arrays);
  hor  : consecutive blocks of --hor copies (from 62_hor_calls; --hor > 0 is required, 2026-11-01 D10); a block is formed only from
         copies that are ADJACENT in the genome (distance between neighbouring copies <= --adj-tol bp, default 20; before: by copy index,
         so a block could span a missing or junk copy); the phase (0..hor-1) is the one with the highest mean identity between
         neighbouring blocks;
  bed  : units from a BED (e.g. 61_units_<G>.bed), sequence taken from the genome in the array orientation (strand column); for 61 units
         (run + following interruption) the satellite / interruption split is in BED columns 7-8, so the TE share of a unit is known.
Tracks per unit, in genomic order:
  exact      : 1 if an identical unit sequence occurs elsewhere in the array;
  exact_win  : number of units with an exact match among this unit and the next --win - 1 (window of --win units); win_units = units
               actually in the window (< --win for the last --win - 1 units: truncated, 2026-11-01 D10; plot63 scales by win_units);
  hap_<K>    : unit haplotype, k-means on centred k-mer (k = --hap-k) frequency profiles, for every K in --hap-K (e.g. broad 2, fine 6);
               the silhouette of every K is written to 63_haps_<name>.tsv (K validation).
Exact-match expectation (D10): the exact-match rate is compared with a column-shuffle null: units of the modal length form a matrix,
every column (position) is permuted independently across units --null times (site frequencies kept, linkage destroyed); exact_rate_null
= mean exact-match rate of the shuffled matrices. A rate at the null means identical units are what independent sites already produce
(short units); only the excess speaks for recent copying.
Self-identity matrix: the array is cut into --bins bins; each bin is a bottom-s MinHash sketch (s = --sketch) of canonical --dot-k-mers;
  Jaccard between all bin pairs (63_dotplot_<name>.npy + bin coordinates).
Outputs (<dest>/): 63_units_<name>.tsv, 63_haps_<name>.tsv, 63_summary_<name>.tsv, 63_dotplot_<name>.npy, 63_dotplot_bins_<name>.tsv.
--sweep writes K x window to <repo>/benchmarks/63_sweep.tsv (or --sweep-out).
usage: 63_expansion_tracks.py --chrom C --start S --end E --name N --genome <fa> --dest <dir> --mode copy|hor|bed
       [--copies 03_monomer_copies.tsv ...] [--hor n] [--adj-tol 20] [--bed units.bed] [--win 10] [--hap-k 5] [--hap-K 2,6] [--bins 300]
       [--dot-k 21] [--sketch 300] [--null 20] [--sweep]
"""
import argparse, collections, csv, hashlib, os, random
import numpy as np
from Bio import Align
csv.field_size_limit(1 << 30)
RC = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(RC)[::-1]
GA = Align.PairwiseAligner(mode='global', match_score=1, mismatch_score=0, open_gap_score=-1, extend_gap_score=-1)


class Fasta:
    def __init__(self, path):
        self.fh = open(path, 'rb'); self.idx = {}
        for l in open(path + '.fai'):
            n, ln, off, lb, lby = l.split('\t')[:5]; self.idx[n] = (int(ln), int(off), int(lb), int(lby))
    def get(self, c, a, b):
        ln, off, lb, lby = self.idx[c]; a = max(0, a); b = min(ln, b)
        s0 = off + (a // lb) * lby + a % lb; s1 = off + (b // lb) * lby + b % lb
        self.fh.seek(s0); return self.fh.read(s1 - s0).replace(b'\n', b'').replace(b'\r', b'').decode().upper()


def h64(x): return int.from_bytes(hashlib.blake2b(x.encode(), digest_size=8).digest(), 'little')


def kprof(s, k):
    c = collections.Counter(s[i:i + k] for i in range(len(s) - k + 1) if 'N' not in s[i:i + k]); return c


def kmeans(X, K, rng, it=50):
    best = None
    for rep in range(5):
        C = X[rng.choice(len(X), K, replace=False)]
        for _ in range(it):
            lab = np.argmin(((X[:, None, :] - C[None]) ** 2).sum(-1), 1)
            C2 = np.array([X[lab == j].mean(0) if (lab == j).any() else C[j] for j in range(K)])
            if np.allclose(C, C2): break
            C = C2
        inertia = ((X - C[lab]) ** 2).sum()
        if best is None or inertia < best[0]: best = (inertia, lab)
    # relabel by size (0 = largest)
    lab = best[1]; order = [j for j, _ in collections.Counter(lab).most_common()]
    return np.array([order.index(j) for j in lab])


def silhouette(X, lab):
    """mean silhouette over points (euclidean); nan when one cluster"""
    if len(set(lab.tolist())) < 2: return float('nan')
    D = np.sqrt(((X[:, None, :] - X[None]) ** 2).sum(-1)); s = []
    for i in range(len(X)):
        own = (lab == lab[i]).copy(); own[i] = False
        a_ = D[i, own].mean() if own.any() else 0.0
        b_ = min(D[i, lab == g].mean() for g in set(lab.tolist()) if g != lab[i])
        s.append((b_ - a_) / max(a_, b_, 1e-12))
    return float(np.mean(s))


def hor_blocks(cp, n, adj_tol, seed):
    """blocks of n genomically adjacent copies (gap between neighbours <= adj_tol), phase chosen by neighbour identity (D10)"""
    adj = [cp[i + 1][0] - cp[i][1] <= adj_tol for i in range(len(cp) - 1)]
    best = None
    for ph in range(n):
        us = []
        i = ph
        while i + n <= len(cp):
            if all(adj[i:i + n - 1]): us.append((cp[i][0], cp[i + n - 1][1], ''.join(x[2] for x in cp[i:i + n]))); i += n
            else: i += 1                                           # skip past the break, keep the phase
        if len(us) < 2: continue
        smp = random.Random(seed).sample(range(len(us) - 1), min(100, len(us) - 1))
        sc = np.mean([GA.score(us[i][2], us[i + 1][2]) / max(len(us[i][2]), len(us[i + 1][2])) for i in smp])
        if best is None or sc > best[0]: best = (sc, ph, us)
    return best


def exact_null(seqs, reps, rng):
    """column-shuffle null of the exact-match rate on units of the modal length"""
    L = collections.Counter(len(s) for s in seqs).most_common(1)[0][0]; M = np.array([list(s) for s in seqs if len(s) == L])
    if len(M) < 2: return float('nan'), 0
    out = []
    for _ in range(reps):
        S = M.copy()
        for j in range(L): S[:, j] = S[rng.permutation(len(S)), j]
        rows = [''.join(r) for r in S]; cnt = collections.Counter(rows); out.append(np.mean([cnt[r] > 1 for r in rows]))
    return float(np.mean(out)), len(M)


def main():
    ap = argparse.ArgumentParser()
    for x in ('--chrom', '--name', '--genome', '--dest', '--mode'): ap.add_argument(x, required=True)
    ap.add_argument('--start', type=int, required=True); ap.add_argument('--end', type=int, required=True)
    ap.add_argument('--copies', nargs='*', default=[]); ap.add_argument('--hor', type=int, default=0); ap.add_argument('--adj-tol', type=int, default=20); ap.add_argument('--bed', default='')
    ap.add_argument('--win', type=int, default=10); ap.add_argument('--hap-k', type=int, default=5); ap.add_argument('--hap-K', default='2,6')
    ap.add_argument('--bins', type=int, default=300); ap.add_argument('--dot-k', type=int, default=21); ap.add_argument('--sketch', type=int, default=300)
    ap.add_argument('--min-identity', type=float, default=60); ap.add_argument('--null', type=int, default=20); ap.add_argument('--seed', type=int, default=63)
    ap.add_argument('--sweep', action='store_true'); ap.add_argument('--sweep-out', default='')
    a = ap.parse_args(); rng = np.random.default_rng(a.seed); os.makedirs(a.dest, exist_ok=True); fa = Fasta(a.genome)
    if a.mode == 'hor' and a.hor < 2: raise SystemExit('--mode hor needs --hor n >= 2 (the 62 call); got %d' % a.hor)
    units = []                                   # (start, end, seq)
    if a.mode in ('copy', 'hor'):
        cp = []
        for p in a.copies:
            with open(p, newline='') as fh:
                rd = csv.reader(fh, delimiter='\t'); h = next(rd); ix = {k: i for i, k in enumerate(h)}
                for row in rd:
                    if row[ix['chrom']] != a.chrom: continue
                    s, e = int(row[ix['genomic_start_0based']]), int(row[ix['genomic_end']])
                    if s >= a.start and e <= a.end and float(row[ix['best_alignment_identity']] or 0) >= a.min_identity:
                        cp.append((s, e, row[ix['sequence']].upper()))
        cp.sort()
        if a.mode == 'copy': units = cp
        else:
            best = hor_blocks(cp, a.hor, a.adj_tol, a.seed)
            if best is None: raise SystemExit('no two adjacent HOR blocks of %d copies (adj-tol %d bp)' % (a.hor, a.adj_tol))
            units = best[2]; print(f'HOR phase {best[1]} (mean neighbour identity {best[0]:.3f}); {len(units)} blocks from {len(cp)} copies')
    else:
        for l in open(a.bed):
            c, s, e, nm, _, st = l.rstrip('\n').split('\t')[:6]; s, e = int(s), int(e)
            if c == a.chrom and s >= a.start and e <= a.end:
                q = fa.get(c, s, e); units.append((s, e, q if st == '+' else rc(q)))
        units.sort()
    if len(units) < 10: raise SystemExit(f'only {len(units)} units')
    seqs = [u[2] for u in units]; cnt = collections.Counter(seqs)
    exact = np.array([cnt[s] > 1 for s in seqs], int)
    win = np.array([exact[i:i + a.win].sum() for i in range(len(exact))]); win_n = np.array([min(a.win, len(exact) - i) for i in range(len(exact))])
    null_rate, null_n = exact_null(seqs, a.null, rng)
    K = [int(x) for x in a.hap_K.split(',')]
    allk = sorted({x for s in seqs for x in kprof(s, a.hap_k)})
    ix = {x: i for i, x in enumerate(allk)}; X = np.zeros((len(seqs), len(allk)))
    for i, s in enumerate(seqs):
        for x, v in kprof(s, a.hap_k).items(): X[i, ix[x]] = v
    X = X / X.sum(1, keepdims=True); X = X - X.mean(0)
    U, S_, Vt = np.linalg.svd(X, full_matrices=False); Xp = U[:, :10] * S_[:10]
    haps = {k: kmeans(Xp, k, rng) for k in K}; sil = {k: silhouette(Xp, haps[k]) for k in K}
    with open(os.path.join(a.dest, f'63_units_{a.name}.tsv'), 'w') as fh:
        fh.write('unit\tstart\tend\tlength\texact\texact_win\twin_units\t' + '\t'.join(f'hap_{k}' for k in K) + '\n')
        for i, (s, e, q) in enumerate(units):
            fh.write(f'{i}\t{s}\t{e}\t{len(q)}\t{exact[i]}\t{win[i]}\t{win_n[i]}\t' + '\t'.join(str(haps[k][i]) for k in K) + '\n')
    with open(os.path.join(a.dest, f'63_haps_{a.name}.tsv'), 'w') as fh:
        fh.write('K\tsilhouette\tsizes\n')
        for k in K: fh.write(f'{k}\t{sil[k]:.4f}\t{",".join(str(v) for _, v in sorted(collections.Counter(haps[k].tolist()).items()))}\n')
    with open(os.path.join(a.dest, f'63_summary_{a.name}.tsv'), 'w') as fh:
        fh.write('name\tmode\tunits\texact_rate\texact_rate_null\texact_excess\tnull_units_modal_length\twin\ttruncated_windows\n')
        fh.write(f'{a.name}\t{a.mode}\t{len(units)}\t{exact.mean():.4f}\t{null_rate:.4f}\t{exact.mean() - null_rate:.4f}\t{null_n}\t{a.win}\t{int((win_n < a.win).sum())}\n')
    # self-identity matrix on MinHash sketches
    seq = fa.get(a.chrom, a.start, a.end); B = max(1000, (a.end - a.start) // a.bins); edges = list(range(a.start, a.end, B))
    sk = []
    for b0 in edges:
        s = seq[b0 - a.start:b0 - a.start + B]; r = rc(s); n = len(s); hs = set()
        for i in range(n - a.dot_k + 1):
            x = s[i:i + a.dot_k]
            if 'N' in x: continue
            y = r[n - i - a.dot_k:n - i]; hs.add(h64(min(x, y)))
        sk.append(np.array(sorted(hs)[:a.sketch], dtype=np.uint64))
    M = np.zeros((len(sk), len(sk)), np.float32)
    for i in range(len(sk)):
        for j in range(i, len(sk)):
            if len(sk[i]) == 0 or len(sk[j]) == 0: continue
            u = np.union1d(sk[i], sk[j])[:a.sketch]; inter = np.intersect1d(np.intersect1d(sk[i], sk[j]), u)
            M[i, j] = M[j, i] = len(inter) / len(u)
    np.save(os.path.join(a.dest, f'63_dotplot_{a.name}.npy'), M)
    with open(os.path.join(a.dest, f'63_dotplot_bins_{a.name}.tsv'), 'w') as fh:
        fh.write('bin\tstart\tend\n'); [fh.write(f'{i}\t{b0}\t{min(a.end, b0 + B)}\n') for i, b0 in enumerate(edges)]
    print(f'{a.name}: {len(units)} units, exact-match units {exact.sum()} ({exact.mean():.1%}; column-shuffle null {null_rate:.1%}), '
          f'silhouette {", ".join(f"K{k} {sil[k]:.2f}" for k in K)}, bins {len(edges)} x {B} bp')
    if a.sweep:                                  # K x window sweep (D09)
        SW = []
        for k in sorted(set(K) | {2, 3, 4, 6, 8}):
            if k < len(seqs): SW.append(dict(step='63', name=a.name, K=k, silhouette=round(silhouette(Xp, kmeans(Xp, k, rng)), 4)))
        for wn in (5, 10, 20):
            wv = np.array([exact[i:i + wn].sum() / min(wn, len(exact) - i) for i in range(len(exact))])
            SW.append(dict(step='63', name=a.name, win=wn, median_exact_share=round(float(np.median(wv)), 4), p90_exact_share=round(float(np.percentile(wv, 90)), 4)))
        out = a.sweep_out or os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'benchmarks', '63_sweep.tsv'); os.makedirs(os.path.dirname(out), exist_ok=True)
        keys = list(dict.fromkeys(k for r in SW for k in r))
        with open(out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=keys, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(SW)
        print('sweep ->', out)


if __name__ == '__main__':
    main()
