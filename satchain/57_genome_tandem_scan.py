#!/usr/bin/env python3
"""
57_genome_tandem_scan.py  (v4.6; genome-wide, label-free tandem signal; candidate generator replacing the fixed step-02 threshold)

Every chromosome is tiled with windows at two scales:
  small  --small bp windows, step half a window, lags 2..--small-lag   (micro-, mini- and satellites up to ~1 kb units)
  large  --large bp windows, step half a window, lags --small-lag..--large-lag   (long units, 1-5 kb)
For each window: the k-mer (k = 10) next-recurrence distance histogram minus the same histogram of a base-shuffled copy of the window
(the step-48 excess), smoothed +/- 6 % of the lag; reported: the maximum excess, its lag, the excess at 2x that lag (harmonic support),
the fraction of the window's k-mers that recur within the lag range, the N fraction.
No threshold is applied here: calibration against the clean null and the tiers are a separate step.
Output: one TSV per chromosome group: chrom, scale, start, end, n_frac, max_excess, best_lag, excess_at_2x, recur_frac
usage: 57_genome_tandem_scan.py <genome fasta> <tag> <out tsv> <chrom1,chrom2,...|all> [--small 2000 --large 10000]
"""
import argparse, sys
import numpy as np
CODE = np.full(256, 4, np.int64)
for i, c in enumerate(b'ACGT'): CODE[c] = i; CODE[c + 32] = i
K = 10


def read(fa, want):
    seqs, name, buf = {}, None, []
    for l in open(fa, 'rb'):
        if l.startswith(b'>'):
            if name and (want is None or name in want): seqs[name] = b''.join(buf)
            name = l[1:].split()[0].decode(); buf = []
        elif want is None or name in want: buf.append(l.rstrip())
    if name and (want is None or name in want): seqs[name] = b''.join(buf)
    return seqs


def kcode(c):
    n = len(c) - K + 1
    if n <= 1: return None, None
    h = np.zeros(n, dtype=np.int64); bad = np.zeros(n, dtype=bool)
    for j in range(K): h = h * 4 + np.minimum(c[j:j + n], 3); bad |= c[j:j + n] >= 4
    return h, bad


def hist(c, lo, hi):
    h, bad = kcode(c)
    if h is None: return np.zeros(hi + 1), 1, 0
    pos = np.flatnonzero(~bad)
    if len(pos) < 2: return np.zeros(hi + 1), 1, 0
    hv = h[pos]; o = np.lexsort((pos, hv)); hs, ps = hv[o], pos[o]
    d = (ps[1:] - ps[:-1])[hs[1:] == hs[:-1]]; rec = ((d >= lo) & (d <= hi)).sum()
    d = d[d <= hi]
    return np.bincount(d, minlength=hi + 1)[:hi + 1].astype(float), len(pos), rec


def smooth(h, nv, rel=0.06):
    cs = np.concatenate([[0.0], np.cumsum(h)]); L = np.arange(len(h)); w = np.maximum(1, np.round(rel * L)).astype(int)
    return (cs[np.clip(L + w + 1, 0, len(h))] - cs[np.clip(L - w, 0, len(h))]) / max(1, nv)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('fasta'); ap.add_argument('tag'); ap.add_argument('out'); ap.add_argument('chroms')
    ap.add_argument('--small', type=int, default=2000); ap.add_argument('--large', type=int, default=10000)
    ap.add_argument('--small-lag', type=int, default=1000); ap.add_argument('--large-lag', type=int, default=5000)
    a = ap.parse_args(); rng = np.random.default_rng(57)
    seqs = read(a.fasta, None if a.chroms == 'all' else set(a.chroms.split(',')))
    with open(a.out, 'w') as fh:
        fh.write('chrom\tscale\tstart\tend\tn_frac\tmax_excess\tbest_lag\texcess_at_2x\trecur_frac\n')
        for chrom, s in seqs.items():
            c = CODE[np.frombuffer(s, dtype=np.uint8)]
            for scale, W, lo, hi in (('small', a.small, 2, a.small_lag), ('large', a.large, a.small_lag, a.large_lag)):
                hi_ = min(hi, W // 2)
                for st in range(0, max(1, len(c) - W + 1), W // 2):
                    w = c[st:st + W]; nf = float((w >= 4).mean())
                    if nf > 0.5: fh.write(f'{chrom}\t{scale}\t{st}\t{st + len(w)}\t{nf:.3f}\tNA\tNA\tNA\tNA\n'); continue
                    h, nv, rec = hist(w, lo, hi_); h0, nv0, _ = hist(rng.permutation(w), lo, hi_)
                    ex = smooth(h, nv) - smooth(h0, nv0); ex[:max(2, lo)] = -1
                    j = int(np.argmax(ex)); e2 = float(ex[2 * j]) if 2 * j < len(ex) else float('nan')
                    fh.write(f'{chrom}\t{scale}\t{st}\t{st + len(w)}\t{nf:.3f}\t{ex[j]:.4f}\t{j}\t{e2:.4f}\t{rec / max(1, nv):.4f}\n')
            print(a.tag, chrom, 'done', flush=True)


if __name__ == '__main__':
    main()
