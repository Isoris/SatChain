#!/usr/bin/env python3
"""benchmark_band_test.py (2026-11-01; ships the band test whose numbers the README quoted but whose script was not in the repo, audit HB01).
Question: are the v4.5 scan windows just below the old lag-identity threshold (0.34) tandem on the sequence? Bands of lag_identity_max from
02_scan_windows.tsv.gz; --n windows per band (seeded); for each window the k-mer (k = 10) next-recurrence excess over a base-shuffled copy
(as step 57/58), its best lag, and the tandem fraction at excess >= 0.01 / 0.02 / 0.05.
Null = --n random 10-kb genome windows OUTSIDE every step-02 interval, every catalogue run and every TRF call with unit >= 7 bp (the
2026-10 version drew the null genome-wide, so 46-61 % of 'random' windows exceeded 0.01). The summary reports the null exceedance at the
same cuts, so a band's tandem share can be read against it. Writes <out>.tsv (per window) and <out>_summary.tsv.
usage: benchmark_band_test.py <pipeline out of one subgenome> <genome fasta> <tag> <catalogue tsv|-> <TRF bed.gz|-> <out prefix> [n=300]"""
import gzip, sys, random, collections, csv, bisect, re
import numpy as np
OUT, FA, TAG, CAT, TRF, DEST = sys.argv[1:7]; N = int(sys.argv[7]) if len(sys.argv) > 7 else 300
rng = random.Random(2); NR = np.random.default_rng(2)
CODE = np.full(256, 4, np.int64)
for i, c in enumerate(b'ACGT'): CODE[c] = i; CODE[c + 32] = i
BANDS = [(0, 0.25), (0.25, 0.30), (0.30, 0.32), (0.32, 0.34), (0.34, 0.36), (0.36, 0.40), (0.40, 1.01)]
win = collections.defaultdict(list)
with gzip.open(f'{OUT}/02_scan_windows.tsv.gz', 'rt') as fh:
    next(fh)
    for l in fh:
        c, s, e, p, mx, *_ = l.split('\t')
        if mx == 'NA' or p == 'NA': continue
        mx = float(mx)
        for b in BANDS:
            if b[0] <= mx < b[1]: win[b].append((c, int(s), int(e), int(p), mx)); break
mask = collections.defaultdict(list)                                  # excluded from the null
for r in csv.DictReader(open(f'{OUT}/02_satellite_arrays.tsv'), delimiter='\t'): mask[r['chrom']].append((int(r['start_0based']), int(r['end'])))
if CAT != '-':
    for r in csv.DictReader(open(CAT), delimiter='\t'): mask[r['chrom']].append((int(r['start0']), int(r['end'])))
if TRF != '-':
    for l in gzip.open(TRF, 'rt'):
        f = l.rstrip('\n').split('\t'); m = re.match(r'(?:(\d+)bp|([ACGTN]+))x([\d.]+)', f[3])
        if m and (int(m.group(1)) if m.group(1) else len(m.group(2))) >= 7: mask[f[0]].append((int(f[1]), int(f[2])))
for c in mask:
    v = sorted(mask[c]); m = []
    for s, e in v:
        if m and s <= m[-1][1]: m[-1] = (m[-1][0], max(m[-1][1], e))
        else: m.append((s, e))
    mask[c] = m
starts = {c: [x[0] for x in v] for c, v in mask.items()}
def masked(c, s, e):
    v = mask.get(c, []); i = bisect.bisect_left(starts.get(c, []), e)
    return any(x[1] > s for x in v[max(0, i - 3):i])
seqs = {}; name = None; buf = []
for l in open(FA, 'rb'):
    if l.startswith(b'>'):
        if name: seqs[name] = b''.join(buf)
        name = l[1:].split()[0].decode(); buf = []
    else: buf.append(l.rstrip())
seqs[name] = b''.join(buf)


def excess(s, lmax=2000, k=10):
    def rec(b):
        c = CODE[np.frombuffer(b, dtype=np.uint8)]; n = len(c) - k + 1
        h = np.zeros(n, dtype=np.int64); bad = np.zeros(n, dtype=bool)
        for j in range(k): h = h * 4 + np.minimum(c[j:j + n], 3); bad |= c[j:j + n] >= 4
        pos = np.flatnonzero(~bad); hv = h[pos]; o = np.lexsort((pos, hv)); hs, ps = hv[o], pos[o]
        d = (ps[1:] - ps[:-1])[hs[1:] == hs[:-1]]; d = d[d <= lmax]
        return np.bincount(d, minlength=lmax + 1)[:lmax + 1].astype(float), max(1, len(pos))
    def sm(h, nv, rel=0.06):
        cs = np.concatenate([[0.0], np.cumsum(h)]); L = np.arange(len(h)); w = np.maximum(1, np.round(rel * L)).astype(int)
        return (cs[np.clip(L + w + 1, 0, len(h))] - cs[np.clip(L - w, 0, len(h))]) / nv
    h, n = rec(s); arr = np.frombuffer(s, dtype=np.uint8); ok = CODE[arr] < 4; sh = arr.copy(); sh[ok] = NR.permutation(arr[ok]); h0, n0 = rec(bytes(sh))   # permute non-N only (CB15)
    ex = sm(h, n) - sm(h0, n0); ex[:7] = -1                                                                                                  # lags >= 7 as in 58
    j = int(np.argmax(ex)); return float(ex[j]), j


rows = []
for b in BANDS:
    for c, s, e, p, mx in rng.sample(win[b], min(N, len(win[b]))):
        sq = seqs[c][s:e]
        if sq.upper().count(b'N') > 0.1 * len(sq): continue
        ex, L = excess(sq); rows.append((f'{b[0]:.2f}-{min(b[1], 1):.2f}', len(win[b]), c, s, e, p, mx, ex, L))
chroms = [c for c in seqs]; lens = [len(seqs[c]) for c in chroms]; got = 0; tried = 0
while got < N and tried < 100 * N:
    tried += 1; c = rng.choices(chroms, weights=lens)[0]; s = rng.randrange(0, len(seqs[c]) - 10000); sq = seqs[c][s:s + 10000]
    if sq.upper().count(b'N') > 0.1 * len(sq) or masked(c, s, s + 10000): continue
    ex, L = excess(sq); rows.append(('random_outside_arrays', 0, c, s, s + 10000, 0, float('nan'), ex, L)); got += 1
with open(DEST + '.tsv', 'w') as fh:
    fh.write('band\tband_windows_total\tchrom\tstart\tend\tscan_period\tlag_identity_max\tmax_excess\tbest_lag\n')
    for r in rows: fh.write('\t'.join(map(str, r)) + '\n')
by = collections.defaultdict(list)
for r in rows: by[r[0]].append(r)
null = np.array([x[7] for x in by['random_outside_arrays']]); p95 = float(np.quantile(null, 0.95)) if len(null) else float('nan')
with open(DEST + '_summary.tsv', 'w') as fh:
    fh.write('tag\tband\twindows_in_band\tsampled\ttandem_at_0.01\ttandem_at_0.02\ttandem_at_0.05\ttandem_above_null_p95\tmedian_excess\tnull_p95\n')
    for k_, v in by.items():
        e = np.array([x[7] for x in v])
        fh.write(f'{TAG}\t{k_}\t{v[0][1]}\t{len(v)}\t{(e >= 0.01).mean():.3f}\t{(e >= 0.02).mean():.3f}\t{(e >= 0.05).mean():.3f}\t{(e > p95).mean():.3f}\t{np.median(e):.4f}\t{p95:.4f}\n')
print(open(DEST + '_summary.tsv').read())
