#!/usr/bin/env python3
"""val58: independent support of the NEW scan-57 candidates (outside every step-02 array), per unit class:
  TRF support    a TRF call overlapping the candidate with period within 10 % of the candidate lag (or a 2x / 3x multiple / divisor)
                 and >= 3 copies
  composition    on a seeded sample of --n candidates per class: k-mer recurrence excess at the candidate lag against a
                 dinucleotide-block shuffle (non-overlapping 2-mers permuted) instead of a base shuffle; supported if > 0.02
usage: val58_new.py <tag> <cal58 dir> <tier> <func genome folder> <pipeline out> <genome fasta> [n]"""
import sys, gzip, bisect, collections, csv, random, re
import numpy as np
TAG, CAL, TIER, FUNC, OUT, FA = sys.argv[1:7]; N = int(sys.argv[7]) if len(sys.argv) > 7 else 200
rng = random.Random(58); NR = np.random.default_rng(58)
a02 = collections.defaultdict(list)
for r in csv.DictReader(open(f'{OUT}/fClaHyb_{TAG}_satellites/02_satellite_arrays.tsv'), delimiter='\t'): a02[r['chrom']].append((int(r['start_0based']), int(r['end'])))
for v in a02.values(): v.sort()
def ov(d, c, s, e):
    v = d.get(c, []); i = bisect.bisect_left(v, (s - 3000000, 0)); t = 0
    for x in v[i:]:
        if x[0] >= e: break
        t += max(0, min(e, x[1]) - max(s, x[0]))
    return t
trf = collections.defaultdict(list)
for l in gzip.open(f'{FUNC}/_TR/TRF/fClaHyb_{TAG}.TRF.bed.gz', 'rt'):
    f = l.rstrip('\n').split('\t'); m = re.match(r'(?:(\d+)bp|([ACGTN]+))x([\d.]+)', f[3])
    if not m: continue
    per = int(m.group(1)) if m.group(1) else len(m.group(2)); trf[f[0]].append((int(f[1]), int(f[2]), per, float(m.group(3))))
for v in trf.values(): v.sort()
def trf_support(c, s, e, lag):
    v = trf.get(c, []); i = bisect.bisect_left(v, (s - 200000, 0, 0, 0))
    for x in v[i:]:
        if x[0] >= e: break
        if x[1] <= s or x[3] < 3: continue
        r = max(x[2], lag) / min(x[2], lag); k = round(r)
        if k <= 3 and abs(r - k) <= 0.1 * k: return True
    return False
fai = {l.split()[0]: list(map(int, l.split()[1:5])) for l in open(FA + '.fai')}
fh = open(FA, 'rb')
def fetch(c, s, e):
    L, off, lb, lB = fai[c]; fh.seek(off + (s // lb) * lB + s % lb); return fh.read(int((e - s) * lB / lb) + 2 * lB).replace(b'\n', b'')[:e - s]
CODE = np.full(256, 4, np.int64)
for i, ch in enumerate(b'ACGT'): CODE[ch] = i; CODE[ch + 32] = i
def rec_at(c, L, k=10):
    n = len(c) - k + 1
    if n < 2: return 0.0
    h = np.zeros(n, dtype=np.int64); bad = np.zeros(n, dtype=bool)
    for j in range(k): h = h * 4 + np.minimum(c[j:j + n], 3); bad |= c[j:j + n] >= 4
    pos = np.flatnonzero(~bad); hv = h[pos]; o = np.lexsort((pos, hv)); hs, ps = hv[o], pos[o]
    d = (ps[1:] - ps[:-1])[hs[1:] == hs[:-1]]
    return float(((d >= 0.94 * L) & (d <= 1.06 * L)).sum()) / max(1, len(pos))
def chain_at(c, L, k=10, tol=0.06):
    """fraction of k-mer positions p whose k-mer recurs at ~p+L AND ~p+2L (>= 3 tandem copies); vs a dinucleotide shuffle"""
    n = len(c) - k + 1
    if n < 2: return 0.0
    h = np.zeros(n, dtype=np.int64); bad = np.zeros(n, dtype=bool)
    for j in range(k): h = h * 4 + np.minimum(c[j:j + n], 3); bad |= c[j:j + n] >= 4
    pos = np.flatnonzero(~bad); hv = h[pos]; o = np.lexsort((pos, hv)); hs, ps = hv[o], pos[o]
    same = hs[1:] == hs[:-1]; d = ps[1:] - ps[:-1]
    ok1 = same & (d >= (1 - tol) * L) & (d <= (1 + tol) * L)
    two = ok1[:-1] & ok1[1:]                       # consecutive occurrences of one k-mer, both steps ~L
    return float(two.sum()) / max(1, len(pos))


def dinuc_shuffle(c):
    n = len(c) // 2 * 2; b = c[:n].reshape(-1, 2); return np.concatenate([b[NR.permutation(len(b))].ravel(), c[n:]])
classes = collections.defaultdict(list)
for l in open(f'{CAL}/58_candidates_p{TIER}.bed'):
    c, s, e, nm = l.split('\t')[:4]; s, e = int(s), int(e); lag = int(nm.split('|')[1][3:])
    out = (e - s) - ov(a02, c, s, e)
    if out <= 0 or lag < 7: continue
    cls = '7-39' if lag <= 39 else '40-60' if lag <= 60 else '61-1000' if lag <= 1000 else '>1000'
    classes[cls].append((c, s, e, lag, out))
print(f'{TAG} tier p{TIER}: class\tn\tMb\tTRF-supported (n / bp share)\tdinucleotide-shuffle supported (sample)')
for cls in ('7-39', '40-60', '61-1000', '>1000'):
    v = classes[cls]
    if not v: continue
    t = [x for x in v if trf_support(x[0], x[1], x[2], x[3])]
    smp = rng.sample(v, min(N, len(v))); ok = 0; ch = 0; both = 0
    for c, s, e, lag, _ in smp:
        sq = CODE[np.frombuffer(fetch(c, s, min(e, s + 20000)), dtype=np.uint8)]
        o1 = (rec_at(sq, lag) - rec_at(dinuc_shuffle(sq), lag)) > 0.02; o2 = (chain_at(sq, lag) - chain_at(dinuc_shuffle(sq), lag)) > 0.01
        ok += o1; ch += o2; both += o1 and o2 and trf_support(c, s, e, lag)
    print(f'{TAG}\t{cls}\t{len(v)}\t{sum(x[4] for x in v) / 1e6:.2f}\t{len(t) / len(v):.2f} / {sum(x[4] for x in t) / max(1, sum(x[4] for x in v)):.2f}\t{ok / len(smp):.2f}	3-copy chain {ch / len(smp):.2f}	chain+TRF {both / len(smp):.2f} (n={len(smp)})')
