#!/usr/bin/env python3
"""val58: independent support of the NEW scan-57 candidates (outside every step-02 array), per unit class:
  TRF support    2026-11-01 (CB03): compatible TRF calls = period >= --min-per (7 bp), >= 3 copies, period within 10 % of the
                 candidate lag or of a 2x / 3x multiple / divisor; support = union of the compatible calls covers >= --trf-cover
                 (default 0.5) of the candidate. The pre-2026-11-01 flag (any compatible call, any period, >= 1 bp overlap) is
                 still returned side by side as 'trf_asis'.
  composition    on a seeded sample of --n candidates per class: k-mer recurrence excess at the candidate lag against a
                 dinucleotide-block shuffle (non-overlapping 2-mers permuted) instead of a base shuffle; supported if > 0.02
k-mers are forward-strand 10-mers (not canonical; copies of one array lie on one strand).
2026-11-01 (CA10 / CB16): the helpers are plain functions (58b and 58c import this module instead of exec()-ing its source) and
every annotation file is resolved with satpaths.func_path (raises when missing).
usage: val58_new.py <tag> <cal58 dir> <tier> <func genome folder> <pipeline out> <genome fasta> [n] [--trf-cover 0.5] [--min-per 7]"""
import sys, os, gzip, bisect, collections, csv, random, re
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from satpaths import func_path
CODE = np.full(256, 4, np.int64)
for i, ch in enumerate(b'ACGT'): CODE[ch] = i; CODE[ch + 32] = i
CLASSES = ('7-39', '40-60', '61-1000', '>1000')


def unit_class(lag): return '7-39' if lag <= 39 else '40-60' if lag <= 60 else '61-1000' if lag <= 1000 else '>1000'


def load_a02(out, tag):
    """step-02 arrays per chromosome, sorted (start, end)"""
    a02 = collections.defaultdict(list)
    for r in csv.DictReader(open(f'{out}/fClaHyb_{tag}_satellites/02_satellite_arrays.tsv'), delimiter='\t'):
        a02[r['chrom']].append((int(r['start_0based']), int(r['end'])))
    for v in a02.values(): v.sort()
    return a02


def ov(d, c, s, e):
    """bp of [s, e) covered by the (possibly overlapping) intervals of d[c]; summed, not a union"""
    v = d.get(c, []); i = bisect.bisect_left(v, (s - 3000000, 0)); t = 0
    for x in v[i:]:
        if x[0] >= e: break
        t += max(0, min(e, x[1]) - max(s, x[0]))
    return t


def load_trf(func, tag):
    """TRF calls per chromosome: (start, end, period, copies), sorted; from the full TRF bed (name 'NNbpxC' or 'MOTIFxC')"""
    trf = collections.defaultdict(list)
    for l in gzip.open(func_path(func, tag, '_TR', f'TRF/fClaHyb_{tag}.TRF.bed.gz'), 'rt'):
        f = l.rstrip('\n').split('\t'); m = re.match(r'(?:(\d+)bp|([ACGTN]+))x([\d.]+)', f[3])
        if not m: continue
        per = int(m.group(1)) if m.group(1) else len(m.group(2)); trf[f[0]].append((int(f[1]), int(f[2]), per, float(m.group(3))))
    for v in trf.values(): v.sort()
    return trf


def compatible_period(per, lag):
    r = max(per, lag) / min(per, lag); k = round(r)
    return k <= 3 and abs(r - k) <= 0.1 * k


def trf_support(trf, c, s, e, lag, min_per=7):
    """-> (trf_asis, cover): trf_asis = pre-2026-11-01 flag (any compatible call with >= 3 copies overlapping by >= 1 bp, any period);
    cover = fraction of [s, e) covered by the union of compatible calls with period >= min_per and >= 3 copies (CB03)."""
    v = trf.get(c, []); i = bisect.bisect_left(v, (s - 200000, 0, 0, 0)); asis = False; pieces = []
    for x in v[i:]:
        if x[0] >= e: break
        if x[1] <= s or x[3] < 3 or not compatible_period(x[2], lag): continue
        asis = True
        if x[2] >= min_per: pieces.append((max(s, x[0]), min(e, x[1])))
    pieces.sort(); cov = 0; cur = None
    for ps, pe in pieces:
        if cur and ps <= cur[1]: cur[1] = max(cur[1], pe)
        else:
            if cur: cov += cur[1] - cur[0]
            cur = [ps, pe]
    if cur: cov += cur[1] - cur[0]
    return asis, cov / max(1, e - s)


def make_fetch(fa):
    """fetch(c, s, e) -> bytes of the genome slice [s, e) through the .fai (line width from the index; LF genome)"""
    fai = {l.split()[0]: list(map(int, l.split()[1:5])) for l in open(fa + '.fai')}
    fh = open(fa, 'rb')
    def fetch(c, s, e):
        L, off, lb, lB = fai[c]; fh.seek(off + (s // lb) * lB + s % lb)
        return fh.read(int((e - s) * lB / lb) + 2 * lB).replace(b'\n', b'').replace(b'\r', b'')[:e - s]
    return fetch


def _kmers(c, k):
    n = len(c) - k + 1
    h = np.zeros(n, dtype=np.int64); bad = np.zeros(n, dtype=bool)
    for j in range(k): h = h * 4 + np.minimum(c[j:j + n], 3); bad |= c[j:j + n] >= 4
    pos = np.flatnonzero(~bad); hv = h[pos]; o = np.lexsort((pos, hv)); return hv[o], pos[o], len(pos)


def rec_at(c, L, k=10, tol=0.06):
    """fraction of k-mer positions whose k-mer next recurs at ~L (within tol)"""
    if len(c) - k + 1 < 2: return 0.0
    hs, ps, npos = _kmers(c, k)
    d = (ps[1:] - ps[:-1])[hs[1:] == hs[:-1]]
    return float(((d >= (1 - tol) * L) & (d <= (1 + tol) * L)).sum()) / max(1, npos)


def chain_at(c, L, k=10, tol=0.06, count=False):
    """fraction of k-mer positions p whose k-mer recurs at ~p+L AND ~p+2L (>= 3 tandem copies); vs a dinucleotide shuffle.
    count=True returns (fraction, absolute number of chained positions, number of k-mer positions)  (2026-11-01, CB11)"""
    if len(c) - k + 1 < 2: return (0.0, 0, 0) if count else 0.0
    hs, ps, npos = _kmers(c, k)
    same = hs[1:] == hs[:-1]; d = ps[1:] - ps[:-1]
    ok1 = same & (d >= (1 - tol) * L) & (d <= (1 + tol) * L)
    two = ok1[:-1] & ok1[1:]                       # consecutive occurrences of one k-mer, both steps ~L
    n2 = int(two.sum())
    return (n2 / max(1, npos), n2, npos) if count else n2 / max(1, npos)


def dinuc_shuffle(c, nr):
    n = len(c) // 2 * 2; b = c[:n].reshape(-1, 2); return np.concatenate([b[nr.permutation(len(b))].ravel(), c[n:]])


def load_classes(cal, tier, a02):
    """58 candidates of the tier outside step-02 arrays with lag >= 7, by unit class: (chrom, s, e, lag, out_bp)"""
    classes = collections.defaultdict(list)
    for l in open(f'{cal}/58_candidates_p{tier}.bed'):
        c, s, e, nm = l.split('\t')[:4]; s, e = int(s), int(e); lag = int(nm.split('|')[1][3:])
        out = (e - s) - ov(a02, c, s, e)
        if out <= 0 or lag < 7: continue
        classes[unit_class(lag)].append((c, s, e, lag, out))
    return classes


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    for x in ('tag', 'cal', 'tier', 'func', 'out', 'fa'): ap.add_argument(x)
    ap.add_argument('n', nargs='?', type=int, default=200)
    ap.add_argument('--trf-cover', type=float, default=0.5); ap.add_argument('--min-per', type=int, default=7)
    A = ap.parse_args(); TAG, CAL, TIER, FUNC, OUT, FA, N = A.tag, A.cal, A.tier, A.func, A.out, A.fa, A.n
    TRF_COVER = A.trf_cover; MIN_PER = A.min_per
    rng = random.Random(58); NR = np.random.default_rng(58)
    a02 = load_a02(OUT, TAG); trf = load_trf(FUNC, TAG); fetch = make_fetch(FA); classes = load_classes(CAL, TIER, a02)
    print(f'{TAG} tier p{TIER}: class\tn\tMb\tTRF as coded (n / bp share)\tTRF cover >= {TRF_COVER} (n / bp share)\tmedian cover'
          f'\tdinucleotide-shuffle supported (sample)')
    for cls in CLASSES:
        v = classes[cls]
        if not v: continue
        sup = [trf_support(trf, x[0], x[1], x[2], x[3], MIN_PER) for x in v]
        t = [x for x, (a, _) in zip(v, sup) if a]; tc = [x for x, (_, cv) in zip(v, sup) if cv >= TRF_COVER]
        bp = max(1, sum(x[4] for x in v))
        smp = rng.sample(v, min(N, len(v))); ok = 0; ch = 0; both = 0
        for c, s, e, lag, _ in smp:
            sq = CODE[np.frombuffer(fetch(c, s, min(e, s + 20000)), dtype=np.uint8)]
            o1 = (rec_at(sq, lag) - rec_at(dinuc_shuffle(sq, NR), lag)) > 0.02; o2 = (chain_at(sq, lag) - chain_at(dinuc_shuffle(sq, NR), lag)) > 0.01
            ok += o1; ch += o2; both += o1 and o2 and trf_support(trf, c, s, e, lag, MIN_PER)[1] >= TRF_COVER
        print(f'{TAG}\t{cls}\t{len(v)}\t{sum(x[4] for x in v) / 1e6:.2f}\t{len(t) / len(v):.2f} / {sum(x[4] for x in t) / bp:.2f}'
              f'\t{len(tc) / len(v):.2f} / {sum(x[4] for x in tc) / bp:.2f}\t{np.median([cv for _, cv in sup]):.2f}'
              f'\t{ok / len(smp):.2f}\t3-copy chain {ch / len(smp):.2f}\tchain+TRFcover {both / len(smp):.2f} (n={len(smp)})')
