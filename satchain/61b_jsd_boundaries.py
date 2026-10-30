#!/usr/bin/env python3
"""
61b_jsd_boundaries.py  (SatChain v4.7 draft; composition-change boundaries with a fixed-midpoint Jensen-Shannon divergence)

For every run edge that faces an interruption (OUT = end of a run, IN = start of the next run of the same array, read in genomic order; arrays
never chained across an assembly gap), candidate boundary positions x are scanned within +-S bp of the catalogue edge. At each x the midpoint is
FIXED: a left window [x - w, x) and a right window [x, x + w) of equal width w. P, Q = frequencies of overlapping k-mers in the two windows (same
strand). JSD(P, Q) = H(M) - (H(P) + H(Q)) / 2 with M = (P + Q) / 2, log2, so 0 <= JSD <= 1; D = sqrt(JSD) is a metric (Jensen-Shannon
distance) and is what is reported.
Reading D: inside a satellite the two windows are copies -> D low; inside unique sequence they share few k-mers -> D high at any k. A boundary
is therefore a STEP, not a peak: at a run end b, D ramps from the satellite level to the spacer level over [b - w, b] (the right window enters
the spacer) and stays high; at a run start b it stays high until b and ramps down over [b, b + w]. The edge = the end of the ramp on the
plateau side, fitted by least squares (low level, ramp of width w with shape f^gamma, gamma 0.25 / 0.5 / 1, high level; both directions tried, best SSE). argmax(D) is kept as
a column for comparison only (in unique sequence it lands anywhere on the plateau).
Stop at gaps: a position x is not scored when either window contains a non-ACGT base or touches an assembly gap / contig join (gaps.py).
Null: D at random x deep inside the same array's runs (both windows satellite, >= w + S from any run edge), per (k, w). A junction is a
composition boundary when the fitted high level > null p99 AND the step (high - low) > the null spread (p99 - p50); else 'no_change'.
Window width w = m x the array's unit length (catalogue unit_bp; --w-units m, default 1,2,3): a window shorter than one unit holds no full
copy, so even inside the satellite the two windows differ and D is high (tested: w = 100 bp on a 220-bp unit fails). --w (bp) is used only
when unit_bp is missing. Scan range S = max(--search, max w + 200).
Swept: k (--k), w (--w-units). Per junction: edge per setting, spread of the called edges across settings (stability), offset from the catalogue edge.
Outputs (<dest>/): 61b_jsd_junctions_<tag>.tsv (junction x k x w), 61b_jsd_stability_<tag>.tsv (per junction), 61b_jsd_null_<tag>.tsv,
61b_jsd_summary_<tag>.tsv, optional 61b_jsd_profiles_<tag>.tsv (D along x for --profile-arrays).
usage: 61b_jsd_boundaries.py --catalogue 59_catalogue_<G>.tsv --genome <G>.fa --gaps <G>.gaps.bed --dest <dir> --tag <G>
                             [--arrays id1,id2] [--min-runs 2] [--max-gap 50000] [--k 3,4,5,6] [--w 100,200,400] [--search 600] [--null 200]
"""
import argparse, collections, csv, os, random, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gaps import Gaps
csv.field_size_limit(1 << 30)
CODE = np.full(256, -1, dtype=np.int64)
for i, ch in enumerate(b'ACGT'): CODE[ch] = i; CODE[ch + 32] = i


class Fasta:
    def __init__(self, path):
        self.fh = open(path, 'rb'); self.idx = {}
        for l in open(path + '.fai'):
            n, ln, off, lb, lby = l.split('\t')[:5]; self.idx[n] = (int(ln), int(off), int(lb), int(lby))
    def get(self, c, a, b):
        L, off, lb, lby = self.idx[c]; a = max(0, a); b = min(L, b)
        if b <= a: return b''
        self.fh.seek(off + (a // lb) * lby + a % lb)
        return self.fh.read((b - a) + ((b - a) // lb + 2) * (lby - lb)).replace(b'\n', b'').replace(b'\r', b'')[:b - a]


def kmer_codes(seq, k):
    """code of the k-mer starting at each position; -1 when it contains a non-ACGT base"""
    s = CODE[np.frombuffer(seq, dtype=np.uint8)]; n = len(s) - k + 1
    if n <= 0: return np.zeros(0, dtype=np.int64)
    code = np.zeros(n, dtype=np.int64); bad = np.zeros(n, dtype=bool)
    for j in range(k):
        x = s[j:j + n]; bad |= x < 0; code = code * 4 + np.maximum(x, 0)
    code[bad] = -1
    return code


def _flog(c): return c * np.log2(c) if c > 0 else 0.0


def jsd_scan(seq, k, w, xs):
    """sqrt(JSD) between k-mer spectra of [x - w, x) and [x, x + w) for each x (positions in seq); nan when a window has a bad k-mer.
    Incremental: counts and the sums S = sum c log2 c of the left, right and pooled (M) spectra are updated as x moves by one base, so each
    step costs O(1) instead of O(4^k). H(P) = log2 N - S_P / N; H(M) = log2 2N - S_M / 2N (N = k-mers per window)."""
    codes = kmer_codes(seq, k); n = len(codes); xs = np.asarray(xs, dtype=np.int64); out = np.full(len(xs), np.nan)
    if n == 0 or len(xs) == 0: return out
    badcum = np.concatenate([[0], np.cumsum(codes < 0)]); N = w - k + 1
    cl = collections.Counter(); cr = collections.Counter(); sl = sr = sm = 0.0; prev = None
    def add(cnt, other, c, d):                      # change count of code c in one window by d; return deltas of (own sum, pooled sum)
        o, m = cnt[c], cnt[c] + other[c]; cnt[c] = o + d
        return _flog(o + d) - _flog(o), _flog(m + d) - _flog(m)
    for i, x in enumerate(xs):
        l0, l1, r0, r1 = x - w, x - k + 1, x, x + w - k + 1
        if l0 < 0 or r1 > n or N <= 0: prev = None; continue
        if prev is None or x != prev + 1:           # (re)initialise
            cl = collections.Counter(codes[l0:l1].tolist()); cr = collections.Counter(codes[r0:r1].tolist())
            sl = sum(_flog(v) for c, v in cl.items() if c >= 0); sr = sum(_flog(v) for c, v in cr.items() if c >= 0)
            sm = sum(_flog(cl[c] + cr[c]) for c in set(cl) | set(cr) if c >= 0)
        else:                                       # slide by one: left loses codes[l0 - 1], gains codes[l1 - 1]; right loses codes[r0 - 1], gains codes[r1 - 1]
            for cnt, other, c, d, side in ((cl, cr, codes[l0 - 1], -1, 'l'), (cl, cr, codes[l1 - 1], +1, 'l'), (cr, cl, codes[r0 - 1], -1, 'r'), (cr, cl, codes[r1 - 1], +1, 'r')):
                if c < 0: cnt[c] += d; continue
                a_, m_ = add(cnt, other, c, d)
                if side == 'l': sl += a_
                else: sr += a_
                sm += m_
        prev = x
        if badcum[l1] - badcum[l0] or badcum[r1] - badcum[r0]: continue
        hp = np.log2(N) - sl / N; hq = np.log2(N) - sr / N; hm = np.log2(2 * N) - sm / (2 * N)
        out[i] = np.sqrt(max(0.0, hm - (hp + hq) / 2))
    return out


def ramp_fit(x, d, w):
    """least-squares fit of D(x) = lo + (hi - lo) * f(x)^gamma, f = ramp of width w ending (up) or starting (down) at b; gamma 0.25 / 0.5 / 1.
    Vectorised over candidate b (every 4 bp, then every bp around the best). Returns b, lo, hi, sse, shape."""
    ok = ~np.isnan(d); x, d = x[ok], d[ok]
    if len(x) < 10: return None
    def fit(bs):
        best = None
        for direc in ('up', 'down'):
            F0 = np.clip((x[None, :] - (bs[:, None] - w)) / w, 0, 1) if direc == 'up' else np.clip((bs[:, None] + w - x[None, :]) / w, 0, 1)
            for gam in (0.25, 0.5, 1.0):
                F = F0 ** gam; Gm = 1 - F
                a11, a12, a22 = (Gm * Gm).sum(1), (Gm * F).sum(1), (F * F).sum(1); y1, y2 = Gm @ d, F @ d
                det = a11 * a22 - a12 ** 2; good = np.abs(det) > 1e-9
                lo = np.where(good, (a22 * y1 - a12 * y2) / np.where(good, det, 1), np.nan); hi = np.where(good, (a11 * y2 - a12 * y1) / np.where(good, det, 1), np.nan)
                sse = ((d[None, :] - lo[:, None] * Gm - hi[:, None] * F) ** 2).sum(1)
                sse[~good | ~(hi > lo)] = np.inf
                j = int(np.argmin(sse))
                if np.isfinite(sse[j]) and (best is None or sse[j] < best[3]): best = (int(bs[j]), float(lo[j]), float(hi[j]), float(sse[j]), f'{direc}^{gam}')
        return best
    best = fit(x[::4].astype(float))
    if best: best = fit(np.arange(best[0] - 4, best[0] + 5, dtype=float)) or best
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--catalogue', required=True); ap.add_argument('--genome', required=True); ap.add_argument('--gaps', default='')
    ap.add_argument('--scan-n', action='store_true', help='also take N runs from the FASTA as gaps')
    ap.add_argument('--dest', required=True); ap.add_argument('--tag', required=True); ap.add_argument('--arrays', default='')
    ap.add_argument('--min-runs', type=int, default=2); ap.add_argument('--max-gap', type=int, default=50000)
    ap.add_argument('--k', default='3,4,5,6'); ap.add_argument('--w-units', default='1,2,3'); ap.add_argument('--w', default='200,400,600')
    ap.add_argument('--search', type=int, default=600)
    ap.add_argument('--null', type=int, default=200); ap.add_argument('--profile-arrays', default=''); ap.add_argument('--seed', type=int, default=61)
    a = ap.parse_args(); os.makedirs(a.dest, exist_ok=True); rng = random.Random(a.seed)
    KS = [int(x) for x in a.k.split(',')]; WU = [float(x) for x in a.w_units.split(',')]; WBP = [int(x) for x in a.w.split(',')]
    gaps = Gaps(a.gaps or None, a.genome if a.scan_n else None)
    fa = Fasta(a.genome); keep = set(a.arrays.split(',')) if a.arrays else None; prof = set(a.profile_arrays.split(',')) if a.profile_arrays else set()
    runs = collections.defaultdict(list)
    for r in csv.DictReader(open(a.catalogue), delimiter='\t'):
        if r.get('element_class', 'satellite') != 'satellite': continue
        if keep and r['array_id'] not in keep: continue
        runs[r['array_id']].append(r)
    chains, n_gap_breaks = [], 0
    for aid, rs in runs.items():
        rs.sort(key=lambda r: int(r['start0'])); cur = [rs[0]]
        for r in rs[1:]:
            d0, d1 = int(cur[-1]['end']), int(r['start0'])
            if gaps.between(r['chrom'], d0, d1): n_gap_breaks += 1
            if 0 < d1 - d0 <= a.max_gap and not gaps.between(r['chrom'], d0, d1): cur.append(r)
            else:
                if len(cur) >= a.min_runs: chains.append((aid, cur))
                cur = [r]
        if len(cur) >= a.min_runs: chains.append((aid, cur))
    J, NUL, PRO = [], [], []
    for ci, (aid, rs) in enumerate(chains):
        chrom = rs[0]['chrom']
        u = int(float(rs[0].get('unit_bp') or 0))
        WS = [int(round(m * u)) for m in WU] if u >= 20 else WBP; W = max(WS); S = max(a.search, W + 200)
        # null positions: deep inside runs
        inner = [(int(r['start0']) + W + S, int(r['end']) - W - S) for r in rs if int(r['end']) - int(r['start0']) > 2 * (W + S) + 50]
        null = {}
        if inner:
            pts = [rng.randint(lo, hi) for _ in range(a.null) for lo, hi in [rng.choice(inner)]]
            for k in KS:
                for w in WS:
                    vals = []
                    for p in pts:
                        if gaps.between(chrom, p - w, p + w): continue
                        seq = fa.get(chrom, p - w, p + w); v = jsd_scan(seq, k, w, [w])[0]
                        if not np.isnan(v): vals.append(v)
                    null[(k, w)] = (np.percentile(vals, 99), np.percentile(vals, 50)) if len(vals) >= 20 else (np.nan, np.nan)
                    NUL.append(dict(chain=f'{aid}#{ci}', array_id=aid, unit_bp=u, k=k, w=w, n=len(vals), p50=round(np.median(vals), 4) if vals else '', p99=round(null[(k, w)][0], 4) if vals else ''))
        for i, (r1, r2) in enumerate(zip(rs, rs[1:])):
            for jt, edge, run in (('OUT', int(r1['end']), r1['run_id']), ('IN', int(r2['start0']), r2['run_id'])):
                lo, hi = edge - S - W, edge + S + W
                seq = fa.get(chrom, lo, hi)
                for k in KS:
                    for w in WS:
                        xs = np.arange(edge - S, edge + S + 1) - lo
                        # stop at gaps: drop x whose windows touch a gap / join
                        ok = np.array([not gaps.between(chrom, lo + x - w, lo + x + w) for x in xs[::10]])
                        ok = np.repeat(ok, 10)[:len(xs)]
                        d = jsd_scan(seq, k, w, xs); d[~ok] = np.nan
                        if np.all(np.isnan(d)):
                            J.append(dict(chain=f'{aid}#{ci}', array_id=aid, family=r1.get('family', ''), chrom=chrom, run_id=run, gap_index=i, type=jt, k=k, w=w,
                                          catalogue_edge=edge, jsd_edge='', offset='', peak_D='', null_p99='', call='gap_or_N')); continue
                        am = int(np.nanargmax(d)); thr, med = null.get((k, w), (np.nan, np.nan)); fit = ramp_fit(xs.astype(float), d, w)
                        if fit is None: call, bb, lv, hv, dr = 'too_few_positions', None, np.nan, np.nan, ''
                        else:
                            bb, lv, hv, _, dr = fit
                            call = 'no_null' if np.isnan(thr) else ('boundary' if hv > thr and hv - lv > thr - med else 'no_change')
                        J.append(dict(chain=f'{aid}#{ci}', array_id=aid, family=r1.get('family', ''), chrom=chrom, run_id=run, gap_index=i, type=jt, k=k, w=w,
                                      w_units=round(w / u, 2) if u else '', catalogue_edge=edge, jsd_edge=lo + int(bb) if bb is not None else '', offset=lo + int(bb) - edge if bb is not None else '',
                                      ramp=dr, D_low=round(lv, 4), D_high=round(hv, 4), argmax_edge=lo + int(xs[am]),
                                      null_p99=round(thr, 4) if not np.isnan(thr) else '', scored_frac=round(float(np.mean(~np.isnan(d))), 3), call=call))
                        if aid in prof:
                            for x, v in zip(xs[::5], d[::5]): PRO.append(dict(array_id=aid, run_id=run, type=jt, k=k, w=w, pos=lo + int(x), D=round(float(v), 4) if not np.isnan(v) else ''))
        if ci % 50 == 0: print(f'{ci}/{len(chains)} chains', flush=True)
    def write(name, rows):
        if not rows: return
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with open(os.path.join(a.dest, f'61b_jsd_{name}_{a.tag}.tsv'), 'w', newline='') as fh:
            wr = csv.DictWriter(fh, fieldnames=keys, delimiter='\t', lineterminator='\n'); wr.writeheader(); wr.writerows(rows)
    write('junctions', J); write('null', NUL); write('profiles', PRO)
    # stability per junction: called edges across settings
    st = collections.defaultdict(list)
    for r in J:
        if r['call'] == 'boundary': st[(r['chain'], r['run_id'], r['type'])].append(r['jsd_edge'])
    allj = {(r['chain'], r['run_id'], r['type']): r for r in J}
    nset = len(KS) * len(WS); ST = []
    for key, r in allj.items():
        e = st.get(key, [])
        ST.append(dict(chain=key[0], array_id=r['array_id'], run_id=key[1], type=key[2], catalogue_edge=r['catalogue_edge'], settings=nset, settings_boundary=len(e),
                       median_jsd_edge=int(np.median(e)) if e else '', spread_bp=int(max(e) - min(e)) if e else '', median_offset=int(np.median(e) - r['catalogue_edge']) if e else ''))
    write('stability', ST)
    calls = collections.Counter(r['call'] for r in J)
    SM = [dict(tag=a.tag, chains=len(chains), chain_links_refused_at_gaps=n_gap_breaks, junctions=len(ST), junction_settings=len(J), **{f'calls_{k}': v for k, v in calls.items()},
               junctions_boundary_all_settings=sum(1 for x in ST if x['settings_boundary'] == nset),
               median_spread_bp=float(np.median([x['spread_bp'] for x in ST if x['spread_bp'] != ''])) if any(x['spread_bp'] != '' for x in ST) else '')]
    write('summary', SM); print(SM[0])


if __name__ == '__main__':
    main()