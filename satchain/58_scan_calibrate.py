#!/usr/bin/env python3
"""
58_scan_calibrate.py  (v4.6; calibration and validation of the genome-wide tandem scan 57)

Null: ordinary genome windows (per scale) with no N, outside step-02 arrays and TRF calls with unit >= 7 bp; TEs and SSRs allowed
(a TE-free null keeps only the least repetitive windows and is biased low). Best lag < 7 bp counts as 0 (micro texture).
Tiers: p95, p99, p99.9 of the null max_excess (per scale). Tandem window: max_excess > tier.
Candidate regions: consecutive tandem windows (same scale, overlapping or adjacent) whose best lags agree within 10 % (or are 2x / 3x
multiples), merged; a region keeps its median lag and the max excess.
Validation:
  recall    bp of the validated HMM runs (56, p90) and of the step-36 runs covered by candidates, per tier
  new       candidate bp outside every step-02 array, by unit class (2-6, 7-39, 40-60, 61-1000, > 1000 bp) and TE overlap share
  scales    small (2 kb) vs large (10 kb) contributions
Outputs (<dest>/): 58_null_quantiles.tsv, 58_candidates_<tier>.bed, 58_validation.tsv
usage: 58_scan_calibrate.py <tag> <scan57 dir> <func genome folder> <pipeline out> <hmm56 dir> <dest>
"""
import sys, os, glob, gzip, collections, csv, bisect
import numpy as np
TAG, SC, FUNC, OUT, HMM, DEST = sys.argv[1:7]; os.makedirs(DEST, exist_ok=True)
opn = lambda p: gzip.open(p, 'rt') if p.endswith('.gz') else open(p)


def intervals(paths, kind):
    d = collections.defaultdict(list)
    for p in paths:
        if not os.path.exists(p): print('missing', p); continue
        for l in opn(p):
            if l.startswith(('#', 'track')) or not l.strip(): continue
            f = l.rstrip('\n').split('\t')
            if kind == 'gff':
                if len(f) < 9 or f[2] in ('repeat_region', 'target_site_duplication', 'long_terminal_repeat'): continue
                d[f[0]].append((int(f[3]) - 1, int(f[4])))
            elif kind == 'tsv02':
                if f[0] == 'array_id': continue
                d[f[1]].append((int(f[2]), int(f[3])))
            else: d[f[0]].append((int(f[1]), int(f[2])))
    out = {}
    for c, v in d.items():
        v.sort(); m = []
        for s, e in v:
            if m and s <= m[-1][1]: m[-1][1] = max(m[-1][1], e)
            else: m.append([s, e])
        a = np.array(m, dtype=np.int64); out[c] = (a[:, 0], a[:, 1])
    return out


def overlap_bp(iv, c, s, e):
    if c not in iv: return 0
    S, E = iv[c]; i = np.searchsorted(E, s, 'right'); t = 0
    while i < len(S) and S[i] < e: t += min(e, E[i]) - max(s, S[i]); i += 1
    return t


G = TAG; F = lambda *p: os.path.join(FUNC, *p)
te = intervals([F('_repeats_EDTA', f'fClaHyb_{G}.TEanno.sorted.gff3.gz')], 'gff')
trf = intervals([F('_TR', 'TRF', f'fClaHyb_{G}.TRF.mini_7-60bp.bed.gz'), F('_TR', 'TRF', f'fClaHyb_{G}.TRF.sat_gt60bp.bed.gz')], 'bed')   # known tandem, unit >= 7 bp
ssr = intervals(glob.glob(F('_TR', 'GMATA_SSR', f'fClaHyb_{G}.SSR.*.bed.gz')), 'bed')
a02 = intervals([os.path.join(OUT, f'fClaHyb_{G}_satellites', '02_satellite_arrays.tsv')], 'tsv02')
hmm = intervals([os.path.join(HMM, '56_runs_validated_q0.9.bed')], 'bed')
r36 = intervals([os.path.join(OUT, f'fClaHyb_{G}_satellites', '36_characterized_intervals.bed')], 'bed')
rows = collections.defaultdict(list)
for p in sorted(glob.glob(os.path.join(SC, f'scan57_{G}_job*.tsv'))):
    for l in open(p):
        if l.startswith('chrom'): continue
        c, sc, s, e, nf, ex, lag, e2, rf = l.rstrip('\n').split('\t')
        if ex == 'NA': continue
        rows[sc].append((c, int(s), int(e), float(nf), float(ex), int(lag)))
q = {}
with open(os.path.join(DEST, '58_null_quantiles.tsv'), 'w') as fh:
    fh.write('scale\tn_windows\tn_clean_null\tp95\tp99\tp99.9\n')
    for sc, R in rows.items():
        # null = ordinary genome: no N, outside step-02 arrays and TRF calls with unit >= 7 bp; TEs and SSRs allowed (excluding them
        # would keep only the least repetitive windows and bias the null low). Windows whose best lag is < 7 bp count as 0 (micro texture).
        null = [(x[4] if x[5] >= 7 else 0.0) for x in R if x[3] == 0 and not (overlap_bp(trf, x[0], x[1], x[2]) or overlap_bp(a02, x[0], x[1], x[2]))]
        v = np.array(null); q[sc] = {t: float(np.quantile(v, t)) for t in (0.95, 0.99, 0.999)}
        fh.write(f'{sc}\t{len(R)}\t{len(null)}\t{q[sc][0.95]:.4f}\t{q[sc][0.99]:.4f}\t{q[sc][0.999]:.4f}\n')
print(open(os.path.join(DEST, '58_null_quantiles.tsv')).read())


def compatible(a, b):
    if a == 0 or b == 0: return False
    r = max(a, b) / min(a, b); k = round(r)
    return abs(r - k) <= 0.1 * k and k <= 3


val = []
for tier in (0.95, 0.99, 0.999):
    regs = []
    for sc, R in rows.items():
        thr = q[sc][tier]; cur = None
        for c, s, e, nf, ex, lag in sorted(R):
            if ex <= thr or lag < 7: continue          # satellite / minisatellite candidates only; 2-6 bp left to TRF / GMATA
            if cur and cur[0] == c and s <= cur[2] and compatible(lag, cur[4][-1]):
                cur[2] = max(cur[2], e); cur[3] = max(cur[3], ex); cur[4].append(lag)
            else:
                if cur: regs.append((sc,) + tuple(cur))
                cur = [c, s, e, ex, [lag]]
        if cur: regs.append((sc,) + tuple(cur))
    # merge small + large scale regions into one candidate set
    cand = collections.defaultdict(list)
    for sc, c, s, e, ex, lags in regs: cand[c].append((s, e, sc, ex, int(np.median(lags))))
    with open(os.path.join(DEST, f'58_candidates_p{tier}.bed'), 'w') as fh:
        for c, v in cand.items():
            for s, e, sc, ex, lag in sorted(v): fh.write(f'{c}\t{s}\t{e}\t{sc}|lag{lag}|ex{ex:.3f}\t0\t.\n')
    civ = intervals([os.path.join(DEST, f'58_candidates_p{tier}.bed')], 'bed')
    tot = lambda iv: sum(int((E - S).sum()) for S, E in iv.values())
    def cover(a, b):
        t = 0
        for c, (S, E) in a.items():
            for s, e in zip(S, E): t += overlap_bp(b, c, int(s), int(e))
        return t
    cbp = tot(civ); new = collections.Counter(); newte = collections.Counter()
    for c, v in cand.items():
        for s, e, sc, ex, lag in v:
            out_bp = (e - s) - overlap_bp(a02, c, s, e)
            if out_bp <= 0: continue
            cls = '2-6' if lag <= 6 else '7-39' if lag <= 39 else '40-60' if lag <= 60 else '61-1000' if lag <= 1000 else '>1000'
            tb = overlap_bp(te, c, s, e)
            if tb > 0.5 * (e - s): cls = 'TE_derived'
            new[cls] += out_bp; newte[cls] += min(out_bp, tb)
    val.append([tier, round(cbp / 1e6, 2), round(cover(hmm, civ) / max(1, tot(hmm)), 3), round(cover(r36, civ) / max(1, tot(r36)), 3),
                round(cover(a02, civ) / max(1, tot(a02)), 3), round(sum(new.values()) / 1e6, 2)] +
               [f'{new[k] / 1e6:.2f}Mb' for k in ('7-39', '40-60', '61-1000', '>1000', 'TE_derived')])
with open(os.path.join(DEST, '58_validation.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['tier', 'candidate_Mb', 'recall_validated_hmm_bp', 'recall_step36_bp', 'share_of_02_interval_bp', 'new_outside_02_Mb',
                'new_7-39bp', 'new_40-60bp', 'new_61-1000bp', 'new_gt1000bp', 'new_TE_derived'])
    w.writerows(val)
for v in val: print(TAG, *v)
