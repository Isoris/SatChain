#!/usr/bin/env python3
"""
58_scan_calibrate.py  (v4.6; calibration and validation of the genome-wide tandem scan 57)

Null: ordinary genome windows (per scale) with no N, outside step-02 arrays and TRF calls with unit >= 7 bp; TEs and SSRs allowed
(a TE-free null keeps only the least repetitive windows and is biased low).
  2026-11-01 (CB04): the tiers are quantiles of the null windows whose best lag is >= 7 bp (the only windows that can become
  candidates). The pre-2026-11-01 coding (best lag < 7 counted as 0, which made p99.9 ~p99.7 of the lag >= 7 windows) is still
  written as the *_coded0 columns of 58_null_quantiles.tsv; --null coded0 restores it for the tiers.
Tiers: p95, p99, p99.9 of the null max_excess (per scale). Tandem window: max_excess > tier and best lag >= 7 (57 scans lags 2-1,000;
2-6 bp are left to TRF / GMATA).
Candidate regions: consecutive tandem windows (same scale, overlapping or adjacent) whose best lags agree within 10 % (or are 2x / 3x
multiples), merged; 2026-11-01 (CB05) the region's unit is the lag of the window with the maximum excess (the median lag of the
members is written beside it: a median over mixed multiples was a unit present in no window). 2026-11-01 (CB06): small- and
large-scale regions that overlap and have compatible lags are merged into one candidate (scale small+large) before the beds are
written; the per-region sum and the union Mb are both reported. Overlapping regions with incompatible lags stay separate (counted).
Validation:
  recall    bp of the validated HMM runs (56, p90) and of the step-36 runs covered by candidates, per tier (NA when the file is absent)
  new       candidate bp outside every step-02 array, by unit class (7-39, 40-60, 61-1000, > 1000 bp); 2026-11-01 (CB20) TE_derived
            when EDTA covers > 50 % of the part OUTSIDE the step-02 arrays (the part that is counted)
  scales    small (2 kb) vs large (10 kb) contributions
Annotation files are resolved with satpaths.func_path (2026-11-01, CB16): a missing mask raises.
Outputs (<dest>/): 58_null_quantiles.tsv, 58_candidates_<tier>.bed (name: scale|lag<unit>|ex<max excess>|medlag<median>|n<windows>),
58_validation.tsv
usage: 58_scan_calibrate.py <tag> <scan57 dir> <func genome folder> <pipeline out> <hmm56 dir> <dest> [--null lag7|coded0]
"""
import sys, os, glob, gzip, collections, csv, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from satpaths import func_path
opn = lambda p: gzip.open(p, 'rt') if p.endswith('.gz') else open(p)
MIN_LAG = 7


def merge_sorted(v):
    """[(s, e), ...] -> merged, sorted list of [s, e]"""
    v.sort(); m = []
    for s, e in v:
        if m and s <= m[-1][1]: m[-1][1] = max(m[-1][1], e)
        else: m.append([s, e])
    return m


def intervals(paths, kind, required=True):
    """merged intervals per chromosome from bed / gff / 02 tables -> {chrom: (starts, ends)}; a missing file raises unless
    required=False (then it is skipped with a warning; used only for the recall files, which are validation not masks)"""
    d = collections.defaultdict(list)
    for p in paths:
        if not os.path.exists(p):
            if required: raise FileNotFoundError(f'58: required interval file missing: {p}')
            print(f'    WARNING: optional file missing, skipped: {p}'); continue
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
        a = np.array(merge_sorted(v), dtype=np.int64); out[c] = (a[:, 0], a[:, 1])
    return out


def overlap_bp(iv, c, s, e):
    if c not in iv: return 0
    S, E = iv[c]; i = np.searchsorted(E, s, 'right'); t = 0
    while i < len(S) and S[i] < e: t += min(e, E[i]) - max(s, S[i]); i += 1
    return t


def subtract(iv, c, s, e):
    """pieces of [s, e) not covered by iv[c]"""
    if c not in iv: return [(s, e)]
    S, E = iv[c]; i = np.searchsorted(E, s, 'right'); out = []; cur = s
    while i < len(S) and S[i] < e:
        if S[i] > cur: out.append((cur, int(S[i])))
        cur = max(cur, int(E[i])); i += 1
    if cur < e: out.append((cur, e))
    return out


def compatible(a, b):
    if a == 0 or b == 0: return False
    r = max(a, b) / min(a, b); k = round(r)
    return abs(r - k) <= 0.1 * k and k <= 3


def regions_of_scale(R, thr, sc):
    """consecutive tandem windows of one scale with compatible lags -> [chrom, s, e, max_ex, unit(lag of max-excess window),
    median lag, n windows, scale]"""
    regs = []; cur = None
    for c, s, e, nf, ex, lag in sorted(R):
        if ex <= thr or lag < MIN_LAG: continue          # satellite / minisatellite candidates only; 2-6 bp left to TRF / GMATA
        if cur and cur[0] == c and s <= cur[2] and compatible(lag, cur[4][-1]):
            cur[2] = max(cur[2], e); cur[3].append(ex); cur[4].append(lag)
        else:
            if cur: regs.append(_close(cur, sc))
            cur = [c, s, e, [ex], [lag]]
    if cur: regs.append(_close(cur, sc))
    return regs


def _close(cur, sc):
    c, s, e, exs, lags = cur; j = int(np.argmax(exs))
    return [c, s, e, float(exs[j]), int(lags[j]), int(np.median(lags)), len(lags), sc]


def union_regions(regs):
    """CB06: merge overlapping regions (any scale) with compatible units; returns (merged list, n overlapping pairs left)"""
    by = collections.defaultdict(list)
    for r in regs: by[r[0]].append(r)
    out = []; left = 0
    for c, v in by.items():
        v.sort(key=lambda r: (r[1], r[2])); cur = None
        for r in v:
            if cur and r[1] < cur[2]:
                if compatible(r[4], cur[4]):
                    if r[3] > cur[3]: cur[3], cur[4], cur[5] = r[3], r[4], r[5]
                    cur[2] = max(cur[2], r[2]); cur[6] += r[6]
                    cur[7] = cur[7] if cur[7] == r[7] else 'small+large'
                    continue
                left += 1
            if cur: out.append(cur)
            cur = list(r)
        if cur: out.append(cur)
    return out, left


def main():
    ap = argparse.ArgumentParser()
    for a in ('tag', 'scan57_dir', 'func', 'out', 'hmm56_dir', 'dest'): ap.add_argument(a)
    ap.add_argument('--null', choices=('lag7', 'coded0'), default='lag7', help='tiers from lag >= 7 null windows (default) or the pre-2026-11-01 zero-coded null')
    a = ap.parse_args(); TAG, SC, FUNC, OUT, HMM, DEST = a.tag, a.scan57_dir, a.func, a.out, a.hmm56_dir, a.dest
    os.makedirs(DEST, exist_ok=True); G = TAG
    te = intervals([func_path(FUNC, G, '_repeats_EDTA', f'fClaHyb_{G}.TEanno.sorted.gff3.gz')], 'gff')
    trf = intervals([func_path(FUNC, G, '_TR', f'TRF/fClaHyb_{G}.TRF.mini_7-60bp.bed.gz'),
                     func_path(FUNC, G, '_TR', f'TRF/fClaHyb_{G}.TRF.sat_gt60bp.bed.gz')], 'bed')   # known tandem, unit >= 7 bp
    # GMATA SSR beds are not a mask (SSRs are allowed in the null); the former unused load was removed 2026-11-01
    a02 = intervals([os.path.join(OUT, f'fClaHyb_{G}_satellites', '02_satellite_arrays.tsv')], 'tsv02')
    hmm = intervals([os.path.join(HMM, '56_runs_validated_q0.9.bed')], 'bed', required=False)      # CB17: 56 may run after 58
    r36 = intervals([os.path.join(OUT, f'fClaHyb_{G}_satellites', '36_characterized_intervals.bed')], 'bed', required=False)
    rows = collections.defaultdict(list)
    files = sorted(glob.glob(os.path.join(SC, f'scan57_{G}_job*.tsv')))
    if not files: raise FileNotFoundError(f'no scan57_{G}_job*.tsv in {SC}')
    for p in files:
        for l in open(p):
            if l.startswith('chrom'): continue
            c, sc, s, e, nf, ex, lag, e2, rf = l.rstrip('\n').split('\t')
            if ex == 'NA': continue
            rows[sc].append((c, int(s), int(e), float(nf), float(ex), int(lag)))
    q = {}; TIERS = (0.95, 0.99, 0.999)
    with open(os.path.join(DEST, '58_null_quantiles.tsv'), 'w') as fh:
        fh.write('scale\tn_windows\tn_clean_null\tn_clean_null_lag7\tp95_coded0\tp99_coded0\tp99.9_coded0\tp95\tp99\tp99.9\tnull_used\n')
        for sc, R in rows.items():
            # null = ordinary genome: no N, outside step-02 arrays and TRF calls with unit >= 7 bp; TEs and SSRs allowed (excluding them
            # would keep only the least repetitive windows and bias the null low).
            clean = [x for x in R if x[3] == 0 and not (overlap_bp(trf, x[0], x[1], x[2]) or overlap_bp(a02, x[0], x[1], x[2]))]
            coded0 = np.array([(x[4] if x[5] >= MIN_LAG else 0.0) for x in clean])      # pre-2026-11-01: best lag < 7 -> 0
            lag7 = np.array([x[4] for x in clean if x[5] >= MIN_LAG])                   # CB04: lag >= 7 windows only
            q0 = {t: float(np.quantile(coded0, t)) for t in TIERS}; q7 = {t: float(np.quantile(lag7, t)) for t in TIERS}
            q[sc] = q7 if a.null == 'lag7' else q0
            fh.write(f'{sc}\t{len(R)}\t{len(clean)}\t{len(lag7)}\t' + '\t'.join(f'{q0[t]:.4f}' for t in TIERS) + '\t'
                     + '\t'.join(f'{q7[t]:.4f}' for t in TIERS) + f'\t{a.null}\n')
    print(open(os.path.join(DEST, '58_null_quantiles.tsv')).read())

    val = []
    for tier in TIERS:
        regs = []
        for sc, R in rows.items(): regs += regions_of_scale(R, q[sc][tier], sc)
        sum_bp = sum(r[2] - r[1] for r in regs)
        cand, n_left = union_regions(regs)                     # CB06: small + large regions of compatible lag -> one candidate
        union_bp = sum(r[2] - r[1] for r in cand)
        with open(os.path.join(DEST, f'58_candidates_p{tier}.bed'), 'w') as fh:
            for c, s, e, ex, lag, med, nw, sc in sorted(cand): fh.write(f'{c}\t{s}\t{e}\t{sc}|lag{lag}|ex{ex:.3f}|medlag{med}|n{nw}\t0\t.\n')
        civ = intervals([os.path.join(DEST, f'58_candidates_p{tier}.bed')], 'bed')
        tot = lambda iv: sum(int((E - S).sum()) for S, E in iv.values())
        def cover(x, b):
            t = 0
            for c, (S, E) in x.items():
                for s, e in zip(S, E): t += overlap_bp(b, c, int(s), int(e))
            return t
        cbp = tot(civ); new = collections.Counter(); newte = collections.Counter(); new_sum = 0
        for c, s, e, ex, lag, med, nw, sc in cand:
            pieces = subtract(a02, c, s, e); out_bp = sum(pe - ps for ps, pe in pieces)
            if out_bp <= 0: continue
            new_sum += out_bp
            cls = '2-6' if lag <= 6 else '7-39' if lag <= 39 else '40-60' if lag <= 60 else '61-1000' if lag <= 1000 else '>1000'
            tb = sum(overlap_bp(te, c, ps, pe) for ps, pe in pieces)       # CB20: TE share of the part outside step-02
            if tb > 0.5 * out_bp: cls = 'TE_derived'
            new[cls] += out_bp; newte[cls] += tb
        new_union = sum((e - s) - overlap_bp(a02, c, int(s), int(e)) for c, (S, E) in civ.items() for s, e in zip(S, E))
        rec_h = round(cover(hmm, civ) / tot(hmm), 3) if hmm else 'NA'; rec_36 = round(cover(r36, civ) / tot(r36), 3) if r36 else 'NA'
        val.append([tier, round(sum_bp / 1e6, 2), round(union_bp / 1e6, 2), round(cbp / 1e6, 2), len(regs), len(cand), n_left, rec_h, rec_36,
                    round(cover(a02, civ) / max(1, tot(a02)), 3), round(new_sum / 1e6, 2), round(new_union / 1e6, 2)] +
                   [f'{new[k] / 1e6:.2f}Mb' for k in ('7-39', '40-60', '61-1000', '>1000', 'TE_derived')])
        print(f'{TAG} p{tier}: {len(regs)} per-scale regions ({sum_bp / 1e6:.2f} Mb sum) -> {len(cand)} candidates ({union_bp / 1e6:.2f} Mb union; '
              f'{n_left} overlapping pairs with incompatible lags left); new outside 02: sum {new_sum / 1e6:.2f} / union {new_union / 1e6:.2f} Mb')
    with open(os.path.join(DEST, '58_validation.tsv'), 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['tier', 'region_Mb_sum', 'candidate_Mb_union', 'candidate_Mb', 'n_regions', 'n_candidates', 'n_overlapping_incompatible',
                    'recall_validated_hmm_bp', 'recall_step36_bp', 'share_of_02_interval_bp', 'new_outside_02_Mb_sum', 'new_outside_02_Mb',
                    'new_7-39bp', 'new_40-60bp', 'new_61-1000bp', 'new_gt1000bp', 'new_TE_derived'])
        w.writerows(val)
    for v in val: print(TAG, *v)


if __name__ == '__main__':
    main()
