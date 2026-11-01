#!/usr/bin/env python3
"""58c_new_arrays_02format.py  (2026-11-01, audit CB02 / CB12): the 02-format candidate table that 03_monomer_copies.py reads, built
from the 58b bed. The script that produced new_candidates_02format_{Gar,Mac}.tsv for the v4.6 catalogue was lost (scratchpad only);
its rule was inferred from the tables it left (chk58b.py [F]: 5,145 / 5,145 Gar and 4,214 / 4,214 Mac rows reproduced) and is coded
here explicitly:
  --route chain          keep new_array only (>= 3-copy chain excess; what docs/STEP_METHOD step 1 described)
  --route chain_or_trf   keep new_array OR (unit 7-39 bp AND TRF support)   <- what was run for v4.6 (default)
Tier column kept_by: new_array / new_pair+TRF / weak+TRF (route chain_or_trf; the two TRF tiers rest on the HMM + TRF alone).
TRF support (--trf-field): 'cover' = the 58b TRFcov field >= --trf-cover (CB03, compatible calls with period >= 7 bp covering the
candidate; needs a bed written by 58b >= 2026-11-01), 'asis' = the 58b TRF flag as written (pre-2026-11-01 beds: any compatible
call overlapping by >= 1 bp). The default is 'cover'; the script raises when the bed lacks TRFcov, so a reproduction of the v4.6
tables passes --trf-field asis.
Output columns = new_candidates_02format_<G>.tsv (array_id <G>_NEWnnnnn numbered along chrom / start, chrom, start_0based, end,
array_length_bp, period_mode_bp = 58 unit lag, candidate_tier = 58b tier, unit_class, TE_fraction, TRF_support, source) plus kept_by
and trf_cover. Overlapping kept intervals are counted and reported (03 decomposes each interval independently, CB06).
usage: 58c_new_arrays_02format.py <tag> <58b bed> <out tsv> [--route chain|chain_or_trf] [--trf-field cover|asis]
       [--trf-cover 0.5] [--tier 0.999]"""
import sys, argparse, collections, csv


def parse_bed(path):
    """58b bed -> rows (chrom, s, e, tier, lag, cls, te, trf_flag, trf_cover or None)"""
    rows = []
    for l in open(path):
        f = l.rstrip('\n').split('\t'); c, s, e = f[0], int(f[1]), int(f[2]); p = f[3].split('|')
        tier, lag, cls, te, trf = p[0], int(p[1][4:]), p[2], float(p[3][2:]), int(p[4][3:])
        cov = float(p[5][6:]) if len(p) > 5 and p[5].startswith('TRFcov') else None
        rows.append((c, s, e, tier, lag, cls, te, trf, cov))
    return rows


def keep_rule(tier, cls, trf, route):
    """-> kept_by tier or None"""
    if tier == 'new_array': return 'new_array'
    if route == 'chain_or_trf' and cls == '7-39' and trf: return f'{tier}+TRF'
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('tag'); ap.add_argument('bed'); ap.add_argument('out')
    ap.add_argument('--route', choices=('chain', 'chain_or_trf'), default='chain_or_trf')
    ap.add_argument('--trf-field', choices=('cover', 'asis'), default='cover')
    ap.add_argument('--trf-cover', type=float, default=0.5); ap.add_argument('--tier', default='0.999')
    a = ap.parse_args()
    rows = parse_bed(a.bed)
    if a.trf_field == 'cover' and any(r[8] is None for r in rows):
        raise SystemExit(f'58c: --trf-field cover needs the TRFcov field of a 58b bed written >= 2026-11-01 ({a.bed} lacks it); '
                         f'use --trf-field asis to apply the 58b flag as written')
    src = f'scan57_p{float(a.tier) * 100:g}'
    kept = []
    for c, s, e, tier, lag, cls, te, trf, cov in sorted(rows):
        sup = (cov >= a.trf_cover) if a.trf_field == 'cover' else bool(trf)
        kb = keep_rule(tier, cls, sup, a.route)
        if kb: kept.append((c, s, e, tier, lag, cls, te, int(sup), cov, kb))
    n_ovl = sum(1 for p, q in zip(kept, kept[1:]) if p[0] == q[0] and q[1] < p[2])
    stats = collections.Counter(); bp = collections.Counter()
    with open(a.out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['array_id', 'chrom', 'start_0based', 'end', 'array_length_bp', 'period_mode_bp', 'candidate_tier', 'unit_class',
                    'TE_fraction', 'TRF_support', 'source', 'kept_by', 'trf_cover'])
        for i, (c, s, e, tier, lag, cls, te, sup, cov, kb) in enumerate(kept, 1):
            w.writerow([f'{a.tag}_NEW{i:05d}', c, s, e, e - s, lag, tier, cls, f'{te:.2f}', sup, src, kb, 'NA' if cov is None else f'{cov:.2f}'])
            stats[kb] += 1; bp[kb] += e - s
    print(f'{a.tag}: {len(rows)} 58b candidates -> {len(kept)} kept (route {a.route}, TRF field {a.trf_field}'
          f'{f" >= {a.trf_cover}" if a.trf_field == "cover" else ""}); ' + ', '.join(f'{k} {stats[k]} ({bp[k] / 1e6:.2f} Mb)' for k in sorted(stats))
          + f'; {n_ovl} kept intervals overlap the previous kept interval (CB06); wrote {a.out}')


if __name__ == '__main__':
    main()
