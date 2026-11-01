#!/usr/bin/env python3
"""58b: final list of NEW tandem candidates (scan 57, tier p99.9, outside every step-02 array, lag >= 7 bp).
Each candidate (first --test-kb kb tested, default 20): recurrence excess at its lag and 3-copy chain fraction, both against a
dinucleotide-block shuffle;
  new_array   chain excess > --chain-cut (0.01; >= 3 tandem copies)     new_pair  recurrence excess > --rec-cut (0.02) but no chain
  weak        neither
plus TRF support and TE overlap share (EDTA). 2026-11-01:
  CB03  TRF support is written twice: trf_asis (pre-2026-11-01 flag: any compatible call, any period, >= 1 bp overlap) and
        trf_cover (fraction of the candidate covered by compatible calls with period >= 7 bp); the bed TRF field is
        trf_cover >= --trf-cover (0.5), the TRFcov field the fraction.
  CB09  the EDTA share is computed on MERGED EDTA intervals (58.intervals); nested INT / LTR features no longer sum past 1.
  CB11  58b_tier_sweep.tsv: tiers for every combination of chain cut 0.005 / 0.01 / 0.02, recurrence cut 0.005 / 0.01 / 0.02 and
        test length 10 / 20 / 40 kb (n, Mb, TRF-covered share per class and tier); the chained-position count is written too.
  CA10  helpers are imported from val58_new (no exec of truncated source); CB16 annotation paths via satpaths (raise when missing).
Outputs: 58b_new_candidates.bed (BED9, coloured by tier; name tier|unit|class|TE|TRF|TRFcov), 58b_new_candidates.tsv (all
statistics per candidate: the input of 58c), 58b_tier_sweep.tsv, and a summary on stdout.
usage: 58b_new_candidates.py <tag> <cal58 dir> <func genome folder> <pipeline out> <genome fasta> [--tier 0.999] [--test-kb 20]
       [--chain-cut 0.01] [--rec-cut 0.02] [--trf-cover 0.5] [--no-sweep]"""
import sys, os, collections, csv, argparse, importlib.util
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import val58_new as V
from satpaths import func_path
_s = importlib.util.spec_from_file_location('scan_calibrate_58', os.path.join(HERE, '58_scan_calibrate.py'))
M58 = importlib.util.module_from_spec(_s); _s.loader.exec_module(M58)        # intervals() = merged intervals, overlap_bp()
COL = {'new_array': '31,119,180', 'new_pair': '235,104,52', 'weak': '180,178,169'}
SWEEP_CUTS = (0.005, 0.01, 0.02); SWEEP_KB = (10, 20, 40)


def tier_of(ch, r, chain_cut, rec_cut): return 'new_array' if ch > chain_cut else 'new_pair' if r > rec_cut else 'weak'


def main():
    ap = argparse.ArgumentParser()
    for x in ('tag', 'cal', 'func', 'out', 'fa'): ap.add_argument(x)
    ap.add_argument('--tier', default='0.999'); ap.add_argument('--test-kb', type=int, default=20)
    ap.add_argument('--chain-cut', type=float, default=0.01); ap.add_argument('--rec-cut', type=float, default=0.02)
    ap.add_argument('--trf-cover', type=float, default=0.5); ap.add_argument('--no-sweep', action='store_true')
    a = ap.parse_args(); TAG, CAL, FUNC, OUT, FA = a.tag, a.cal, a.func, a.out, a.fa
    NR = np.random.default_rng(58)
    te = M58.intervals([func_path(FUNC, TAG, '_repeats_EDTA', f'fClaHyb_{TAG}.TEanno.sorted.gff3.gz')], 'gff')   # CB09: merged
    a02 = V.load_a02(OUT, TAG); trf = V.load_trf(FUNC, TAG); fetch = V.make_fetch(FA); classes = V.load_classes(CAL, a.tier, a02)
    kbs = sorted(set(SWEEP_KB) | {a.test_kb}) if not a.no_sweep else [a.test_kb]
    tot = collections.Counter(); n = collections.Counter(); trf_ok = collections.Counter(); trf_as = collections.Counter()
    sweep = collections.defaultdict(lambda: [0, 0, 0])        # (kb, chain_cut, rec_cut, cls, tier) -> [n, bp, trf covered]
    with open(f'{CAL}/58b_new_candidates.bed', 'w') as bed, open(f'{CAL}/58b_new_candidates.tsv', 'w', newline='') as tsv:
        w = csv.writer(tsv, delimiter='\t', lineterminator='\n')
        w.writerow(['chrom', 'start_0based', 'end', 'tier', 'lag', 'unit_class', 'out_of_02_bp', 'te_share_merged', 'trf_asis', 'trf_cover',
                    'trf_support', 'rec_excess', 'chain_excess', 'n_chained_positions', 'n_kmer_positions', 'test_bp'])
        for cls, v in classes.items():
            for c, s, e, lag, out in v:
                stat = {}
                for kb in kbs:
                    sq = V.CODE[np.frombuffer(fetch(c, s, min(e, s + kb * 1000)), dtype=np.uint8)]; sh = V.dinuc_shuffle(sq, NR)
                    r = V.rec_at(sq, lag) - V.rec_at(sh, lag)
                    chf, n2, npos = V.chain_at(sq, lag, count=True); ch = chf - V.chain_at(sh, lag)
                    stat[kb] = (r, ch, n2, npos, len(sq))
                r, ch, n2, npos, L = stat[a.test_kb]
                tier = tier_of(ch, r, a.chain_cut, a.rec_cut)
                asis, cov = V.trf_support(trf, c, s, e, lag); tr = cov >= a.trf_cover
                tef = M58.overlap_bp(te, c, s, e) / (e - s)
                k = (cls, tier, 'TE' if tef > 0.5 else 'nonTE'); tot[k] += out; n[k] += 1; trf_ok[k] += tr; trf_as[k] += asis
                bed.write(f'{c}\t{s}\t{e}\t{tier}|unit{lag}|{cls}|TE{tef:.2f}|TRF{int(tr)}|TRFcov{cov:.2f}\t0\t.\t{s}\t{e}\t{COL[tier]}\n')
                w.writerow([c, s, e, tier, lag, cls, out, f'{tef:.3f}', int(asis), f'{cov:.3f}', int(tr), f'{r:.4f}', f'{ch:.4f}', n2, npos, L])
                if not a.no_sweep:
                    for kb in SWEEP_KB:
                        r_, ch_ = stat[kb][:2]
                        for cc in SWEEP_CUTS:
                            for rc in SWEEP_CUTS:
                                x = sweep[(kb, cc, rc, cls, tier_of(ch_, r_, cc, rc))]; x[0] += 1; x[1] += out; x[2] += tr
    if not a.no_sweep:
        with open(f'{CAL}/58b_tier_sweep.tsv', 'w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t', lineterminator='\n')
            w.writerow(['test_kb', 'chain_cut', 'rec_cut', 'unit_class', 'tier', 'n', 'Mb_outside_02', 'trf_cover_supported_share'])
            for key in sorted(sweep): w.writerow(list(key) + [sweep[key][0], f'{sweep[key][1] / 1e6:.3f}', f'{sweep[key][2] / sweep[key][0]:.3f}'])
    print(f'{TAG}: class\ttier\tTE\tn\tMb\tTRF cover>={a.trf_cover}\tTRF as coded')
    for k in sorted(tot): print(f'{TAG}\t' + '\t'.join(k) + f'\t{n[k]}\t{tot[k] / 1e6:.2f}\t{trf_ok[k] / n[k]:.2f}\t{trf_as[k] / n[k]:.2f}')
    arr = sum(v for k, v in tot.items() if k[1] == 'new_array'); arr_nte = sum(v for k, v in tot.items() if k[1] == 'new_array' and k[2] == 'nonTE')
    print(f'{TAG}: new arrays (>= 3 copies) {arr / 1e6:.2f} Mb, of which non-TE {arr_nte / 1e6:.2f} Mb')


if __name__ == '__main__':
    main()
