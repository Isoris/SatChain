#!/usr/bin/env python3
"""
62_hor_detect.py  (SatChain v4.7 draft; higher-order repeat detection on the copy chain, two independent signals)

Per catalogue run (>= --min-copies confident copies, in genomic order):
1 Copy-chain lag identity. For lags n = 1..--max-lag, identity between copy i and copy i+n (global alignment score proxy:
  match 1, mismatch 0, gap -1, divided by the longer length), median over <= --pairs sampled i. A HOR of n monomers gives a peak at
  lag n (and its multiples) above lag 1. Null: the same statistic after shuffling copy order within the run (--perm times).
  HOR call: n* = smallest lag >= 2 with excess(n) = id(n) - id(1) >= --min-excess, id(n) above the null p99, and id(2n) also above
  the null p99 when 2n <= max-lag (periodicity check).
2 k-mer interval spectrum (NTRprism idea). For k in --ks, the distance between consecutive occurrences of each canonical k-mer in the
  run sequence; histogram of distances 1..--max-int bp, normalised. Peaks reported: top 3 distances with their share. With growing k,
  monomer-wide k-mers drop out and position-specific ones remain, so a HOR length should take over from the monomer length.
Outputs (<dest>/): 62_lag_identity_<tag>.tsv (run, lag, id, null_p50, null_p99), 62_kmer_spectrum_peaks_<tag>.tsv, 62_hor_calls_<tag>.tsv.
usage: 62_hor_detect.py --catalogue <59 catalogue> --copies <03 copies> [...] --genome <fa> --dest <dir> --tag <G> [--runs id,...]
       [--min-copies 40] [--max-lag 40] [--pairs 150] [--perm 20] [--min-excess 0.03] [--ks 6,10,15,21,31] [--max-int 20000]
"""
import argparse, collections, csv, os, random
import numpy as np
from Bio import Align
csv.field_size_limit(1 << 30)
RC = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(RC)[::-1]
GA = Align.PairwiseAligner(mode='global', match_score=1, mismatch_score=0, open_gap_score=-1, extend_gap_score=-1)
def ident(a, b): return GA.score(a, b) / max(len(a), len(b), 1)


class Fasta:
    def __init__(self, path):
        self.fh = open(path, 'rb'); self.idx = {}
        for l in open(path + '.fai'):
            n, ln, off, lb, lby = l.split('\t')[:5]; self.idx[n] = (int(ln), int(off), int(lb), int(lby))
    def get(self, c, a, b):
        ln, off, lb, lby = self.idx[c]; a = max(0, a); b = min(ln, b)
        s0 = off + (a // lb) * lby + a % lb; s1 = off + (b // lb) * lby + b % lb
        self.fh.seek(s0); return self.fh.read(s1 - s0).replace(b'\n', b'').replace(b'\r', b'').decode().upper()


def lag_profile(seqs, max_lag, pairs, rng):
    out = []
    for n in range(1, max_lag + 1):
        idx = list(range(len(seqs) - n))
        if not idx: break
        smp = rng.sample(idx, min(pairs, len(idx)))
        out.append(float(np.median([ident(seqs[i], seqs[i + n]) for i in smp])))
    return out


def spectrum(seq, k, max_int):
    last = {}; h = collections.Counter(); n = len(seq)
    r = rc(seq)
    for i in range(n - k + 1):
        a = seq[i:i + k]
        if 'N' in a: continue
        b = r[n - i - k:n - i]; c = a if a < b else b
        if c in last:
            d = i - last[c]
            if d <= max_int: h[d] += 1
        last[c] = i
    tot = sum(h.values()) or 1
    return {d: v / tot for d, v in h.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--catalogue', required=True); ap.add_argument('--copies', nargs='+', required=True); ap.add_argument('--genome', required=True)
    ap.add_argument('--dest', required=True); ap.add_argument('--tag', required=True); ap.add_argument('--runs', default='')
    ap.add_argument('--min-copies', type=int, default=40); ap.add_argument('--max-lag', type=int, default=40); ap.add_argument('--pairs', type=int, default=150)
    ap.add_argument('--perm', type=int, default=20); ap.add_argument('--min-excess', type=float, default=0.03)
    ap.add_argument('--ks', default='6,10,15,21,31'); ap.add_argument('--max-int', type=int, default=20000); ap.add_argument('--seed', type=int, default=62)
    a = ap.parse_args(); rng = random.Random(a.seed); os.makedirs(a.dest, exist_ok=True)
    KS = [int(x) for x in a.ks.split(',')]; want = set(a.runs.split(',')) if a.runs else None
    runs = [r for r in csv.DictReader(open(a.catalogue), delimiter='\t') if r.get('element_class', 'satellite') == 'satellite'
            and int(r['n_conf']) >= a.min_copies and (want is None or r['run_id'] in want)]
    byc = collections.defaultdict(list)
    for r in runs: byc[r['chrom']].append((int(r['start0']), int(r['end']), r['run_id']))
    seqs = collections.defaultdict(list)
    for p in a.copies:
        with open(p, newline='') as fh:
            rd = csv.reader(fh, delimiter='\t'); h = next(rd); ix = {k: i for i, k in enumerate(h)}
            for row in rd:
                c = row[ix['chrom']]
                if c not in byc: continue
                s, e = int(row[ix['genomic_start_0based']]), int(row[ix['genomic_end']])
                for rs, re_, rid in byc[c]:
                    if s >= rs and e <= re_:
                        if float(row[ix['best_alignment_identity']] or 0) >= 60: seqs[rid].append((s, row[ix['sequence']].upper()))
                        break
    fa = Fasta(a.genome); L1, P1, C1 = [], [], []
    for r in runs:
        rid = r['run_id']; cs = [x for _, x in sorted(seqs.get(rid, []))]
        if len(cs) < a.min_copies: continue
        prof = lag_profile(cs, a.max_lag, a.pairs, rng)
        nulls = []
        for _ in range(a.perm):
            sh = cs[:]; rng.shuffle(sh); nulls.append(lag_profile(sh, min(5, a.max_lag), a.pairs // 2, rng))
        nv = np.array([v for n_ in nulls for v in n_]); p50, p99 = float(np.percentile(nv, 50)), float(np.percentile(nv, 99))
        for n, v in enumerate(prof, 1): L1.append(dict(run_id=rid, lag=n, identity=round(v, 4), null_p50=round(p50, 4), null_p99=round(p99, 4)))
        hor = ''
        for n in range(2, len(prof) + 1):
            ok = prof[n - 1] - prof[0] >= a.min_excess and prof[n - 1] > p99
            if ok and 2 * n <= len(prof): ok = prof[2 * n - 1] > p99 and prof[2 * n - 1] - prof[0] >= a.min_excess / 2
            if ok: hor = n; break
        seq = fa.get(r['chrom'], int(r['start0']), int(r['end']))
        peaks = {}
        for k in KS:
            sp = spectrum(seq, k, a.max_int); top = sorted(sp.items(), key=lambda x: -x[1])[:3]
            peaks[k] = top
            for rank, (d, v) in enumerate(top, 1): P1.append(dict(run_id=rid, k=k, rank=rank, interval_bp=d, share=round(v, 4)))
        C1.append(dict(run_id=rid, chrom=r['chrom'], start0=r['start0'], end=r['end'], family=r.get('name_v46') or r['family'], unit_bp=r['unit_bp'],
                       copies=len(cs), id_lag1=round(prof[0], 3), best_lag=int(np.argmax(prof)) + 1, id_best=round(max(prof), 3), null_p99=round(p99, 3),
                       hor_monomers=hor, hor_bp_estimate=(int(hor) * int(r['unit_bp']) if hor else ''),
                       **{f'k{k}_top_interval': peaks[k][0][0] if peaks[k] else '' for k in KS}))
        print(rid, 'lag1 %.3f best lag %d (%.3f) HOR %s | top interval by k: %s' % (prof[0], C1[-1]['best_lag'], max(prof), hor or '-',
              ', '.join(f'k{k}:{peaks[k][0][0]}' for k in KS if peaks[k])), flush=True)
    for name, rows in (('lag_identity', L1), ('kmer_spectrum_peaks', P1), ('hor_calls', C1)):
        if not rows: continue
        with open(os.path.join(a.dest, f'62_{name}_{a.tag}.tsv'), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(rows)


if __name__ == '__main__':
    main()
