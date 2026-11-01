#!/usr/bin/env python3
"""
62_hor_detect.py  (SatChain v4.7 draft; higher-order repeat detection on the copy chain, two independent signals)

Per catalogue run (>= --min-copies copies with identity >= --min-identity, default 80 = the catalogue's confident class, complete copies
only unless --include-truncated; 2026-11-01, D07 / D14; before: identity >= 60, i.e. junk tiles entered), in genomic order:
1 Copy-chain lag identity. For lags n = 1..max_lag, identity between copy i and copy i+n (global alignment score proxy:
  match 1, mismatch 0, gap -1, divided by the longer length), median over <= --pairs sampled i; the pair count per lag is written.
  max_lag = min(--max-lag, copies // 3) per run (D16: lag n rests on copies - n pairs). A HOR of n monomers gives a peak at lag n and
  its multiples above lag 1.
  Null (2026-11-01, D03): under a random copy order every lag draws its pairs from the same pool of pair identities, so the null is built
  from a pool of <= --pool random copy pairs: --perm draws (default 1000) of max_lag medians of --pairs identities each give (a) the
  per-lag null p99 of a median and (b) the MAX statistic max_n [id(n) - id(1)] over all lags; p_max = share of draws whose max excess >=
  the observed best excess; q_max = Benjamini-Hochberg over runs.
  HOR call: n* = smallest lag 2 <= n <= max_lag // 3 with excess(n) = id(n) - id(1) >= --min-excess, id(n) above the per-lag null p99,
  AND id(2n) AND id(3n) above the null p99 with excess >= --min-excess / 2 (periodicity; before 2026-11-01 the 2n check was skipped when
  2n > max_lag and there was no 3n check), AND q_max <= --fdr. hor_monomers is written only for such runs; best_lag / id_best stay for all.
  Strand (D07): reverse_fraction and G3 are copied from the catalogue; strand_flag = 'mixed_strand' when 0.2 < reverse_fraction < 0.8
  (a strand flip inside the run; 03's `sequence` is oriented per copy, so a flip cannot create a lag peak, but the run is flagged).
2 k-mer interval spectrum (NTRprism idea). For k in --ks, the distance between consecutive occurrences of each canonical k-mer in the
  run sequence; histogram of distances 1..--max-int bp, normalised. Peaks reported: top 3 distances with their share. With growing k,
  monomer-wide k-mers drop out and position-specific ones remain, so a HOR length should take over from the monomer length.
Sweep (--sweep, D09): min_identity 60 / 65 / 80 x pairs 75 / 150 x min_excess 0.02 / 0.03 / 0.05 x perm 200 / 1000 -> runs tested, HOR
calls (raw and after BH) per setting, written to <repo>/benchmarks/62_sweep.tsv (or --sweep-out).
Outputs (<dest>/): 62_lag_identity_<tag>.tsv (run, lag, pairs, id, null_p50, null_p99), 62_kmer_spectrum_peaks_<tag>.tsv, 62_hor_calls_<tag>.tsv.
usage: 62_hor_detect.py --catalogue <59 catalogue> --copies <03 copies> [...] --genome <fa> --dest <dir> --tag <G> [--runs id,...]
       [--min-copies 40] [--max-lag 40] [--pairs 150] [--perm 1000] [--pool 2000] [--min-excess 0.03] [--fdr 0.05] [--min-identity 80]
       [--include-truncated] [--ks 6,10,15,21,31] [--max-int 20000] [--sweep] [--sweep-out path]
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
    """median identity per lag 1..max_lag and the number of pairs used per lag"""
    out, npairs = [], []
    for n in range(1, max_lag + 1):
        idx = list(range(len(seqs) - n))
        if not idx: break
        smp = rng.sample(idx, min(pairs, len(idx)))
        out.append(float(np.median([ident(seqs[i], seqs[i + n]) for i in smp]))); npairs.append(len(smp))
    return out, npairs


def pair_pool(seqs, pool, rng):
    """identities of <= pool random copy pairs (i < j): the exchangeable null of a shuffled copy order"""
    n = len(seqs); allp = n * (n - 1) // 2
    if allp <= pool: prs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    else: prs = {(min(i, j), max(i, j)) for i, j in zip(rng.choices(range(n), k=pool * 2), rng.choices(range(n), k=pool * 2)) if i != j}; prs = list(prs)[:pool]
    return np.array([ident(seqs[i], seqs[j]) for i, j in prs])


def null_draws(pool_ids, npairs, perm, seed):
    """perm draws of a lag profile under the null: medians of npairs[n] identities from the pool per lag.
    Returns per-lag (p50, p99) of the median and the null distribution of max_n [id(n) - id(1)]."""
    g = np.random.default_rng(seed); L = len(npairs); meds = np.empty((perm, L))
    for n in range(L): meds[:, n] = np.median(g.choice(pool_ids, size=(perm, npairs[n]), replace=True), axis=1)
    p50, p99 = np.percentile(meds, 50, axis=0), np.percentile(meds, 99, axis=0)
    mx = (meds[:, 1:] - meds[:, [0]]).max(1) if L > 1 else np.zeros(perm)
    return p50, p99, mx


def call_hor(prof, p99, min_excess, max_lag):
    """smallest n with excess >= min_excess, id(n), id(2n), id(3n) > p99 and n <= max_lag // 3 (D03); None when no lag qualifies"""
    for n in range(2, max_lag // 3 + 1):
        if n > len(prof) or 3 * n > len(prof): break
        ok = prof[n - 1] - prof[0] >= min_excess and prof[n - 1] > p99[n - 1]
        for m in (2 * n, 3 * n):
            ok = ok and prof[m - 1] > p99[m - 1] and prof[m - 1] - prof[0] >= min_excess / 2
        if ok: return n
    return None


def bh(p):
    """Benjamini-Hochberg q-values"""
    p = np.asarray(p, float); n = len(p)
    if n == 0: return p
    o = np.argsort(p); q = np.empty(n); prev = 1.0
    for r in range(n - 1, -1, -1):
        prev = min(prev, p[o[r]] * n / (r + 1)); q[o[r]] = prev
    return q


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
    ap.add_argument('--perm', type=int, default=1000); ap.add_argument('--pool', type=int, default=2000); ap.add_argument('--min-excess', type=float, default=0.03)
    ap.add_argument('--fdr', type=float, default=0.05); ap.add_argument('--min-identity', type=float, default=80); ap.add_argument('--include-truncated', action='store_true')
    ap.add_argument('--ks', default='6,10,15,21,31'); ap.add_argument('--max-int', type=int, default=20000); ap.add_argument('--seed', type=int, default=62)
    ap.add_argument('--sweep', action='store_true'); ap.add_argument('--sweep-out', default='')
    a = ap.parse_args(); rng = random.Random(a.seed); os.makedirs(a.dest, exist_ok=True)
    KS = [int(x) for x in a.ks.split(',')]; want = set(a.runs.split(',')) if a.runs else None
    runs = [r for r in csv.DictReader(open(a.catalogue), delimiter='\t') if r.get('element_class', 'satellite') == 'satellite'
            and int(r['n_conf']) >= a.min_copies and (want is None or r['run_id'] in want)]
    byc = collections.defaultdict(list)
    for r in runs: byc[r['chrom']].append((int(r['start0']), int(r['end']), r['run_id']))
    allc = collections.defaultdict(list)                 # run -> (start, identity, complete, sequence); filtered per setting
    for p in a.copies:
        with open(p, newline='') as fh:
            rd = csv.reader(fh, delimiter='\t'); h = next(rd); ix = {k: i for i, k in enumerate(h)}
            for row in rd:
                c = row[ix['chrom']]
                if c not in byc: continue
                s, e = int(row[ix['genomic_start_0based']]), int(row[ix['genomic_end']])
                for rs, re_, rid in byc[c]:
                    if s >= rs and e <= re_:
                        comp = row[ix['complete_or_truncated']] == 'complete' if 'complete_or_truncated' in ix else True
                        allc[rid].append((s, float(row[ix['best_alignment_identity']] or 0), comp, row[ix['sequence']].upper())); break
    def copies_of(rid, min_id, incl_trunc):
        return [x for _, i, comp, x in sorted(allc.get(rid, [])) if i >= min_id and (comp or incl_trunc)]

    def analyse(cs, pairs, perm, seed):
        """profile, pairs per lag, per-lag null (p50, p99), null max statistic, observed max excess"""
        r_ = random.Random(seed); max_lag = min(a.max_lag, len(cs) // 3)
        prof, npairs = lag_profile(cs, max_lag, pairs, r_)
        pool = pair_pool(cs, a.pool, r_); p50, p99, mx = null_draws(pool, npairs, perm, seed)
        obs = max(prof[1:]) - prof[0] if len(prof) > 1 else 0.0
        p_max = (1 + int((mx >= obs).sum())) / (1 + perm)
        return prof, npairs, p50, p99, mx, obs, p_max, max_lag

    fa = Fasta(a.genome); L1, P1, C1 = [], [], []; res = {}
    for r in runs:
        rid = r['run_id']; cs = copies_of(rid, a.min_identity, a.include_truncated)
        if len(cs) < a.min_copies: continue
        prof, npairs, p50, p99, mx, obs, p_max, max_lag = analyse(cs, a.pairs, a.perm, a.seed)
        res[rid] = (prof, p99, p_max, max_lag)
        for n, v in enumerate(prof, 1): L1.append(dict(run_id=rid, lag=n, pairs=npairs[n - 1], identity=round(v, 4), null_p50=round(float(p50[n - 1]), 4), null_p99=round(float(p99[n - 1]), 4)))
        seq = fa.get(r['chrom'], int(r['start0']), int(r['end']))
        peaks = {}
        for k in KS:
            sp = spectrum(seq, k, a.max_int); top = sorted(sp.items(), key=lambda x: -x[1])[:3]
            peaks[k] = top
            for rank, (d, v) in enumerate(top, 1): P1.append(dict(run_id=rid, k=k, rank=rank, interval_bp=d, share=round(v, 4)))
        rf = float(r.get('reverse_fraction') or 0)
        C1.append(dict(run_id=rid, chrom=r['chrom'], start0=r['start0'], end=r['end'], family=r.get('name_v46') or r['family'], unit_bp=r['unit_bp'],
                       copies=len(cs), max_lag=max_lag, id_lag1=round(prof[0], 3), best_lag=int(np.argmax(prof)) + 1, id_best=round(max(prof), 3),
                       best_excess=round(obs, 4), null_max_p99=round(float(np.percentile(mx, 99)), 4), p_max=round(p_max, 4), q_max='',
                       reverse_fraction=rf, G3=r.get('G3', ''), strand_flag='mixed_strand' if 0.2 < rf < 0.8 else '',
                       hor_candidate=call_hor(prof, p99, a.min_excess, max_lag) or '', hor_monomers='', hor_bp_estimate='',
                       **{f'k{k}_top_interval': peaks[k][0][0] if peaks[k] else '' for k in KS}))
        print(rid, 'lag1 %.3f best lag %d (%.3f) p_max %.3f candidate %s | top interval by k: %s' % (prof[0], C1[-1]['best_lag'], max(prof), p_max,
              C1[-1]['hor_candidate'] or '-', ', '.join(f'k{k}:{peaks[k][0][0]}' for k in KS if peaks[k])), flush=True)
    if C1:                                             # BH over runs; the call needs the lag rule AND q_max <= fdr
        q = bh([c['p_max'] for c in C1])
        for c, qv in zip(C1, q):
            c['q_max'] = round(float(qv), 4)
            if c['hor_candidate'] != '' and qv <= a.fdr:
                c['hor_monomers'] = c['hor_candidate']; c['hor_bp_estimate'] = int(c['hor_candidate']) * int(c['unit_bp'])
        print(f'runs {len(C1)}, HOR calls {sum(1 for c in C1 if c["hor_monomers"] != "")} (candidates before FDR {sum(1 for c in C1 if c["hor_candidate"] != "")}), '
              f'mixed-strand flags {sum(1 for c in C1 if c["strand_flag"])}')
    for name, rows in (('lag_identity', L1), ('kmer_spectrum_peaks', P1), ('hor_calls', C1)):
        if not rows: continue
        with open(os.path.join(a.dest, f'62_{name}_{a.tag}.tsv'), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
    if a.sweep:
        SW = []
        for mi in (60, 65, 80):
            for pairs in (75, 150):
                per = {}
                for r in runs:
                    cs = copies_of(r['run_id'], mi, a.include_truncated)
                    if len(cs) >= a.min_copies: per[r['run_id']] = (cs, analyse(cs, pairs, 1000, a.seed))
                for perm in (200, 1000):
                    for me in (0.02, 0.03, 0.05):
                        cand, pm = {}, []
                        for rid, (cs, (prof, npairs, p50, p99, mx, obs, p_max, max_lag)) in per.items():
                            mxp = mx[:perm]; pm.append((1 + int((mxp >= obs).sum())) / (1 + perm))
                            cand[rid] = call_hor(prof, p99, me, max_lag)          # per-lag p99 from the 1000 draws; perm only subsets the max statistic
                        q = bh(pm) if pm else []
                        calls = sum(1 for (rid, c), qv in zip(cand.items(), q) if c and qv <= a.fdr)
                        SW.append(dict(step='62', tag=a.tag, min_identity=mi, pairs=pairs, perm=perm, min_excess=me, runs=len(per),
                                       candidates=sum(1 for c in cand.values() if c), hor_calls_fdr=calls))
                        print(SW[-1], flush=True)
        out = a.sweep_out or os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'benchmarks', '62_sweep.tsv')
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(SW[0].keys()), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(SW)
        print('sweep ->', out)


if __name__ == '__main__':
    main()
