#!/usr/bin/env python3
"""val56_runfilter.py  (2026-11-01 rewrite: HB08 / CA34 / CA47 / HB19) per-run sequence check of copy-chain runs against a background null.
RUNS (--runs): a 56 BED (name = array|run|n_conf|n_tiles|mean_post) or a catalogue / runs-of-record TSV (run_id, chrom, start0, end
  [, unit_bp]); one interface for the existing (56) and the new-source (56new / 58b) runs, which the pre-2026-11-01 version never saw.
SCORE per run = bp-weighted mean over the track bins overlapping the run (weight = overlap bp, so a flanking bin counts by its overlap):
  --scan57 given : step-57 max_excess at ANY lag (windows of the chosen --scan57-scale), i.e. the run's own unit is scored whatever the 02
                   period; lag_matches_unit = the bp-weighted modal best_lag within 6 % of unit_bp, 2 x or 0.5 x unit_bp (needs unit_bp)
  --signal given : the step-48 1-kb tandem signal (k-mer recurrence excess at the 02 ARRAY period): a run whose true unit differs from its
                   02 period can fail for the wrong reason (48_summary unit_mismatch), so --scan57 is preferred when it exists.
  n_bins is written per run; a run with fewer than --min-bins overlapping bins is flagged few_bins (decided on its weighted mean anyway).
BACKGROUND: the same track, bins whose midpoint lies OUTSIDE ALL candidate intervals (step-02 intervals, step-36 runs, the scored runs,
  every --exclude BED: other HMM run sets, new candidates), every --subsample-th bin (default 1 = all bins, as val56_signal). The
  pre-2026-11-01 null (inside 02 intervals but in no run) held uncalled tandem sequence (CA34). Threshold = quantile q of the background;
  the share of exactly-zero background bins is reported because it decides how trivial a low q is (HB08).
OUTPUTS (--dest): 56_runs_validated_q<q>.bed per q, val56_runs_<tag>.tsv (score, n_bins, lag, flags, kept_q<q>), val56_benchmark_<tag>.tsv
  (q, threshold, runs kept / total, Mb kept / dropped, Mb shared with 36, fraction of 36 bp kept, Mb added, background bins, zero share):
  the table every README / METHODS val56 number must be read from (HB19).
usage: val56_runfilter.py --tag G --runs <bed|tsv> --bed36 <36 bed> --a02 <02 tsv> --dest <dir> (--signal <48 bedgraph> | --scan57 <57 tsv>)
       [--exclude <bed> ...] [--q 0.8,0.9,0.95] [--min-bins 2] [--subsample 1] [--scan57-scale small|large|both]"""
import argparse, bisect, collections, csv, os, sys


def bed(p, name=True):
    d = collections.defaultdict(list)
    for l in open(p):
        if l.startswith(('#', 'track')) or not l.strip(): continue
        f = l.rstrip('\n').split('\t'); d[f[0]].append((int(f[1]), int(f[2]), f[3] if name and len(f) > 3 else f'{f[0]}:{f[1]}-{f[2]}'))
    for v in d.values(): v.sort()
    return d


def read_runs(p):
    """56 BED or runs TSV -> {chrom: [(s, e, id, unit_bp|None)]}"""
    head = open(p).readline()
    d = collections.defaultdict(list)
    if 'run_id' in head.split('\t'):
        for r in csv.DictReader(open(p), delimiter='\t'):
            u = r.get('unit_bp'); d[r['chrom']].append((int(r['start0']), int(r['end']), r['run_id'], float(u) if u not in (None, '', 'NA') else None))
    else:
        for c, v in bed(p).items():
            for s, e, nm in v: d[c].append((s, e, nm, None))
    for v in d.values(): v.sort()
    return d


def read_track(a):
    """{chrom: sorted [(s, e, value, lag)]}, maxlen"""
    t = collections.defaultdict(list); ml = 0
    if a.scan57:
        for r in csv.DictReader(open(a.scan57), delimiter='\t'):
            if a.scan57_scale != 'both' and r['scale'] != a.scan57_scale: continue
            if r['max_excess'] == 'NA' or r['best_lag'] == 'NA': continue                       # N-rich windows (57 writes NA)
            s, e = int(r['start']), int(r['end']); t[r['chrom']].append((s, e, float(r['max_excess']), int(float(r['best_lag'])))); ml = max(ml, e - s)
    else:
        for l in open(a.signal):
            if l.startswith('track') or not l.strip(): continue
            c, s, e, v = l.split()[:4]; t[c].append((int(s), int(e), float(v), None)); ml = max(ml, int(e) - int(s))
    for v in t.values(): v.sort()
    return t, ml


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tag', required=True); ap.add_argument('--runs', required=True); ap.add_argument('--bed36', required=True)
    ap.add_argument('--a02', required=True); ap.add_argument('--dest', required=True)
    ap.add_argument('--signal', help='step-48 1-kb tandem-signal bedgraph (at the 02 period)'); ap.add_argument('--scan57', help='step-57 scan tsv (excess at any lag)')
    ap.add_argument('--scan57-scale', choices=['small', 'large', 'both'], default='both')
    ap.add_argument('--exclude', nargs='*', default=[], help='more candidate-interval BEDs kept out of the background')
    ap.add_argument('--q', default='0.80,0.90,0.95'); ap.add_argument('--min-bins', type=int, default=2); ap.add_argument('--subsample', type=int, default=1); ap.add_argument('--allow-untested', action='store_true', help='keep (flagged) the runs of chromosomes the track does not cover instead of stopping')
    a = ap.parse_args()
    if bool(a.signal) == bool(a.scan57): ap.error('give exactly one of --signal (48 track) or --scan57 (57 scan)')
    os.makedirs(a.dest, exist_ok=True); Q = [float(x) for x in a.q.split(',')]
    runs = read_runs(a.runs); R36 = bed(a.bed36); track, ML = read_track(a)
    cand = collections.defaultdict(list)                                   # every candidate interval: 02, 36, runs, --exclude
    for r in csv.DictReader(open(a.a02), delimiter='\t'): cand[r['chrom']].append((int(r['start_0based']), int(r['end'])))
    for d in (R36, runs) + tuple(bed(p) for p in a.exclude):
        for c, v in d.items(): cand[c] += [(x[0], x[1]) for x in v]
    for c in cand:                                                         # merge
        m = []
        for s, e in sorted(cand[c]):
            if m and s <= m[-1][1]: m[-1] = (m[-1][0], max(m[-1][1], e))
            else: m.append((s, e))
        cand[c] = m
    def inside(c, x):
        v = cand.get(c, []); i = bisect.bisect_right(v, (x, 1 << 62)) - 1
        return i >= 0 and v[i][0] <= x < v[i][1]
    def cov(d, c, s, e):
        v = d.get(c, []); i = bisect.bisect_left(v, (s - 5_000_000, 0, '')); t = 0
        for x in v[i:]:
            if x[0] >= e: break
            t += max(0, min(e, x[1]) - max(s, x[0]))
        return t
    def score(c, s, e):
        v = track.get(c, []); i = bisect.bisect_left(v, (s - ML, 0, 0.0, None)); w = 0; acc = 0.0; n = 0; lags = collections.Counter()
        for x in v[i:]:
            if x[0] >= e: break
            o = min(e, x[1]) - max(s, x[0])
            if o <= 0: continue
            acc += o * x[2]; w += o; n += 1
            if x[3] is not None: lags[x[3]] += o
        return (acc / w if w else 0.0), n, (lags.most_common(1)[0][0] if lags else None)
    bg = []; nbg_all = 0
    for c, v in track.items():
        for s, e, x, _ in v[::a.subsample]:
            nbg_all += 1
            if not inside(c, (s + e) // 2): bg.append(x)
    bg.sort()
    if not bg: sys.exit('ERROR no background bins outside the candidate intervals')
    zero = sum(1 for x in bg if x == 0) / len(bg)
    missing = sorted(c for c in runs if c not in track)                   # 2026-11-01: a chromosome absent from the track must never drop its runs silently
    if missing and not a.allow_untested: sys.exit(f'ERROR no track bins on {len(missing)} chromosome(s) carrying runs: {missing}; give a complete track or --allow-untested (runs kept, flagged untested)')
    if missing: print(f'WARNING {len(missing)} chromosome(s) without track bins, their runs are kept as untested: {missing}', file=sys.stderr)
    rows = []
    for c, v in runs.items():
        for s, e, nm, u in v:
            sc, n, lag = score(c, s, e)
            if c in missing: sc = float('inf')                            # untested: kept at every q
            lm = 'NA' if lag is None or not u else ('yes' if any(abs(lag - k * u) <= 0.06 * k * u for k in (1, 2, 0.5)) else 'no')
            rows.append(dict(run_id=nm, chrom=c, start0=s, end=e, bp=e - s, score=round(sc, 4) if sc != float('inf') else 'untested', n_bins=n, few_bins='yes' if n < a.min_bins else 'no', untested='yes' if c in missing else 'no',
                             modal_lag=lag if lag is not None else 'NA', unit_bp=u if u else 'NA', lag_matches_unit=lm, shared36_bp=cov(R36, c, s, e)))
    tot36 = sum(e - s for v in R36.values() for s, e, _ in v); ntot = len(rows); bench = []
    for q in Q:
        thr = bg[min(len(bg) - 1, int(q * len(bg)))]; ok = lambda r: r['untested'] == 'yes' or r['score'] > thr; kept = [r for r in rows if ok(r)]
        for r in rows: r[f'kept_q{q}'] = 'yes' if ok(r) else 'no'
        bp = sum(r['bp'] for r in kept); sh = sum(r['shared36_bp'] for r in kept); lost = sum(r['bp'] for r in rows) - bp
        with open(os.path.join(a.dest, f'56_runs_validated_q{q}.bed'), 'w') as fh:
            for r in kept: fh.write(f"{r['chrom']}\t{r['start0']}\t{r['end']}\t{r['run_id']}\t0\t.\n")
        bench.append(dict(tag=a.tag, track='scan57_any_lag' if a.scan57 else 'signal48_02period', q=q, threshold=round(thr, 4), runs_kept=len(kept), runs_total=ntot,
                          kept_Mb=round(bp / 1e6, 3), dropped_Mb=round(lost / 1e6, 3), shared_with_36_Mb=round(sh / 1e6, 3),
                          frac_36_bp_kept=round(sh / tot36, 4) if tot36 else 'NA', added_beyond_36_Mb=round((bp - sh) / 1e6, 3),
                          runs_few_bins=sum(1 for r in rows if r['few_bins'] == 'yes'), runs_untested=sum(1 for r in rows if r['untested'] == 'yes'), background_bins=len(bg), background_zero_share=round(zero, 3),
                          background_rule='outside_all_candidates', subsample=a.subsample, min_bins=a.min_bins))
        print(f"{a.tag}\tq{q} threshold={thr:.4f}\truns kept {len(kept)}/{ntot}\tkept {bp / 1e6:.2f} Mb (dropped {lost / 1e6:.2f})\t"
              f"shared with 36 {sh / 1e6:.2f} Mb = {sh / tot36 if tot36 else float('nan'):.3f} of step-36 bp\tadded beyond 36 {(bp - sh) / 1e6:.2f} Mb")
    print(f'{a.tag}\tbackground: {len(bg)} bins outside all candidate intervals (of {nbg_all} sampled), zero share {zero:.3f}; '
          f'runs on < {a.min_bins} bins: {sum(1 for r in rows if r["few_bins"] == "yes")}/{ntot}')
    with open(os.path.join(a.dest, f'val56_runs_{a.tag}.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
    with open(os.path.join(a.dest, f'val56_benchmark_{a.tag}.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(bench[0].keys()), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(bench)


if __name__ == '__main__':
    main()
