#!/usr/bin/env python3
"""val56: sequence support of the copy-chain HMM runs (56) vs the step-36 runs, from the 1-kb tandem-signal track (step-48 k-mer
recurrence excess at the array unit; 0 = no tandem structure). Classes of 1-kb bins inside step-02 intervals:
  shared (in a 36 run and an HMM run) | hmm_only | step36_only | neither (inside a 02 interval, in no run)
  plus outside_all (2026-11-01, CA34 / HB08: bins outside every 02 interval and every run = the val56_runfilter background)
Reported per class: n bins, mean / median signal, fraction of bins with signal > 0 and > the 'neither' and 'outside_all' 90th percentiles;
the p80 / p90 / p95 of both backgrounds and their share of exactly-zero bins. Every <subsample>-th bin is used (default 1 = all bins; pass
the same value to val56_runfilter --subsample so both scripts see one background).
usage: val56_signal.py <tag> <tandem_signal bedgraph> <56 runs bed> <36 bed> <02 arrays tsv> [subsample]"""
import sys, bisect, collections, csv, statistics as st
TAG, SIG, H, R36, A02 = sys.argv[1:6]; SUB = int(sys.argv[6]) if len(sys.argv) > 6 else 1
def bed(p):
    d = collections.defaultdict(list)
    for l in open(p):
        if l.startswith(('#', 'track')): continue
        f = l.split('\t'); d[f[0]].append((int(f[1]), int(f[2])))
    for v in d.values(): v.sort()
    return d
def inside(d, c, x):
    v = d.get(c, []); i = bisect.bisect_right(v, (x, 1 << 62)) - 1
    return any(v[j][0] <= x < v[j][1] for j in range(max(0, i - 2), i + 1)) if i >= 0 else False
h, r36 = bed(H), bed(R36)
a02 = collections.defaultdict(list)
for r in csv.DictReader(open(A02), delimiter='\t'): a02[r['chrom']].append((int(r['start_0based']), int(r['end'])))
for v in a02.values(): v.sort()
cls = collections.defaultdict(list)
for k_, l in enumerate(x for x in open(SIG) if not x.startswith('track')):
    if k_ % SUB: continue
    c, s, e, v = l.split()[:4]; m = (int(s) + int(e)) // 2
    a, b = inside(r36, c, m), inside(h, c, m)
    if not inside(a02, c, m):
        if not a and not b: cls['outside_all'].append(float(v))
        continue
    cls['shared' if a and b else 'step36_only' if a else 'hmm_only' if b else 'neither'].append(float(v))
qt = lambda v, q: sorted(v)[min(len(v) - 1, int(q * len(v)))] if v else 0
p90 = qt(cls['neither'], 0.9); o90 = qt(cls['outside_all'], 0.9)
for k in ('neither', 'outside_all'):
    v = cls[k]; print(f'{TAG}: {k} background (subsample {SUB}): bins {len(v)}, zero share {sum(x == 0 for x in v) / len(v) if v else float("nan"):.3f}, '
                      f'p80 {qt(v, 0.8):.4f} p90 {qt(v, 0.9):.4f} p95 {qt(v, 0.95):.4f}')
for k in ('shared', 'hmm_only', 'step36_only', 'neither', 'outside_all'):
    v = cls[k]
    if not v: print(TAG, k, 'none'); continue
    print(f'{TAG}\t{k:12s}\tbins {len(v):7d}\tmean {st.mean(v):.4f}\tmedian {st.median(v):.4f}\t>0 {sum(x > 0 for x in v) / len(v):.2f}\t'
          f'>neither_p90 {sum(x > p90 for x in v) / len(v):.2f}\t>outside_p90 {sum(x > o90 for x in v) / len(v):.2f}')
