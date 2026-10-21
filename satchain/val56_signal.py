#!/usr/bin/env python3
"""val56: sequence support of the copy-chain HMM runs (56) vs the step-36 runs, from the 1-kb tandem-signal track (step-48 k-mer
recurrence excess at the array unit; 0 = no tandem structure). Classes of 1-kb bins inside step-02 intervals:
  shared (in a 36 run and an HMM run) | hmm_only | step36_only | neither (inside a 02 interval, in no run)
Reported per class: n bins, mean / median signal, fraction of bins with signal > 0 and > the 'neither' 90th percentile.
usage: val56_signal.py <tag> <tandem_signal bedgraph> <56 runs bed> <36 bed> <02 arrays tsv>"""
import sys, bisect, collections, csv, statistics as st
TAG, SIG, H, R36, A02 = sys.argv[1:6]
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
for l in open(SIG):
    if l.startswith('track'): continue
    c, s, e, v = l.split()[:4]; m = (int(s) + int(e)) // 2
    if not inside(a02, c, m): continue
    a, b = inside(r36, c, m), inside(h, c, m)
    cls['shared' if a and b else 'step36_only' if a else 'hmm_only' if b else 'neither'].append(float(v))
p90 = sorted(cls['neither'])[int(0.9 * len(cls['neither']))] if cls['neither'] else 0
print(f'{TAG}: neither-class p90 signal = {p90:.4f}')
for k in ('shared', 'hmm_only', 'step36_only', 'neither'):
    v = cls[k]
    if not v: print(TAG, k, 'none'); continue
    print(f'{TAG}\t{k:12s}\tbins {len(v):7d}\tmean {st.mean(v):.4f}\tmedian {st.median(v):.4f}\t>0 {sum(x > 0 for x in v) / len(v):.2f}\t>neither_p90 {sum(x > p90 for x in v) / len(v):.2f}')
