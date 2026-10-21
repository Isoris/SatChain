#!/usr/bin/env python3
"""Per-run sequence check of the 56 HMM runs: mean 1-kb tandem signal over the run vs the 'neither' background quantile (p80 / p90 / p95).
Writes 56_runs_validated_q<q>.bed and prints runs / Mb kept and the bp overlap with step 36. usage: val56_runfilter.py <tag> <signal> <56 bed> <36 bed> <02 tsv> <dest>"""
import sys, bisect, collections, csv
TAG, SIG, H, R36, A02, DEST = sys.argv[1:7]
sig = collections.defaultdict(list)
for l in open(SIG):
    if l.startswith('track'): continue
    c, s, e, v = l.split()[:4]; sig[c].append((int(s), int(e), float(v)))
for v in sig.values(): v.sort()
def bed(p):
    d = collections.defaultdict(list)
    for l in open(p):
        if l.startswith(('#', 'track')): continue
        f = l.rstrip('\n').split('\t'); d[f[0]].append((int(f[1]), int(f[2]), f[3]))
    for v in d.values(): v.sort()
    return d
H_, R = bed(H), bed(R36)
a02 = collections.defaultdict(list)
for r in csv.DictReader(open(A02), delimiter='\t'): a02[r['chrom']].append((int(r['start_0based']), int(r['end'])))
def cov(d, c, s, e):
    v = d.get(c, []); i = bisect.bisect_left(v, (s - 2000000, 0, '')); t = 0
    for x in v[i:]:
        if x[0] >= e: break
        t += max(0, min(e, x[1]) - max(s, x[0]))
    return t
def run_sig(c, s, e):
    v = sig.get(c, []); i = bisect.bisect_left(v, (s - 1000, 0, 0)); vals = []
    for x in v[i:]:
        if x[0] >= e: break
        if x[1] > s: vals.append(x[2])
    return sum(vals) / len(vals) if vals else 0.0
# background: 1-kb bins inside 02 intervals that are in no run
inrun = lambda c, m: cov(H_, c, m, m + 1) or cov(R, c, m, m + 1)
bg = []
for c, v in sig.items():
    A = sorted(a02.get(c, []))
    for s, e, x in v[::3]:
        m = (s + e) // 2
        if any(a <= m < b for a, b in A) and not inrun(c, m): bg.append(x)
bg.sort()
for q in (0.80, 0.90, 0.95):
    thr = bg[int(q * len(bg))]; kept = []; lost = 0
    for c, v in H_.items():
        for s, e, nm in v:
            if run_sig(c, s, e) > thr: kept.append((c, s, e, nm))
            else: lost += e - s
    bp = sum(e - s for c, s, e, _ in kept); sh = sum(cov(R, c, s, e) for c, s, e, _ in kept)
    tot36 = sum(e - s for v in R.values() for s, e, _ in v)
    with open(f'{DEST}/56_runs_validated_q{q}.bed', 'w') as fh:
        for c, s, e, nm in kept: fh.write(f'{c}\t{s}\t{e}\t{nm}\t0\t.\n')
    print(f'{TAG}\tbackground q{q}={thr:.4f}\truns kept {len(kept)}/{sum(len(v) for v in H_.values())}\tkept {bp / 1e6:.2f} Mb (dropped {lost / 1e6:.2f})\t'
          f'shared with 36 {sh / 1e6:.2f} Mb = {sh / tot36:.3f} of step-36 bp\tadded beyond 36 {(bp - sh) / 1e6:.2f} Mb')
