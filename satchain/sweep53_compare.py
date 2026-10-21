#!/usr/bin/env python3
"""Compare step-53 runs of a one-at-a-time parameter sweep with the default run.
Per setting: accepted boundaries (n), recall / precision of the default's boundaries (match within max(5 kb, the array's
default window)), bp-weighted agreement of segment profiles with the default (1-kb bins), class-merge decisions that differ.
usage: sweep53_compare.py <sweep dir with one sub-folder per setting, 'default' among them> <out tsv>"""
import csv, collections, os, sys
D, OUT = sys.argv[1:3]
rd = lambda p: list(csv.DictReader(open(p), delimiter='\t')) if os.path.exists(p) else []


def load(run):
    p = lambda f: os.path.join(D, run, f)
    b = collections.defaultdict(list)
    for r in rd(p('53_change_points.tsv')):
        if r['status'] == 'accepted': b[r['array_id']].append(int(r['position']))
    win = {r['array_id']: int(r['window_bp']) for r in rd(p('53_array_summary.tsv'))}
    bins = {}
    for r in rd(p('53_segments.tsv')):
        for x in range(int(r['start0']) // 1000, int(r['end']) // 1000): bins[(r['array_id'], x)] = r['profile']
    merges = {(r['array_id'], r['class_a'], r['class_b']): r['decision'] for r in rd(p('53_class_merges.tsv')) if r['class_b'] != 'NA'}
    prof = collections.Counter(bins.values())
    return b, win, bins, merges, prof


ref = load('default'); rows = []
for run in sorted(os.listdir(D)):
    if not os.path.isdir(os.path.join(D, run)): continue
    b, win, bins, merges, prof = load(run)
    hit = tot = 0; ph = pt = 0
    for aid, xs in ref[0].items():
        tol = max(5000, ref[1].get(aid, 0))
        for x in xs: tot += 1; hit += any(abs(x - y) <= tol for y in b.get(aid, []))
    for aid, ys in b.items():
        tol = max(5000, ref[1].get(aid, 0))
        for y in ys: pt += 1; ph += any(abs(x - y) <= tol for x in ref[0].get(aid, []))
    common = set(bins) & set(ref[2]); agree = sum(bins[k] == ref[2][k] for k in common)
    md = sum(1 for k, v in merges.items() if ref[3].get(k) not in (None, v))
    rows.append([run, sum(map(len, b.values())), round(hit / tot, 3) if tot else 'NA', round(ph / pt, 3) if pt else 'NA',
                 round(agree / len(common), 3) if common else 'NA', md, len(merges)] +
                [round(prof[p] / max(1, sum(prof.values())), 3) for p in ('homogeneous_monomeric', 'variant_mixed', 'alternating', 'interleaved_families', 'unresolved')])
with open(OUT, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['setting', 'n_accepted_boundaries', 'recall_of_default', 'precision_vs_default', 'profile_agreement_bp', 'merge_decisions_changed', 'n_pairs_tested',
                'frac_homogeneous', 'frac_variant_mixed', 'frac_alternating', 'frac_interleaved', 'frac_unresolved'])
    w.writerows(rows)
for r in rows: print('\t'.join(map(str, r)))
