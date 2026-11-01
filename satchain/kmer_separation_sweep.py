#!/usr/bin/env python3
"""Sensitivity sweep: at which k-mer length and diagnostic cutoff do two groups of copies separate?
For each pair of groups (step-03 classes, or the two domains either side of a step-53 boundary), counts k-mers present in
>= p_in of one group's copies and <= p_out of the other's, for k = 5..20 and four (p_in, p_out) cutoffs, and the same count
after shuffling the group labels (--perms permutations; group sizes kept). A k-mer is counted once per copy (presence).
Also reports the best single-k-mer classifier accuracy (balanced) per k.
2026-11-01 (CB13 / CA36): --perms default 50 (was 5: null_max_total was 0 in every benchmark row); the stderr line prints the
0.9 / 0.1 cut (it printed the 0.8 / 0.2 row); --k-report <k> prints, per comparison, the decision at that k (diagnostic k-mers at the
0.9 / 0.1 cut vs the null maximum and the METHODS rule '>= 5 near-fixed k-mers and above the null') and the k range where the rule
holds, and writes it to <out prefix>.k_report.tsv; the merge decision depends on k (benchmark Gar_SAT01307: 3-11 k-mers for k 5-14,
0-1 for k >= 18), so the k used must be stated.
usage: kmer_separation_sweep.py <03_monomer_copies.tsv> <array_id> <boundary_bp|-> <out prefix> [--perms 50] [--k-report 10]"""
import csv, collections, random, sys, argparse
csv.field_size_limit(1 << 30)
ap = argparse.ArgumentParser()
for x in ('copies', 'aid', 'bnd', 'out'): ap.add_argument(x)
ap.add_argument('--perms', type=int, default=50); ap.add_argument('--k-report', type=int, default=10)
ap.add_argument('--min-diag', type=int, default=5, help='METHODS rule: diagnostic k-mers needed at the 0.9 / 0.1 cut')
ARGS = ap.parse_args(); COPIES, AID, BND, OUT = ARGS.copies, ARGS.aid, ARGS.bnd, ARGS.out
KS = [5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20]
CUTS = [(0.9, 0.1), (0.8, 0.2), (0.7, 0.3), (0.6, 0.4)]
PERMS = ARGS.perms

rows = []
with open(COPIES, newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if r['array_id'] == AID and r['complete_or_truncated'] == 'complete' and float(r['best_alignment_identity'] or 0) >= 80:
            rows.append((r['monomer_class'], int(r['genomic_start_0based']), int(r['genomic_end']), r['sequence'].upper()))
cls = collections.Counter(c for c, *_ in rows)
top = [c for c, _ in cls.most_common(3)]
pairs = [(f'{a} vs {b}', [i for i, r in enumerate(rows) if r[0] == a], [i for i, r in enumerate(rows) if r[0] == b]) for i_, a in enumerate(top) for b in top[i_ + 1:]]
if BND != '-':
    b = int(BND)
    pairs.append(('left vs right of boundary', [i for i, r in enumerate(rows) if r[2] <= b], [i for i, r in enumerate(rows) if r[1] >= b]))


def counts(sets, idx):
    c = collections.Counter()
    for i in idx: c.update(sets[i])
    return c


def diag(ca, na, cb, nb, pin, pout):
    da = sum(1 for k, v in ca.items() if v >= pin * na and cb.get(k, 0) <= pout * nb)
    db = sum(1 for k, v in cb.items() if v >= pin * nb and ca.get(k, 0) <= pout * na)
    return da, db


def best_acc(ca, na, cb, nb):
    best = 0.5
    for k in set(ca) | set(cb):
        fa, fb = ca.get(k, 0) / na, cb.get(k, 0) / nb
        best = max(best, (fa + 1 - fb) / 2, (fb + 1 - fa) / 2)
    return best


rng = random.Random(53); out = []
for k in KS:
    sets = [{s[i:i + k] for i in range(len(s) - k + 1)} for _, _, _, s in rows]
    for name, A, B in pairs:
        ca, cb = counts(sets, A), counts(sets, B); na, nb = len(A), len(B)
        acc = best_acc(ca, na, cb, nb)
        null = collections.defaultdict(list); nacc = []
        pool = A + B
        for _ in range(PERMS):
            rng.shuffle(pool); pa, pb = pool[:na], pool[na:]
            qa, qb = counts(sets, pa), counts(sets, pb)
            for cut in CUTS: null[cut].append(sum(diag(qa, na, qb, nb, *cut)))
            nacc.append(best_acc(qa, na, qb, nb))
        for cut in CUTS:
            da, db = diag(ca, na, cb, nb, *cut)
            out.append([name, na, nb, k, cut[0], cut[1], da, db, da + db, max(null[cut]), round(acc, 3), round(max(nacc), 3)])
        print(name, k, 'cut 0.9/0.1:', out[-len(CUTS)][6:], file=sys.stderr)       # 2026-11-01: the 0.9 / 0.1 row (was the 0.8 / 0.2 row)
with open(OUT + '.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['comparison', 'n_a', 'n_b', 'k', 'p_in', 'p_out', 'diag_a', 'diag_b', 'diag_total', 'null_max_total', 'best_kmer_balanced_acc', 'null_max_acc'])
    w.writerows(out)
print('wrote', OUT + '.tsv')
# 2026-11-01 (CB13): k-report, the decision at --k-report and the k range where the METHODS rule holds (0.9 / 0.1 cut)
with open(OUT + '.k_report.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['comparison', 'k_report', 'diag_total_at_k', 'null_max_at_k', 'rule_holds_at_k', 'k_where_rule_holds', 'perms'])
    for name, _, _ in pairs:
        rows09 = {r[3]: r for r in out if r[0] == name and r[4] == 0.9}
        holds = [k for k, r in sorted(rows09.items()) if r[8] >= ARGS.min_diag and r[8] > r[9]]
        r = rows09.get(ARGS.k_report)
        at = (r[8], r[9], 'yes' if ARGS.k_report in holds else 'no') if r else ('NA', 'NA', 'NA')
        w.writerow([name, ARGS.k_report, *at, ','.join(map(str, holds)) or 'none', PERMS])
        print(f'k-report {name}: k={ARGS.k_report} diagnostic {at[0]} (null max {at[1]}) rule >= {ARGS.min_diag} & > null: {at[2]}; '
              f'rule holds for k in [{",".join(map(str, holds)) or "none"}]')
print('wrote', OUT + '.k_report.tsv')
