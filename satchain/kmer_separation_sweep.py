#!/usr/bin/env python3
"""Sensitivity sweep: at which k-mer length and diagnostic cutoff do two groups of copies separate?
For each pair of groups (step-03 classes, or the two domains either side of a step-53 boundary), counts k-mers present in
>= p_in of one group's copies and <= p_out of the other's, for k = 5..20 and four (p_in, p_out) cutoffs, and the same count
after shuffling the group labels (--perms permutations; group sizes kept). A k-mer is counted once per copy (presence).
Also reports the best single-k-mer classifier accuracy (balanced) per k.
usage: kmer_separation_sweep.py <03_monomer_copies.tsv> <array_id> <boundary_bp|-> <out prefix>"""
import csv, collections, random, sys
csv.field_size_limit(1 << 30)
COPIES, AID, BND, OUT = sys.argv[1:5]
KS = [5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20]
CUTS = [(0.9, 0.1), (0.8, 0.2), (0.7, 0.3), (0.6, 0.4)]
PERMS = 5

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
        print(name, k, out[-3][6:], file=sys.stderr)
with open(OUT + '.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['comparison', 'n_a', 'n_b', 'k', 'p_in', 'p_out', 'diag_a', 'diag_b', 'diag_total', 'null_max_total', 'best_kmer_balanced_acc', 'null_max_acc'])
    w.writerows(out)
print('wrote', OUT + '.tsv')
