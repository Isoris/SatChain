#!/usr/bin/env python3
"""
65_hor_mining.py  (SatChain v4.7 draft; HOR mining on the copy chain. Monomer alphabet + hierarchical tandem-repeat mining adapted from
HiCAT (Gao et al. 2023): threshold-swept block graph -> Louvain communities -> HTRM -> rotation merge -> coverage score.)

Differences from HiCAT: copies are SatChain confident copies (step 03 + HMM runs), not raw StringDecomposer blocks of a whole region;
every threshold is reported with its HOR coverage (stability table), the max-coverage threshold is only the default pick; --runs may
list several arrays (cross-array mode): one shared alphabet is built from the pooled copies, so HOR patterns are comparable between arrays.

1 Unique copy sequences (in array orientation; strand = majority strand of each run) with counts. Pairwise identity = global alignment
  score (match 1, mismatch 0, gap -1) / longer length, on unique sequences (<= --max-unique; above that, the most frequent ones and the
  rest assigned to their best match).
2 For t in --t-min..--t-max by --t-step: graph (edge if identity >= t, weight = identity x counts), Louvain communities = monomers.
3 Per run: monomer string in genomic order (array orientation). HTRM: repeatedly find the tandem repeat (unit length 2..--max-unit,
  >= --min-copies consecutive copies) covering the most monomers, replace each maximal run of the unit by one new symbol (the unit
  pattern is recorded), until none is found. Units that are rotations of each other are merged (canonical = smallest rotation).
4 HOR = pattern of >= 2 monomers; coverage = monomers inside HOR runs / all monomers; nesting = depth of the symbol hierarchy.
  Rank by coverage (ties: nesting); name R<rank>L<pattern length in monomers>.
Outputs (<dest>/<name>_): thresholds.tsv (t, monomers, HOR coverage, top HOR), hors.tsv (best t: rank, name, pattern, length, copies,
coverage, depth, runs), monomer_track.tsv (copy -> monomer, HOR id at best t).
usage: 65_hor_mining.py --catalogue <59> --copies <03> [...] --runs id1,id2 --dest <dir> --name <n> [--t-min 0.90 --t-max 0.995 --t-step 0.005]
       [--max-unit 40] [--min-copies 2] [--min-identity 60] [--max-unique 1500]
"""
import argparse, collections, csv, os, random
import numpy as np
from Bio import Align
csv.field_size_limit(1 << 30)
RC = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(RC)[::-1]
GA = Align.PairwiseAligner(mode='global', match_score=1, mismatch_score=0, open_gap_score=-1, extend_gap_score=-1)


def louvain(W, seed=0, passes=2):
    """small Louvain (local moving + one aggregation per pass) on a dense symmetric weight matrix; returns community per node"""
    rng = np.random.default_rng(seed); n = len(W); comm = np.arange(n); nodes_of = [[i] for i in range(n)]; A = W.copy()
    for _ in range(passes):
        m = A.sum() / 2
        if m == 0: break
        k = A.sum(1); c = np.arange(len(A)); tot = k.copy(); moved = True
        while moved:
            moved = False
            for i in rng.permutation(len(A)):
                ci = c[i]; tot[ci] -= k[i]
                nb = np.nonzero(A[i])[0]; w = collections.defaultdict(float)
                for j in nb:
                    if j != i: w[c[j]] += A[i, j]
                best, gain = ci, w.get(ci, 0.0) - tot[ci] * k[i] / (2 * m)
                for cj, wij in w.items():
                    g = wij - tot[cj] * k[i] / (2 * m)
                    if g > gain + 1e-12: best, gain = cj, g
                c[i] = best; tot[best] += k[i]
                if best != ci: moved = True
        labs = {x: i for i, x in enumerate(sorted(set(c)))}; c = np.array([labs[x] for x in c])
        new_nodes = [[] for _ in labs]
        for i, ci in enumerate(c): new_nodes[ci] += nodes_of[i]
        if len(new_nodes) == len(A): break
        B = np.zeros((len(labs), len(labs)))
        for i in range(len(A)):
            for j in np.nonzero(A[i])[0]: B[c[i], c[j]] += A[i, j]
        A = B; nodes_of = new_nodes
    for ci, ns in enumerate(nodes_of):
        for x in ns: comm[x] = ci
    return comm


def canonical(unit):
    rots = [tuple(unit[i:] + unit[:i]) for i in range(len(unit))]; return min(rots, key=repr)


def htrm(seq, max_unit, min_copies, patterns):
    """seq: list of symbols (ints or tuples). Returns compressed sequence; patterns: dict symbol -> (unit, depth). Repeats until no TR."""
    while True:
        best = None; n = len(seq)
        for L in range(2, min(max_unit, n // min_copies) + 1):
            i = 0
            while i + L <= n:
                j = i + L
                while j + L <= n and seq[j:j + L] == seq[i:i + L]: j += L
                cop = (j - i) // L
                if cop >= min_copies and len(set(seq[i:i + L])) > 1:
                    cov = j - i
                    if best is None or cov > best[0] or (cov == best[0] and L < best[1]): best = (cov, L, tuple(seq[i:i + L]))
                    i = j
                else: i += 1
        if best is None: return seq
        _, L, unit = best; key = canonical(list(unit))
        sym = ('H', key); depth = 1 + max((patterns[x][1] for x in unit if x in patterns), default=0)
        patterns[sym] = (key, depth)
        out, i = [], 0
        rots = {tuple(key[r:] + key[:r]) for r in range(len(key))}
        while i < len(seq):
            if i + L <= len(seq) and tuple(seq[i:i + L]) in rots:
                j = i + L
                while j + L <= len(seq) and tuple(seq[j:j + L]) == tuple(seq[i:i + L]): j += L
                if (j - i) // L >= min_copies or (out and out[-1] == sym):
                    out.extend([sym] * ((j - i) // L)); i = j; continue
            out.append(seq[i]); i += 1
        # collapse consecutive identical HOR symbols into one occurrence-counted run (keep one symbol per copy of the unit)
        if out == seq: return seq
        seq = out


def expand(sym, patterns):
    if sym in patterns: return sum((expand(x, patterns) for x in patterns[sym][0]), [])
    return [sym]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--catalogue', required=True); ap.add_argument('--copies', nargs='+', required=True); ap.add_argument('--runs', required=True)
    ap.add_argument('--dest', required=True); ap.add_argument('--name', required=True)
    ap.add_argument('--t-min', type=float, default=0.90); ap.add_argument('--t-max', type=float, default=0.995); ap.add_argument('--t-step', type=float, default=0.005)
    ap.add_argument('--max-unit', type=int, default=40); ap.add_argument('--min-copies', type=int, default=2)
    ap.add_argument('--min-identity', type=float, default=60); ap.add_argument('--max-unique', type=int, default=1500); ap.add_argument('--seed', type=int, default=65)
    a = ap.parse_args(); os.makedirs(a.dest, exist_ok=True); want = a.runs.split(',')
    runs = {r['run_id']: r for r in csv.DictReader(open(a.catalogue), delimiter='\t') if r['run_id'] in want}
    span = collections.defaultdict(list)
    for r in runs.values(): span[r['chrom']].append((int(r['start0']), int(r['end']), r['run_id']))
    cp = collections.defaultdict(list)
    for p in a.copies:
        with open(p, newline='') as fh:
            rd = csv.reader(fh, delimiter='\t'); h = next(rd); ix = {k: i for i, k in enumerate(h)}
            for row in rd:
                c = row[ix['chrom']]
                if c not in span or float(row[ix['best_alignment_identity']] or 0) < a.min_identity: continue
                s, e = int(row[ix['genomic_start_0based']]), int(row[ix['genomic_end']])
                for rs, re_, rid in span[c]:
                    if s >= rs and e <= re_: cp[rid].append((s, e, row[ix['strand']], row[ix['sequence']].upper())); break
    chains = {}
    for rid, L in cp.items():
        L.sort(); rev = sum(1 for x in L if x[2] == '-') > len(L) / 2
        seqs = [x[3] for x in L]                      # copy sequences are already in consensus orientation (step 03)
        chains[rid] = (L[::-1] if rev else L, seqs[::-1] if rev else seqs)
    cnt = collections.Counter(s for _, seqs in chains.values() for s in seqs)
    uniq = [s for s, _ in cnt.most_common()]; core = uniq[:a.max_unique]; n = len(core)
    print(f'copies {sum(cnt.values())}, unique {len(uniq)}, aligned {n}', flush=True)
    S = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n): S[i, j] = S[j, i] = GA.score(core[i], core[j]) / max(len(core[i]), len(core[j]))
    assign = {s: i for i, s in enumerate(core)}
    for s in uniq[a.max_unique:]:
        sc = [GA.score(s, c) / max(len(s), len(c)) for c in core[:200]]; assign[s] = int(np.argmax(sc))
    w = np.array([cnt[s] for s in core], float)
    rows, best = [], None
    t = a.t_min
    while t <= a.t_max + 1e-9:
        Wm = np.where(S >= t, S, 0.0) * np.sqrt(np.outer(w, w)); np.fill_diagonal(Wm, 0)
        comm = louvain(Wm, a.seed)
        patterns = {}; total = covered = 0; per_run = {}
        for rid, (L, seqs) in chains.items():
            ms = [int(comm[assign[s]]) for s in seqs]; total += len(ms)
            comp = htrm(ms, a.max_unit, a.min_copies, patterns)
            per_run[rid] = (ms, comp)
            covered += sum(len(expand(x, patterns)) for x in comp if x in patterns)
        cov = covered / max(total, 1)
        hor_counts = collections.Counter()
        for rid, (ms, comp) in per_run.items():
            for x in comp:
                if x in patterns: hor_counts[x] += 1
        top = hor_counts.most_common(1)
        rows.append(dict(t=round(t, 4), monomers=len(set(comm)), hor_coverage=round(cov, 4), n_hor_patterns=len(hor_counts),
                         top_hor_len=len(expand(top[0][0], patterns)) if top else 0, top_hor_copies=top[0][1] if top else 0))
        print(rows[-1], flush=True)
        if best is None or cov > best[0] + 1e-9: best = (cov, t, comm.copy(), patterns, per_run, hor_counts)
        t += a.t_step
    with open(os.path.join(a.dest, f'{a.name}_thresholds.tsv'), 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t'); wr.writeheader(); wr.writerows(rows)
    cov, tb, comm, patterns, per_run, hc = best
    total = sum(len(ms) for ms, _ in per_run.values())
    hs = []
    for sym, k in hc.items():
        L = len(expand(sym, patterns)); cv = k * L / max(total, 1)
        hs.append(dict(pattern='-'.join(str(x) for x in expand(sym, patterns)), length=L, copies=k, coverage=round(cv, 4), depth=patterns[sym][1],
                       runs=sum(1 for _, comp in per_run.values() if sym in comp)))
    hs.sort(key=lambda h: (-h['coverage'], -h['depth']))
    for i, h in enumerate(hs, 1): h['name'] = f"R{i}L{h['length']}"
    with open(os.path.join(a.dest, f'{a.name}_hors.tsv'), 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=['name', 'pattern', 'length', 'copies', 'coverage', 'depth', 'runs'], delimiter='\t'); wr.writeheader(); wr.writerows(hs)
    with open(os.path.join(a.dest, f'{a.name}_monomer_track.tsv'), 'w') as fh:
        fh.write('run_id\tstart\tend\tmonomer\n')
        for rid, (L, seqs) in chains.items():
            for (s, e, st, _), sq in zip(L, seqs): fh.write(f'{rid}\t{s}\t{e}\t{int(comm[assign[sq]])}\n')
    print(f'best t {tb:.3f}: HOR coverage {cov:.3f}; top HORs:', [(h["name"], h["pattern"][:60], h["copies"]) for h in hs[:5]])


if __name__ == '__main__':
    main()
