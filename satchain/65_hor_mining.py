#!/usr/bin/env python3
"""
65_hor_mining.py  (SatChain v4.7 draft; HOR mining on the copy chain. Monomer alphabet + hierarchical tandem-repeat mining adapted from
HiCAT (Gao et al. 2023): threshold-swept block graph -> Louvain communities -> HTRM -> rotation merge -> coverage score.)

Differences from HiCAT: copies are SatChain confident copies (step 03 + HMM runs), not raw StringDecomposer blocks of a whole region;
every threshold is reported with its HOR coverage AND a shuffled-string null (stability table); the threshold is chosen by stability, not
by maximum coverage (2026-11-01, D04); --runs may list several arrays (cross-array mode): one shared alphabet is built from the pooled
copies, so HOR patterns are comparable between arrays.

1 Unique copy sequences (in array orientation; strand = majority strand of each run) with counts. Pairwise identity = global alignment
  score (match 1, mismatch 0, gap -1) / longer length, on unique sequences (<= --max-unique; above that, the most frequent ones form the
  core and the rest are assigned to their best match among ALL core sequences; D18, before: the 200 most frequent only).
2 For t in --t-min..--t-max by --t-step: graph (edge if identity >= t, weight = identity x sqrt(count_i x count_j); D14), Louvain
  communities = monomers.
3 Per run: monomer string in genomic order (array orientation). HTRM: repeatedly find the tandem repeat (unit length 2..--max-unit,
  >= --min-copies consecutive copies; default 3 since 2026-11-01: with 2, any two adjacent distinct monomers repeated twice are a 'HOR',
  coverage 0.66 on random 2-letter strings) covering the most monomers, replace each maximal run of the unit by one new symbol (the unit
  pattern is recorded), until none is found. Units that are rotations of each other are merged (canonical = smallest rotation).
4 HOR = pattern of >= 2 monomers; coverage = monomers inside HOR runs / all monomers; nesting = depth of the symbol hierarchy.
  Null (D04): per t and per run the monomer string is shuffled --null times (same alphabet and composition, no order) and mined the same
  way; null_coverage = mean, null_p95 = 95th percentile of the pooled coverage; excess = coverage - null_coverage.
  Threshold choice (stability): consecutive t with the same (top HOR length, monomer count) form a plateau; the chosen t is the middle of
  the longest plateau whose coverage exceeds the null p95; when no t exceeds its null, no HOR is reported (hors.tsv empty).
  Rank by coverage (ties: nesting); name R<rank>L<pattern length in monomers>.
Sweep (--sweep, D09): min_copies 2 / 3 / 4 x t -> coverage, null, excess, chosen t, written to <repo>/benchmarks/65_sweep.tsv (or --sweep-out).
Outputs (<dest>/<name>_): thresholds.tsv (t, monomers, HOR coverage, null coverage, null p95, excess, top HOR, plateau, chosen),
hors.tsv (chosen t: rank, name, pattern, length, copies, coverage, depth, runs), monomer_track.tsv (copy -> monomer at the chosen t).
usage: 65_hor_mining.py --catalogue <59> --copies <03> [...] --runs id1,id2 --dest <dir> --name <n> [--t-min 0.90 --t-max 0.995 --t-step 0.005]
       [--max-unit 40] [--min-copies 3] [--min-identity 60] [--max-unique 1500] [--null 20] [--sweep] [--sweep-out path]
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


def coverage_of(strings, max_unit, min_copies):
    """mine every monomer string; returns coverage, patterns, per-run (string, compressed), HOR occurrence counts"""
    patterns = {}; total = covered = 0; per_run = {}; hor_counts = collections.Counter()
    for rid, ms in strings.items():
        comp = htrm(list(ms), max_unit, min_copies, patterns); per_run[rid] = (ms, comp); total += len(ms)
        covered += sum(len(expand(x, patterns)) for x in comp if x in patterns)
        for x in comp:
            if x in patterns: hor_counts[x] += 1
    return covered / max(total, 1), patterns, per_run, hor_counts


def null_coverage(strings, max_unit, min_copies, reps, rng):
    """coverage of the same strings with the monomer order shuffled per run (D04), reps times"""
    out = []
    for _ in range(reps):
        sh = {}
        for rid, ms in strings.items():
            s = list(ms); rng.shuffle(s); sh[rid] = s
        out.append(coverage_of(sh, max_unit, min_copies)[0])
    return out


def choose_t(rows):
    """stability rule: longest plateau of consecutive t with the same (top_hor_len, monomers) among rows above their null p95; middle t.
    Returns the chosen row index or None; marks rows['plateau'] with the plateau id."""
    pid = 0; prev = None
    for r in rows:
        key = (r['top_hor_len'], r['monomers'])
        if key != prev: pid += 1; prev = key
        r['plateau'] = pid
    good = collections.defaultdict(list)
    for i, r in enumerate(rows):
        if r['hor_coverage'] > r['null_p95'] and r['top_hor_len'] > 0: good[r['plateau']].append(i)
    if not good: return None
    best = max(good.values(), key=lambda v: (len(v), rows[v[0]]['hor_coverage'] - rows[v[0]]['null_coverage']))
    return best[len(best) // 2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--catalogue', required=True); ap.add_argument('--copies', nargs='+', required=True); ap.add_argument('--runs', required=True)
    ap.add_argument('--dest', required=True); ap.add_argument('--name', required=True)
    ap.add_argument('--t-min', type=float, default=0.90); ap.add_argument('--t-max', type=float, default=0.995); ap.add_argument('--t-step', type=float, default=0.005)
    ap.add_argument('--max-unit', type=int, default=40); ap.add_argument('--min-copies', type=int, default=3); ap.add_argument('--null', type=int, default=20)
    ap.add_argument('--min-identity', type=float, default=60); ap.add_argument('--max-unique', type=int, default=1500); ap.add_argument('--seed', type=int, default=65)
    ap.add_argument('--sweep', action='store_true'); ap.add_argument('--sweep-out', default='')
    a = ap.parse_args(); os.makedirs(a.dest, exist_ok=True); want = a.runs.split(','); rng = random.Random(a.seed)
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
    for s in uniq[a.max_unique:]:                      # D18: best match among ALL core sequences
        sc = [GA.score(s, c) / max(len(s), len(c)) for c in core]; assign[s] = int(np.argmax(sc))
    w = np.array([cnt[s] for s in core], float)
    ts = []; t = a.t_min
    while t <= a.t_max + 1e-9: ts.append(round(t, 4)); t += a.t_step
    comms = {}
    for t in ts:
        Wm = np.where(S >= t, S, 0.0) * np.sqrt(np.outer(w, w)); np.fill_diagonal(Wm, 0)
        comms[t] = louvain(Wm, a.seed)

    def sweep_mc(mc):
        rows, keep = [], {}
        for t in ts:
            comm = comms[t]; strings = {rid: [int(comm[assign[s]]) for s in seqs] for rid, (L, seqs) in chains.items()}
            cov, patterns, per_run, hor_counts = coverage_of(strings, a.max_unit, mc)
            nul = null_coverage(strings, a.max_unit, mc, a.null, rng)
            top = hor_counts.most_common(1)
            rows.append(dict(t=t, min_copies=mc, monomers=len(set(comm)), hor_coverage=round(cov, 4), null_coverage=round(float(np.mean(nul)), 4),
                             null_p95=round(float(np.percentile(nul, 95)), 4), excess=round(cov - float(np.mean(nul)), 4), n_hor_patterns=len(hor_counts),
                             top_hor_len=len(expand(top[0][0], patterns)) if top else 0, top_hor_copies=top[0][1] if top else 0))
            keep[t] = (comm, patterns, per_run, hor_counts)
            print(rows[-1], flush=True)
        ch = choose_t(rows)
        for i, r in enumerate(rows): r['chosen'] = int(i == ch)
        return rows, keep, ch

    rows, keep, ch = sweep_mc(a.min_copies)
    with open(os.path.join(a.dest, f'{a.name}_thresholds.tsv'), 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t', lineterminator='\n'); wr.writeheader(); wr.writerows(rows)
    hs = []
    if ch is None: tb = ts[0]; comm = keep[tb][0]; print('no threshold exceeds its shuffled null: no HOR reported')
    else:
        tb = rows[ch]['t']; comm, patterns, per_run, hc = keep[tb]
        total = sum(len(ms) for ms, _ in per_run.values())
        for sym, k in hc.items():
            L = len(expand(sym, patterns)); cv = k * L / max(total, 1)
            hs.append(dict(pattern='-'.join(str(x) for x in expand(sym, patterns)), length=L, copies=k, coverage=round(cv, 4), depth=patterns[sym][1],
                           runs=sum(1 for _, comp in per_run.values() if sym in comp)))
        hs.sort(key=lambda h: (-h['coverage'], -h['depth']))
        for i, h in enumerate(hs, 1): h['name'] = f"R{i}L{h['length']}"
    with open(os.path.join(a.dest, f'{a.name}_hors.tsv'), 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=['name', 'pattern', 'length', 'copies', 'coverage', 'depth', 'runs'], delimiter='\t', lineterminator='\n'); wr.writeheader(); wr.writerows(hs)
    with open(os.path.join(a.dest, f'{a.name}_monomer_track.tsv'), 'w') as fh:
        fh.write('run_id\tstart\tend\tmonomer\n')
        for rid, (L, seqs) in chains.items():
            for (s, e, st, _), sq in zip(L, seqs): fh.write(f'{rid}\t{s}\t{e}\t{int(comm[assign[sq]])}\n')
    if ch is not None: print(f'chosen t {tb:.3f} (plateau of {sum(1 for r in rows if r["plateau"] == rows[ch]["plateau"])} t values): HOR coverage {rows[ch]["hor_coverage"]:.3f} '
                             f'vs null {rows[ch]["null_coverage"]:.3f}; top HORs:', [(h["name"], h["pattern"][:60], h["copies"]) for h in hs[:5]])
    if a.sweep:
        SW = [dict(step='65', name=a.name, **r) for r in rows]
        for mc in (2, 3, 4):
            if mc == a.min_copies: continue
            SW += [dict(step='65', name=a.name, **r) for r in sweep_mc(mc)[0]]
        out = a.sweep_out or os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'benchmarks', '65_sweep.tsv')
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out, 'w', newline='') as fh:
            wr = csv.DictWriter(fh, fieldnames=list(SW[0].keys()), delimiter='\t', lineterminator='\n'); wr.writeheader(); wr.writerows(SW)
        print('sweep ->', out)


if __name__ == '__main__':
    main()
