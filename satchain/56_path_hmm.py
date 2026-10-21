#!/usr/bin/env python3
"""
56_path_hmm.py  (v4.6 prototype; arrays and runs from the copy chain, replacing the fixed run rules of step 36)

Every step-03 tile of a candidate interval is one observation along the chromosome, in order:
  C = confident copy (complete, identity >= --conf), W = weak copy (identity --weak..conf, or truncated with identity >= --conf),
  J = junk tile (identity < --weak), G = uncovered stretch longer than --gap-units x period (one G per stretch)
Two hidden states: A (array) and B (background). Transitions and emissions are learned by Baum-Welch over all candidate paths of
one subgenome (start from A: C .70 W .20 J .08 G .02 ; B: C .05 W .15 J .50 G .30 ; stay .98). The posterior P(A) of each tile is
computed by forward-backward; an array run = a maximal stretch of tiles with P(A) >= --post (swept), trimmed to its first and last
C tile, kept if it holds >= --min-c confident copies. Short interruptions (TE, degraded copies) stay inside a run when the chain
continues on both sides.
Outputs (<dest>/): 56_hmm_params.tsv, 56_runs.bed (BED6, name = array|run|n_conf|n_tiles|mean_post), 56_vs_36.tsv (bp overlap with the
step-36 runs: shared, HMM-only, 36-only), 56_summary.tsv.
usage: 56_path_hmm.py --out <pipeline out> --dest <dir> [--post 0.5,0.8,0.95]
"""
import argparse, collections, csv, math, os
import numpy as np
csv.field_size_limit(1 << 30)
SYM = {'C': 0, 'W': 1, 'J': 2, 'G': 3}


def paths(out, conf, weak, gapu):
    arr = {r['array_id']: r for r in csv.DictReader(open(os.path.join(out, '02_satellite_arrays.tsv')), delimiter='\t')}
    cop = collections.defaultdict(list)
    with open(os.path.join(out, '03_monomer_copies.tsv'), newline='') as fh:
        rd = csv.reader(fh, delimiter='\t'); h = next(rd); ix = {k: i for i, k in enumerate(h)}
        for row in rd:
            try: idt = float(row[ix['best_alignment_identity']])
            except ValueError: idt = 0.0
            cop[row[ix['array_id']]].append((int(row[ix['genomic_start_0based']]), int(row[ix['genomic_end']]), idt, row[ix['complete_or_truncated']] == 'complete'))
    P = {}
    for aid, cs in cop.items():
        cs.sort(); per = float(arr[aid]['period_mode_bp'] or 0) if aid in arr else 0
        obs, pos = [], []; last = None
        for s, e, idt, comp in cs:
            if last is not None and s - last > gapu * max(1, per): obs.append(3); pos.append((last, s))
            o = 0 if (comp and idt >= conf) else 1 if idt >= weak and (comp or idt >= conf) else 1 if idt >= weak else 2
            obs.append(o); pos.append((s, e)); last = e if last is None else max(last, e)
        P[aid] = (arr[aid]['chrom'] if aid in arr else cs[0][0], np.array(obs), pos)
    return P, arr


def fb(obs, T, E, pi):
    n = len(obs); a = np.zeros((n, 2)); c = np.zeros(n)
    a[0] = pi * E[:, obs[0]]; c[0] = a[0].sum(); a[0] /= c[0]
    for t in range(1, n):
        a[t] = (a[t - 1] @ T) * E[:, obs[t]]; c[t] = a[t].sum(); a[t] /= c[t]
    b = np.ones((n, 2))
    for t in range(n - 2, -1, -1): b[t] = (T @ (E[:, obs[t + 1]] * b[t + 1])) / c[t + 1]
    g = a * b; g /= g.sum(1, keepdims=True)
    return a, b, c, g


def baum_welch(seqs, iters=30):
    T = np.array([[0.98, 0.02], [0.02, 0.98]]); E = np.array([[0.70, 0.20, 0.08, 0.02], [0.05, 0.15, 0.50, 0.30]]); pi = np.array([0.5, 0.5]); ll0 = None
    for it in range(iters):
        Tn = np.zeros((2, 2)); En = np.zeros((2, 4)); pin = np.zeros(2); ll = 0
        for o in seqs:
            if len(o) < 2: continue
            a, b, c, g = fb(o, T, E, pi); ll += np.log(c).sum(); pin += g[0]
            for k in range(4): En[:, k] += g[o == k].sum(0)
            xi = a[:-1, :, None] * T[None] * (E[:, o[1:]].T * b[1:])[:, None, :] / c[1:, None, None]
            Tn += xi.sum(0)
        T = Tn / Tn.sum(1, keepdims=True); E = (En + 1e-6) / (En + 1e-6).sum(1, keepdims=True); pi = pin / pin.sum()
        if E[0, 0] < E[1, 0]: T, E, pi = T[::-1, ::-1], E[::-1], pi[::-1]       # keep state 0 = array (more C)
        if ll0 is not None and abs(ll - ll0) < 1e-4 * abs(ll0): break
        ll0 = ll
    return T, E, pi, ll, it + 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True); ap.add_argument('--dest', required=True)
    ap.add_argument('--conf', type=float, default=80); ap.add_argument('--weak', type=float, default=65); ap.add_argument('--gap-units', type=float, default=0.5)
    ap.add_argument('--post', default='0.5,0.8,0.95'); ap.add_argument('--min-c', type=int, default=4)
    ap.add_argument('--train-n', type=int, default=400, help='arrays used to learn the parameters (seeded sample)')
    a = ap.parse_args(); os.makedirs(a.dest, exist_ok=True)
    P, arr = paths(a.out, a.conf, a.weak, a.gap_units)
    keys = sorted(P); rs = np.random.default_rng(56); tr = rs.choice(len(keys), min(a.train_n, len(keys)), replace=False)
    seqs = [P[keys[i]][1] for i in tr]
    T, E, pi, ll, nit = baum_welch(seqs)
    with open(os.path.join(a.dest, '56_hmm_params.tsv'), 'w') as fh:
        fh.write(f'iterations\t{nit}\nloglik\t{ll:.1f}\nP(A->A)\t{T[0,0]:.5f}\nP(B->B)\t{T[1,1]:.5f}\n')
        for i, st in enumerate('AB'):
            fh.write(f'emission_{st}\t' + '\t'.join(f'{k}:{E[i, j]:.4f}' for k, j in SYM.items()) + '\n')
    # step-36 runs for comparison
    r36 = collections.defaultdict(list)
    p36 = os.path.join(a.out, '36_characterized_intervals.bed')
    for l in (open(p36) if os.path.exists(p36) else []):
        if l.startswith(('#', 'track')): continue
        f = l.split('\t'); r36[f[0]].append((int(f[1]), int(f[2])))
    summ = []
    for post in [float(x) for x in a.post.split(',')]:
        runs = collections.defaultdict(list); bed = open(os.path.join(a.dest, f'56_runs_post{post}.bed'), 'w'); n = 0
        for aid, (chrom, o, pos) in P.items():
            if len(o) < 2: continue
            _, _, _, g = fb(o, T, E, pi); inA = g[:, 0] >= post; i = 0
            while i < len(o):
                if not inA[i]: i += 1; continue
                j = i
                while j + 1 < len(o) and inA[j + 1]: j += 1
                cidx = [k for k in range(i, j + 1) if o[k] == 0]
                if len(cidx) >= a.min_c:
                    s, e = pos[cidx[0]][0], pos[cidx[-1]][1]; n += 1
                    runs[chrom].append((s, e))
                    bed.write(f'{chrom}\t{s}\t{e}\t{aid}|hmm{n}|{len(cidx)}|{j - i + 1}|{g[i:j + 1, 0].mean():.3f}\t0\t.\n')
                i = j + 1
        bed.close()
        def bp(iv):
            tot = 0
            for c, v in iv.items():
                v = sorted(v); cs_, ce = None, None
                for s, e in v:
                    if cs_ is None or s > ce:
                        if cs_ is not None: tot += ce - cs_
                        cs_, ce = s, e
                    else: ce = max(ce, e)
                if cs_ is not None: tot += ce - cs_
            return tot
        def inter(x, y):
            tot = 0
            for c in x:
                X = sorted(x[c]); Y = sorted(y.get(c, [])); j = 0
                for s, e in X:
                    while j < len(Y) and Y[j][1] <= s: j += 1
                    k = j
                    while k < len(Y) and Y[k][0] < e: tot += max(0, min(e, Y[k][1]) - max(s, Y[k][0])); k += 1
            return tot
        bh, b36 = bp(runs), bp(r36); sh = inter(runs, r36)
        summ.append([post, n, round(bh / 1e6, 2), round(b36 / 1e6, 2), round(sh / 1e6, 2), round((bh - sh) / 1e6, 2), round((b36 - sh) / 1e6, 2)])
    with open(os.path.join(a.dest, '56_summary.tsv'), 'w') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['posterior_cut', 'n_hmm_runs', 'hmm_Mb', 'step36_Mb', 'shared_Mb', 'hmm_only_Mb', 'step36_only_Mb']); w.writerows(summ)
    print(open(os.path.join(a.dest, '56_hmm_params.tsv')).read()); [print(s) for s in summ]


if __name__ == '__main__':
    main()
