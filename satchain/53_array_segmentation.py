#!/usr/bin/env python3
"""
53_array_segmentation.py  (v4.5.2; windowed features, local topology and segmentation INSIDE arrays; read-only)

Grain: sliding windows inside each step-02 array, W = max(--min-win-bp, --win-units x period), step W/2
(short-unit arrays get >= 1 kb, long-unit arrays >= 20 units). Each window carries five feature blocks:
  comp    GC, CpG obs/exp, normalised entropy of 1/2/3/5-mers                        (sequence, label-free)
  kmer    canonical 4-mer profile (136) -> first --kmer-pcs principal components per array (sequence, label-free)
  period  k-mer (k=10) next-recurrence: dominant distance (local period, bp), its strength,
          strength at the array period (02), fraction of k-mers that recur at all     (sequence, label-free)
  label   step-03 copies whose midpoint lies in the window: median length / period, length IQR / period,
          mean identity, confident fraction (identity >= 80, complete), dominant-class fraction, n classes,
          fraction of the array's top-3 families                                     (decomposition labels)
  topo    local transition graph of confident copies (class+strand): self-loop fraction, order-1 entropy,
          (AB)n alternation, distinct edges per transition                            (decomposition labels)
Segmentation: per block, binary segmentation on the window series (features z-scored within the array; cost =
within-segment sum of squares; split while gain > --penalty x d x log n; min segment --min-seg windows).
A boundary is ACCEPTED when >= 2 blocks place a change point within +/- 1 window AND at least one of them is
label-free (comp, kmer, period), so wrong monomer calls cannot confirm themselves. Single-block change points are
kept as 'candidate'.
Segments (between accepted boundaries) get a profile from their confident copies:
  unresolved (< 10 confident copies or confident fraction < 0.3) | interleaved_families (>= 2 families each >= 20 %)
  | alternating (alternation >= 0.5) | variant_mixed (class self-loop < 0.6) | homogeneous_monomeric
  plus local period, dominant family and its fraction. HOR / composite calls stay with steps 30/34/47.
Outputs (<out>/53_segmentation/): 53_window_features.tsv.gz, 53_change_points.tsv, 53_segments.tsv,
53_segments.bed (BED9 coloured by profile), 53_array_summary.tsv, 53_summary.tsv.
NOT done here: re-decomposition of each side of a boundary (step 50 local library) - the validation stage.
"""
import argparse, collections, csv, gzip, math, os, sys
import numpy as np

csv.field_size_limit(1 << 30)
ENC = np.full(256, -1, np.int8)
for i, c in enumerate(b'ACGT'): ENC[c] = i; ENC[c + 32] = i
COLOR = {'too_short': '210,210,210', 'homogeneous_monomeric': '31,119,180', 'variant_mixed': '44,160,44', 'alternating': '148,103,189',
         'interleaved_families': '214,39,40', 'unresolved': '150,150,150'}


def read_fasta(path, need):
    seqs, name, buf = {}, None, []
    with open(path, 'rb') as fh:
        for l in fh:
            if l.startswith(b'>'):
                if name in need: seqs[name] = b''.join(buf)
                name = l[1:].split()[0].decode(); buf = []
            elif name in need: buf.append(l.rstrip())
    if name in need: seqs[name] = b''.join(buf)
    return seqs


def kcodes(a, k):
    """base-4 codes of all k-mers; -1 where the k-mer has a non-ACGT base"""
    n = len(a) - k + 1
    if n <= 0: return np.empty(0, np.int64)
    bad = np.convolve((a < 0).astype(np.int32), np.ones(k, np.int32), 'valid') > 0
    x = np.where(a < 0, 0, a).astype(np.int64); c = np.zeros(n, np.int64)
    for j in range(k): c = c * 4 + x[j:j + n]
    c[bad] = -1
    return c


def ent_norm(codes, k):
    v = codes[codes >= 0]
    if len(v) < 2: return np.nan
    _, cnt = np.unique(v, return_counts=True); p = cnt / cnt.sum()
    return float(-(p * np.log2(p)).sum() / math.log2(min(4 ** k, len(v))))


CANON = None


def canon_index():
    global CANON
    if CANON is None:
        idx, m = {}, np.zeros(256, np.int64)
        for c in range(256):
            s = [(c >> (2 * (3 - j))) & 3 for j in range(4)]; r = [3 - b for b in s[::-1]]
            rc = sum(b << (2 * (3 - j)) for j, b in enumerate(r)); key = min(c, rc)
            m[c] = idx.setdefault(key, len(idx))
        CANON = m
    return CANON


def window_seq_features(a, period, k=10):
    f = {}
    n = int((a >= 0).sum()); f['valid_bp'] = n
    if n < 50: return f, None
    g = int(((a == 1) | (a == 2)).sum()); f['GC'] = g / n
    c2 = kcodes(a, 2); v2 = c2[c2 >= 0]
    cg = int((v2 == 6).sum()); fc, fg = (a == 1).sum() / n, (a == 2).sum() / n
    f['CpG_oe'] = (cg / max(1, len(v2))) / (fc * fg) if fc * fg > 0 else np.nan
    for kk in (1, 2, 3, 5): f[f'ent{kk}'] = ent_norm(kcodes(a, kk) if kk > 1 else a[a >= 0].astype(np.int64), kk)
    c4 = kcodes(a, 4); c4 = c4[c4 >= 0]
    prof = np.bincount(canon_index()[c4], minlength=136).astype(float); prof /= max(1, prof.sum())
    ck = kcodes(a, k); pos = np.nonzero(ck >= 0)[0]
    if len(pos) > 20:
        o = np.lexsort((pos, ck[pos])); cs, ps = ck[pos][o], pos[o]
        same = cs[1:] == cs[:-1]; d = (ps[1:] - ps[:-1])[same]
        f['recur_frac'] = len(d) / len(pos)
        if len(d):
            h = np.bincount(np.minimum(d, 20000))
            h[:2] = 0; lp = int(h.argmax()) if h.sum() else 0
            tol = lambda L: (d >= L * 0.94) & (d <= L * 1.06)
            f['local_period'] = lp; f['period_strength'] = float(tol(lp).mean()) if lp else 0.0
            f['strength_at_02'] = float(tol(period).mean()) if period >= 2 else np.nan
    return f, prof


def binseg(X, pen, minseg):
    """indices (window positions) of change points; X n x d, already z-scored"""
    n, d = X.shape
    if n < 2 * minseg: return []
    beta = pen * d * math.log(n); cps = []
    cs = np.vstack([np.zeros(d), np.cumsum(X, 0)]); cs2 = np.concatenate([[0], np.cumsum((X ** 2).sum(1))])
    cost = lambda i, j: (cs2[j] - cs2[i]) - ((cs[j] - cs[i]) ** 2).sum() / (j - i)
    stack = [(0, n)]
    while stack:
        i, j = stack.pop()
        if j - i < 2 * minseg: continue
        base = cost(i, j); best, bt = 0, None
        for t in range(i + minseg, j - minseg + 1):
            gain = base - cost(i, t) - cost(t, j)
            if gain > best: best, bt = gain, t
        if bt is not None and best > beta:
            cps.append(bt); stack += [(i, bt), (bt, j)]
    return sorted(cps)


def zs(M):
    M = np.array(M, float)
    if M.ndim == 1: M = M[:, None]
    col = np.nanmedian(M, 0); M = np.where(np.isnan(M), col, M)
    sd = M.std(0); keep = sd > 1e-9
    return (M[:, keep] - M[:, keep].mean(0)) / sd[keep] if keep.any() else np.zeros((len(M), 0))


def topo(labels):
    tr = [(x, y) for x, y in zip(labels, labels[1:]) if x is not None and y is not None]
    if not tr: return np.nan, np.nan, np.nan, np.nan
    sl = sum(x == y for x, y in tr) / len(tr)
    out = collections.defaultdict(collections.Counter)
    for x, y in tr: out[x][y] += 1
    h = sum(sum(c.values()) / len(tr) * -sum(v / sum(c.values()) * math.log2(v / sum(c.values())) for v in c.values()) for c in out.values())
    tri = [(x, y, z) for x, y, z in zip(labels, labels[1:], labels[2:]) if None not in (x, y, z)]
    alt = sum(x == z and x != y for x, y, z in tri) / len(tri) if tri else np.nan
    return sl, h, alt, len(set(tr)) / len(tr)


def seg_profile(cs):
    conf = [c for c in cs if c['conf']]
    fam = collections.Counter(c['fam'] for c in conf)
    lab = [c['cls'] if c['conf'] else None for c in cs]
    sl, h, alt, _ = topo(lab)
    if len(cs) < 10: p = 'too_short'
    elif len(conf) < 10 or len(conf) / max(1, len(cs)) < 0.3: p = 'unresolved'
    elif sum(1 for f, n in fam.items() if n >= 0.2 * len(conf)) >= 2: p = 'interleaved_families'
    elif alt == alt and alt >= 0.5: p = 'alternating'
    elif sl == sl and sl < 0.6: p = 'variant_mixed'
    else: p = 'homogeneous_monomeric'
    top = fam.most_common(1)[0] if fam else ('NA', 0)
    return p, top[0], round(top[1] / len(conf), 3) if conf else 'NA', len(cs), len(conf), sl, alt


def r(x, n=4): return 'NA' if x is None or (isinstance(x, float) and x != x) else (round(x, n) if isinstance(x, float) else x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True); ap.add_argument('--genome', required=True)
    ap.add_argument('--dest'); ap.add_argument('--arrays', help='comma list or file of array_ids (default all)')
    ap.add_argument('--runs', help='catalogue runs (tsv: run_id, array_id, chrom, start0, end, unit_bp): segment runs instead of step-02 arrays')
    ap.add_argument('--min-win-bp', type=int, default=1000); ap.add_argument('--win-units', type=int, default=20)
    ap.add_argument('--max-win-bp', type=int, default=60000); ap.add_argument('--kmer-pcs', type=int, default=5)
    ap.add_argument('--penalty', type=float, default=3.0); ap.add_argument('--min-seg', type=int, default=3)
    ap.add_argument('--merge-ks', type=lambda s: [int(x) for x in s.split(',')], default=[5, 6, 7, 8, 10, 12, 14, 16, 20])
    ap.add_argument('--merge-pin', type=float, default=0.9); ap.add_argument('--merge-pout', type=float, default=0.1)
    ap.add_argument('--merge-min-diag', type=int, default=5); ap.add_argument('--merge-min-copies', type=int, default=20)
    ap.add_argument('--merge-sample', type=int, default=400)
    a = ap.parse_args()
    P = lambda f: os.path.join(a.out, f); dest = a.dest or P('53_segmentation'); os.makedirs(dest, exist_ok=True)
    arrays = list(csv.DictReader(open(P('02_satellite_arrays.tsv')), delimiter='\t'))
    if a.arrays:
        want = set(open(a.arrays).read().split()) if os.path.exists(a.arrays) else set(a.arrays.split(','))
        arrays = [x for x in arrays if x['array_id'] in want]
    if a.runs:                                   # run mode: each catalogue run is a unit; its copies are the parent array's copies inside it
        arrays = [dict(array_id=r['run_id'], parent=r['array_id'], chrom=r['chrom'], start_0based=r['start0'], end=r['end'], period_mode_bp=r['unit_bp'])
                  for r in csv.DictReader(open(a.runs), delimiter='	')]
    ids = {x.get('parent', x['array_id']) for x in arrays}
    copies = collections.defaultdict(list)
    with open(P('03_monomer_copies.tsv'), newline='') as fh:
        rd = csv.reader(fh, delimiter='\t'); h = next(rd); ix = {k: i for i, k in enumerate(h)}
        for row in rd:
            if row[ix['array_id']] not in ids: continue
            try: idt = float(row[ix['best_alignment_identity']])
            except ValueError: idt = 0.0
            comp = row[ix['complete_or_truncated']] == 'complete'
            s, e = int(row[ix['genomic_start_0based']]), int(row[ix['genomic_end']])
            conf = comp and idt >= 80
            copies[row[ix['array_id']]].append(dict(mid=(s + e) // 2, len=e - s, idt=idt, conf=conf, fam=row[ix['family_id']],
                                                    cls0=row[ix['monomer_class']], strand=row[ix['strand']],
                                                    seq=row[ix['sequence']].upper() if conf else ''))
    for v in copies.values(): v.sort(key=lambda c: c['mid'])
    # class merge rule (sensitivity sweep 2026-10-20): two step-03 classes stay separate only if >= --merge-min-diag k-mers are
    # near-fixed (in >= --merge-pin of one class's copies, <= --merge-pout of the other's) at some k in --merge-ks;
    # otherwise they are merged (union-find) before any label or topology feature is computed. Classes with < --merge-min-copies
    # confident copies are not testable and keep their own label (flag untested).
    mlog = open(os.path.join(dest, '53_class_merges.tsv'), 'w'); mw = csv.writer(mlog, delimiter='\t', lineterminator='\n')
    mw.writerow(['array_id', 'class_a', 'class_b', 'n_a', 'n_b'] + [f'diag_k{k}' for k in a.merge_ks] + ['max_diag', 'best_k', 'decision'])
    rng = np.random.default_rng(53)
    for aid, cs in copies.items():
        by = collections.defaultdict(list)
        for c in cs:
            if c['conf']: by[c['cls0']].append(c['seq'])
        testable = sorted((k for k, v in by.items() if len(v) >= a.merge_min_copies), key=lambda k: -len(by[k]))
        samp = {k: [by[k][i] for i in rng.choice(len(by[k]), min(len(by[k]), a.merge_sample), replace=False)] for k in testable}
        parent = {k: k for k in by}
        def find(x):
            while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
            return x
        kc = {k: {cl: collections.Counter(km for s in samp[cl] for km in {s[i:i + k] for i in range(len(s) - k + 1)}) for cl in testable} for k in a.merge_ks}
        for i, x in enumerate(testable):
            for y in testable[i + 1:]:
                nx, ny = len(samp[x]), len(samp[y]); dk = []
                for k in a.merge_ks:
                    cx, cy = kc[k][x], kc[k][y]
                    d = sum(1 for m, v in cx.items() if v >= a.merge_pin * nx and cy.get(m, 0) <= a.merge_pout * ny) + \
                        sum(1 for m, v in cy.items() if v >= a.merge_pin * ny and cx.get(m, 0) <= a.merge_pout * nx)
                    dk.append(d)
                sep = max(dk) >= a.merge_min_diag
                if not sep: parent[find(y)] = find(x)
                mw.writerow([aid, x, y, len(by[x]), len(by[y])] + dk + [max(dk), a.merge_ks[int(np.argmax(dk))], 'separate' if sep else 'merge'])
        for k in by:
            if k not in testable: mw.writerow([aid, k, 'NA', len(by[k]), 'NA'] + ['NA'] * len(a.merge_ks) + ['NA', 'NA', 'untested_small_class'])
        state = {k: '+'.join(sorted(m for m in by if find(m) == find(k))) for k in by}
        for c in cs:
            c['cls'] = state.get(c['cls0'], c['cls0']) + c['strand']; c['seq'] = ''
    mlog.close()
    print(f'{len(arrays)} arrays, {sum(map(len, copies.values()))} copies', file=sys.stderr)
    genome = read_fasta(a.genome, {x['chrom'] for x in arrays})
    WCOL = ['array_id', 'chrom', 'win', 'start0', 'end', 'valid_bp', 'GC', 'CpG_oe', 'ent1', 'ent2', 'ent3', 'ent5'] + \
           [f'kmer_pc{i + 1}' for i in range(a.kmer_pcs)] + ['local_period', 'period_strength', 'strength_at_02', 'recur_frac',
           'n_copies', 'len_med_rel', 'len_iqr_rel', 'mean_identity', 'confident_frac', 'dom_class_frac', 'n_classes', 'fam1_frac', 'fam2_frac', 'fam3_frac',
           'selfloop', 'trans_entropy', 'alternation', 'edge_diversity']
    BLOCKS = {'comp': ['GC', 'CpG_oe', 'ent1', 'ent2', 'ent3', 'ent5'], 'kmer': [f'kmer_pc{i + 1}' for i in range(a.kmer_pcs)],
              'period': ['local_period', 'period_strength', 'strength_at_02', 'recur_frac'],
              'label': ['len_med_rel', 'len_iqr_rel', 'mean_identity', 'confident_frac', 'dom_class_frac', 'n_classes', 'fam1_frac', 'fam2_frac', 'fam3_frac'],
              'topo': ['selfloop', 'trans_entropy', 'alternation', 'edge_diversity']}
    FREE = {'comp', 'kmer', 'period'}
    wf = gzip.open(os.path.join(dest, '53_window_features.tsv.gz'), 'wt'); ww = csv.writer(wf, delimiter='\t', lineterminator='\n'); ww.writerow(WCOL)
    cpf = open(os.path.join(dest, '53_change_points.tsv'), 'w'); cw = csv.writer(cpf, delimiter='\t', lineterminator='\n')
    cw.writerow(['array_id', 'chrom', 'position', 'window_index', 'blocks', 'n_blocks', 'label_free_support', 'status'])
    sgf = open(os.path.join(dest, '53_segments.tsv'), 'w'); sw = csv.writer(sgf, delimiter='\t', lineterminator='\n')
    sw.writerow(['array_id', 'chrom', 'segment', 'start0', 'end', 'bp', 'n_windows', 'profile', 'dominant_family', 'dominant_family_frac',
                 'n_copies', 'n_confident', 'class_selfloop', 'alternation', 'local_period_median', 'period_02'])
    bed = open(os.path.join(dest, '53_segments.bed'), 'w')
    asum = []; tot = collections.Counter(); totbp = collections.Counter()
    for k_, A in enumerate(arrays):
        aid, chrom = A['array_id'], A['chrom']; s0, e0 = int(A['start_0based']), int(A['end'])
        period = float(A['period_mode_bp'] or 0); cs = copies.get(A.get('parent', aid), [])
        if 'parent' in A: cs = [c for c in cs if s0 <= c['mid'] < e0]
        W = int(min(a.max_win_bp, max(a.min_win_bp, a.win_units * period))); st = max(1, W // 2)
        if e0 - s0 < W: starts = [s0]; W = e0 - s0
        else: starts = list(range(s0, e0 - W + 1, st))
        seq = genome.get(chrom, b'')
        mids = np.array([c['mid'] for c in cs]) if cs else np.empty(0)
        top3 = [f for f, _ in collections.Counter(c['fam'] for c in cs if c['conf']).most_common(3)]
        rows, profs = [], []
        for wi, ws in enumerate(starts):
            we = ws + W; arr = ENC[np.frombuffer(seq[ws:we], np.uint8)] if seq else np.empty(0, np.int8)
            f, prof = window_seq_features(arr, period)
            i0, i1 = np.searchsorted(mids, ws), np.searchsorted(mids, we); wc = cs[i0:i1]
            conf = [c for c in wc if c['conf']]
            if wc:
                L = np.array([c['len'] for c in wc]); q1, q3 = np.percentile(L, [25, 75])
                f.update(n_copies=len(wc), len_med_rel=float(np.median(L)) / max(1, period), len_iqr_rel=float(q3 - q1) / max(1, period),
                         mean_identity=float(np.mean([c['idt'] for c in wc])), confident_frac=len(conf) / len(wc))
                cc = collections.Counter(c['cls'] for c in conf)
                f['dom_class_frac'] = cc.most_common(1)[0][1] / len(conf) if conf else np.nan; f['n_classes'] = len(cc)
                fc = collections.Counter(c['fam'] for c in conf)
                for j in range(3): f[f'fam{j + 1}_frac'] = (fc[top3[j]] / len(conf) if conf and j < len(top3) else 0.0)
                sl, hh, alt, ed = topo([c['cls'] if c['conf'] else None for c in wc])
                f.update(selfloop=sl, trans_entropy=hh, alternation=alt, edge_diversity=ed)
            else: f['n_copies'] = 0
            rows.append((ws, we, f)); profs.append(prof if prof is not None else np.full(136, np.nan))
        P_ = np.array(profs)
        if len(rows) >= 3 and np.isfinite(P_).all(1).sum() >= 3:
            ok = np.isfinite(P_).all(1); M = P_[ok] - P_[ok].mean(0)
            U, Sv, Vt = np.linalg.svd(M, full_matrices=False); pcs = np.full((len(P_), a.kmer_pcs), np.nan)
            m = min(a.kmer_pcs, Vt.shape[0]); pcs[ok, :m] = M @ Vt[:m].T
        else: pcs = np.full((len(P_), a.kmer_pcs), np.nan)
        for wi, (ws, we, f) in enumerate(rows):
            for j in range(a.kmer_pcs): f[f'kmer_pc{j + 1}'] = float(pcs[wi, j]) if np.isfinite(pcs[wi, j]) else np.nan
            ww.writerow([aid, chrom, wi, ws, we] + [r(f.get(c, np.nan)) for c in WCOL[5:]])
        # segmentation
        cps = {}
        if len(rows) >= 2 * a.min_seg:
            for b, cols in BLOCKS.items():
                X = zs([[f.get(c, np.nan) for c in cols] for _, _, f in rows])
                if X.shape[1]: cps[b] = binseg(X, a.penalty, a.min_seg)
        allc = sorted({(t, b) for b, ts in cps.items() for t in ts})
        clusters = []
        for t, b in allc:
            if clusters and t - clusters[-1]['t'][-1] <= 1: clusters[-1]['t'].append(t); clusters[-1]['b'].add(b)
            else: clusters.append({'t': [t], 'b': {b}})
        acc = []
        for c in clusters:
            t = int(round(np.median(c['t']))); ok = len(c['b']) >= 2 and bool(c['b'] & FREE)
            pos = rows[t][0] + (rows[t][1] - rows[t][0]) // 2 if t < len(rows) else rows[-1][1]
            cw.writerow([aid, chrom, pos, t, ','.join(sorted(c['b'])), len(c['b']), ','.join(sorted(c['b'] & FREE)) or 'none', 'accepted' if ok else 'candidate'])
            if ok: acc.append(pos)
        bounds = [s0] + [p for p in acc if s0 < p < e0] + [e0]
        segs = []
        for si in range(len(bounds) - 1):
            ss, se = bounds[si], bounds[si + 1]
            sc = [c for c in cs if ss <= c['mid'] < se]
            prof, df, dff, nc, ncf, sl, alt = seg_profile(sc)
            lp = [f.get('local_period') for ws, we, f in rows if ws >= ss - W // 2 and we <= se + W // 2 and f.get('local_period')]
            nwin = sum(1 for ws, we, f in rows if ss <= (ws + we) // 2 < se)
            segs.append(prof)
            sw.writerow([aid, chrom, si + 1, ss, se, se - ss, nwin, prof, df, dff, nc, ncf, r(sl, 3), r(alt, 3), int(np.median(lp)) if lp else 'NA', period])
            bed.write(f'{chrom}\t{ss}\t{se}\t{aid}|seg{si + 1}|{prof}|{df}\t0\t.\t{ss}\t{se}\t{COLOR[prof]}\n')
            tot[prof] += 1; totbp[prof] += se - ss
        arch = 'single_segment' if len(segs) == 1 else ('same_profile_segments' if len(set(segs)) == 1 else 'multi_profile')
        asum.append([aid, chrom, s0, e0, period, W, len(rows), len(cs), len(acc), sum(len(v) for v in cps.values()), len(segs), arch, ','.join(segs)])
        if k_ % 200 == 0: print(f'{k_}/{len(arrays)} {aid}', file=sys.stderr)
    wf.close(); cpf.close(); sgf.close(); bed.close()
    with open(os.path.join(dest, '53_array_summary.tsv'), 'w') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['array_id', 'chrom', 'start0', 'end', 'period_02', 'window_bp', 'n_windows', 'n_copies', 'n_accepted_boundaries',
                    'n_block_change_points', 'n_segments', 'architecture', 'segment_profiles']); w.writerows(asum)
    arch = collections.Counter(x[11] for x in asum)
    with open(os.path.join(dest, '53_summary.tsv'), 'w') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n'); w.writerow(['kind', 'name', 'n', 'Mb'])
        for p, n in tot.most_common(): w.writerow(['segment_profile', p, n, round(totbp[p] / 1e6, 2)])
        for p, n in arch.most_common(): w.writerow(['array_architecture', p, n, ''])
    for l in open(os.path.join(dest, '53_summary.tsv')): print(l.rstrip())


if __name__ == '__main__':
    main()
