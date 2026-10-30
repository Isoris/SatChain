#!/usr/bin/env python3
"""
64_object_profiles.py  (SatChain v4.7 draft; one sequence profile per object, then separation tests for any grouping)

Objects (all optional, combined):
  copy    step-03 copies inside the selected catalogue runs (identity >= --min-identity); attributes: class, strand, identity,
          index from the run start / end (array orientation not applied: genomic order)
  run     the selected catalogue runs (satellite blocks)
  gap     interruptions between consecutive selected runs of one array (the spacers), gap <= --max-gap
  extra   any BED (name:path, repeatable via --extra): e.g. RTE fragments, flanks, inserts; column 4 = object id, column 6 = strand
Profile per object: length, GC, CpG obs/exp, GC skew, AT skew, normalised Shannon entropy of k-mers k = 1..6, and the canonical
k-mer frequency vector for every k in --ks (stored sparse in 64_kmer_<k>.npz with the object order of 64_objects.tsv).
Separation (--group col[,col...]; columns of 64_objects.tsv or of a --meta table joined on chrom / array_id / run_id):
  for every object type x grouping x k: mean Jensen-Shannon distance within groups vs between groups (sample <= --sample objects per
  group), permutation pseudo-F (--perm label shuffles), leave-one-out nearest-centroid accuracy vs chance (largest group share).
Outputs (<dest>/): 64_objects.tsv, 64_kmer_<k>.npz, 64_separation.tsv.
usage: 64_object_profiles.py --catalogue <59> --copies <03> [...] --genome <fa> --dest <dir> [--family F | --runs-file ids]
       [--extra name:bed ...] [--meta meta.tsv] [--group type,chrom] [--ks 2,3,4,5,6,7,8] [--sample 60] [--perm 199]
"""
import argparse, collections, csv, math, os, random
import numpy as np


def save_sparse(path, r, c, d, shape): np.savez_compressed(path, r=np.array(r, np.int32), c=np.array(c, np.int32), d=np.array(d, np.float32), shape=np.array(shape))


def load_rows(path, rows):
    z = np.load(path); M = np.zeros((len(rows), int(z['shape'][1])), np.float32); pos = {r: i for i, r in enumerate(rows)}
    sel = np.isin(z['r'], rows)
    for r, c, d in zip(z['r'][sel], z['c'][sel], z['d'][sel]): M[pos[r], c] += d
    return M
csv.field_size_limit(1 << 30)
RC = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(RC)[::-1]


class Fasta:
    def __init__(self, path):
        self.fh = open(path, 'rb'); self.idx = {}
        for l in open(path + '.fai'):
            n, ln, off, lb, lby = l.split('\t')[:5]; self.idx[n] = (int(ln), int(off), int(lb), int(lby))
    def get(self, c, a, b):
        ln, off, lb, lby = self.idx[c]; a = max(0, a); b = min(ln, b)
        s0 = off + (a // lb) * lby + a % lb; s1 = off + (b // lb) * lby + b % lb
        self.fh.seek(s0); return self.fh.read(s1 - s0).replace(b'\n', b'').replace(b'\r', b'').decode().upper()


def canon(x): y = rc(x); return x if x < y else y


def profile(s, ks):
    n = len(s); acgt = sum(s.count(b) for b in 'ACGT') or 1
    g, c, a, t = (s.count(b) for b in 'GCAT')
    cpg = s.count('CG'); exp = (c * g) / max(n, 1)
    f = dict(length=n, GC=(g + c) / acgt, CpG_oe=(cpg / exp) if exp else 0.0, GC_skew=(g - c) / max(g + c, 1), AT_skew=(a - t) / max(a + t, 1))
    for k in range(1, 7):
        cnt = collections.Counter(s[i:i + k] for i in range(n - k + 1) if 'N' not in s[i:i + k]); tot = sum(cnt.values())
        H = -sum(v / tot * math.log(v / tot) for v in cnt.values()) if tot else 0.0
        f[f'entropy_k{k}'] = H / math.log(min(4 ** k, max(tot, 2)))
    vec = {}
    for k in ks:
        vec[k] = collections.Counter(canon(s[i:i + k]) for i in range(n - k + 1) if 'N' not in s[i:i + k])
    return f, vec


def jsd(p, q):
    m = 0.5 * (p + q); out = 0.0
    for x in (p, q):
        nz = x > 0; out += 0.5 * np.sum(x[nz] * np.log2(x[nz] / m[nz]))
    return math.sqrt(max(out, 0.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--catalogue', required=True); ap.add_argument('--copies', nargs='*', default=[]); ap.add_argument('--genome', required=True)
    ap.add_argument('--dest', required=True); ap.add_argument('--family', default=''); ap.add_argument('--runs-file', default='')
    ap.add_argument('--extra', action='append', default=[]); ap.add_argument('--meta', default=''); ap.add_argument('--group', default='')
    ap.add_argument('--ks', default='2,3,4,5,6,7,8'); ap.add_argument('--min-identity', type=float, default=60); ap.add_argument('--max-gap', type=int, default=50000)
    ap.add_argument('--min-runs', type=int, default=1); ap.add_argument('--sample', type=int, default=60); ap.add_argument('--perm', type=int, default=199)
    ap.add_argument('--seed', type=int, default=64); ap.add_argument('--loo-block', default='run_id,array_id', help='per grouping order: leave out every object sharing this column with the tested one')
    a = ap.parse_args(); rng = random.Random(a.seed); os.makedirs(a.dest, exist_ok=True); KS = [int(x) for x in a.ks.split(',')]
    fa = Fasta(a.genome)
    want = set(open(a.runs_file).read().split()) if a.runs_file else None
    runs = [r for r in csv.DictReader(open(a.catalogue), delimiter='\t') if (not a.family or r['family'] == a.family) and (want is None or r['run_id'] in want)]
    by = collections.defaultdict(list)
    for r in runs: by[r['array_id']].append(r)
    by = {k: sorted(v, key=lambda r: int(r['start0'])) for k, v in by.items() if len(v) >= a.min_runs}
    objs = []                                                    # (type, id, chrom, start, end, strand, attrs)
    for aid, rs in by.items():
        for i, r in enumerate(rs):
            objs.append(('run', r['run_id'], r['chrom'], int(r['start0']), int(r['end']), '+', dict(array_id=aid, run_id=r['run_id'], run_index=i)))
        for i, (r1, r2) in enumerate(zip(rs, rs[1:])):
            s, e = int(r1['end']), int(r2['start0'])
            if 0 < e - s <= a.max_gap:
                objs.append(('gap', f'{aid}#gap{i}', r1['chrom'], s, e, '+', dict(array_id=aid, run_id=r1['run_id'], run_index=i)))
    if a.copies:
        span = collections.defaultdict(list)
        for aid, rs in by.items():
            for r in rs: span[r['chrom']].append((int(r['start0']), int(r['end']), r['run_id'], aid))
        per_run = collections.defaultdict(list)
        for p in a.copies:
            with open(p, newline='') as fh:
                rd = csv.reader(fh, delimiter='\t'); h = next(rd); ix = {k: i for i, k in enumerate(h)}
                for row in rd:
                    c = row[ix['chrom']]
                    if c not in span or float(row[ix['best_alignment_identity']] or 0) < a.min_identity: continue
                    s, e = int(row[ix['genomic_start_0based']]), int(row[ix['genomic_end']])
                    for rs_, re_, rid, aid in span[c]:
                        if s >= rs_ and e <= re_:
                            per_run[rid].append((s, e, c, aid, row[ix['monomer_class']], row[ix['strand']], float(row[ix['best_alignment_identity']])))
                            break
        for rid, L in per_run.items():
            L.sort(); n = len(L)
            for j, (s, e, c, aid, cl, st, idt) in enumerate(L):
                objs.append(('copy', f'{rid}#c{j}', c, s, e, st, dict(array_id=aid, run_id=rid, run_index='', copy_class=cl, identity=idt, index_from_start=j, index_from_end=n - 1 - j)))
    for x in a.extra:
        name, path = x.split(':', 1)
        for l in open(path):
            f = l.rstrip('\n').split('\t')
            objs.append((name, f[3], f[0], int(f[1]), int(f[2]), f[5] if len(f) > 5 else '+', dict(array_id='', run_id='', run_index='')))
    rows, mats = [], {k: [] for k in KS}; vocab = {k: {} for k in KS}
    for t, oid, c, s, e, st, at in objs:
        seq = fa.get(c, s, e); seq = rc(seq) if st == '-' else seq
        f, vec = profile(seq, KS)
        rows.append(dict(object_type=t, object_id=oid, chrom=c, start=s, end=e, strand=st, **at, **f))
        for k in KS:
            v = vocab[k]; mats[k].append({v.setdefault(x, len(v)): n for x, n in vec[k].items()})
    cols = sorted({k for r in rows for k in r}, key=lambda k: list(rows[0].keys()).index(k) if k in rows[0] else 99)
    with open(os.path.join(a.dest, '64_objects.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore'); w.writeheader(); w.writerows(rows)
    for k in KS:
        r_, c_, d_ = [], [], []
        for i, m in enumerate(mats[k]):
            for j, n in m.items(): r_.append(i); c_.append(j); d_.append(n)
        save_sparse(os.path.join(a.dest, f'64_kmer_{k}.npz'), r_, c_, d_, (len(rows), len(vocab[k])))
    print(f'objects {len(rows)}:', dict(collections.Counter(r['object_type'] for r in rows)), flush=True)
    if not a.group: return
    meta = {}
    if a.meta:
        for r in csv.DictReader(open(a.meta), delimiter='\t'): meta[r[next(iter(r))]] = r
    key0 = next(iter(next(iter(meta.values())))) if meta else None
    out = []
    for grp in a.group.split(','):
        for t in sorted({r['object_type'] for r in rows}):
            idx = [i for i, r in enumerate(rows) if r['object_type'] == t]
            lab = []
            for i in idx:
                r = rows[i]; v = r.get(grp)
                if v in (None, '') and meta: v = meta.get(r.get(key0, ''), {}).get(grp)
                lab.append(v)
            keep = [(i, l) for i, l in zip(idx, lab) if l not in (None, '')]
            groups = collections.defaultdict(list)
            for i, l in keep: groups[l].append(i)
            groups = {g: rng.sample(v, min(a.sample, len(v))) for g, v in groups.items() if len(v) >= 3}
            if len(groups) < 2: continue
            ids = [i for v in groups.values() for i in v]; labs = np.array([g for g, v in groups.items() for _ in v])
            for k in KS:
                M = load_rows(os.path.join(a.dest, f'64_kmer_{k}.npz'), ids); M = M / np.maximum(M.sum(1, keepdims=True), 1)
                n = len(ids); D = np.zeros((n, n))
                for x in range(n):
                    for y in range(x + 1, n): D[x, y] = D[y, x] = jsd(M[x], M[y])
                def pseudoF(L):
                    W = sum(D[np.ix_(L == g, L == g)].sum() / (2 * max((L == g).sum(), 1)) for g in set(L)); T = D.sum() / (2 * n)
                    G = len(set(L)); return ((T - W) / max(G - 1, 1)) / (W / max(n - G, 1)) if W > 0 else float('inf')
                F0 = pseudoF(labs); null = [pseudoF(np.random.default_rng(a.seed + p).permutation(labs)) for p in range(a.perm)]
                same = D[(labs[:, None] == labs[None, :]) & ~np.eye(n, dtype=bool)].mean(); diff = D[labs[:, None] != labs[None, :]].mean()
                blk_col = (a.loo_block.split(',') + [a.loo_block.split(',')[-1]] * 9)[a.group.split(',').index(grp)]
                blk = np.array([str(rows[i].get(blk_col, '')) or str(i) for i in ids]); ok = 0
                for x in range(n):
                    keep_ = blk != blk[x]
                    cands = {g: M[keep_ & (labs == g)].mean(0) for g in set(labs) if (keep_ & (labs == g)).any()}
                    if not cands: continue
                    best = min(cands, key=lambda g: jsd(M[x], cands[g])); ok += best == labs[x]
                cent = cands
                out.append(dict(grouping=grp, object_type=t, k=k, groups=len(cent), objects=n, jsd_within=round(same, 4), jsd_between=round(diff, 4),
                                pseudoF=round(F0, 3), p_perm=round((1 + sum(f >= F0 for f in null)) / (1 + a.perm), 4),
                                loo_accuracy=round(ok / n, 3), loo_block=blk_col, chance=round(max(collections.Counter(labs).values()) / n, 3)))
                print(grp, t, k, out[-1]['jsd_within'], out[-1]['jsd_between'], out[-1]['loo_accuracy'], out[-1]['p_perm'], flush=True)
    if out:
        with open(os.path.join(a.dest, '64_separation.tsv'), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(out)


if __name__ == '__main__':
    main()
