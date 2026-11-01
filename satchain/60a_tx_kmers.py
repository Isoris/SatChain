#!/usr/bin/env python3
"""
60a_tx_kmers.py  (v4.6; k-mer SETS for the satellite transcription test, one subgenome)

60a builds canonical 31-mer sets (min of the forward and reverse-complement code), not a transition model: a reverse-strand unit needs
no handling, and the FASTA writes the canonical form, so the downstream RNA-seq count MUST be canonical too (jellyfish -C or
equivalent); a stranded count against this set loses about half the hits (D15).
Groups:
  satellite families   catalogue satellite families with >= --min-bp in runs (name_v46); canonical 31-mers from their runs (up to
                       --max-seq bp sampled per family), kept when seen >= 2 times inside the family's runs
  gene_array controls  catalogue gene-array families (5S, U1, U2, U6, tRNA ...) the same way
  5S                   31-mers of barrnap 5S genes
  actb                 31-mers of the exons of the actb gene(s) (LiftOn) = housekeeping positive control
  intergenic           31-mers from random positions outside every annotation available (EDTA TEanno, LiftOn genes, catalogue runs,
                       --trf BED), genome count 1 = DNA-contamination / background control (2026-11-01, D15: before, only genome count 1)
Specificity: every query k-mer is counted genome-wide (canonical, exact); a satellite k-mer is kept if its genome count outside its
own family's runs is 0 (strict) or <= 10 % of its total (relaxed; both written); 'inside' is counted on the UNION of the family's runs
(overlapping runs merged, 2026-11-01 D15; before, overlapping runs were counted twice and out = total - inside could reach <= 0).
k-mers present in the 5S genes or in the EDTA TE libraries of both subgenomes are removed from satellite sets.
Paths (2026-11-01, D08): annotation files are resolved with satpaths.func_path (current fClaHyb_<G>_<sub> layout, then the old _<sub>
layout); the other subgenome's annotation folder is given explicitly with --other-func (default: <func> with G swapped, warned).
Outputs (<dest>/): tx_kmers_<G>.fa (one record per k-mer, id = group|tier|k), tx_kmer_table_<G>.tsv (k-mer, group, genome count,
count inside group, tier)
usage: 60a_tx_kmers.py <G> <cat59 dir> <func folder> <genome fasta> <dest> [--other-func DIR] [--min-bp 50000] [--max-seq 400000]
                       [--trf BED]
"""
import sys, os, csv, gzip, collections, random, re, bisect, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from satpaths import func_path
ap = argparse.ArgumentParser()
for x in ('G', 'C', 'FUNC', 'FA', 'DEST'): ap.add_argument(x)
ap.add_argument('--other-func', default=''); ap.add_argument('--min-bp', type=int, default=50000); ap.add_argument('--max-seq', type=int, default=400000)
ap.add_argument('--trf', default='', help='TRF BED of the subgenome (optional background filter)')
A = ap.parse_args(); G, C, FUNC, FA, DEST = A.G, A.C, A.FUNC, A.FA, A.DEST; MINBP = A.min_bp; MAXSEQ = A.max_seq; K = 31; N_PER = 4000
G2 = 'Mac' if G == 'Gar' else 'Gar'; OTHER = A.other_func or FUNC.replace(G, G2)
if not A.other_func: print(f'--other-func not given: using {OTHER}')
os.makedirs(DEST, exist_ok=True); rng = random.Random(60)
CODE = np.full(256, 4, np.int64)
for i, c in enumerate(b'ACGT'): CODE[c] = i; CODE[c + 32] = i
fai = {l.split()[0]: list(map(int, l.split()[1:5])) for l in open(FA + '.fai')}; fh_ = open(FA, 'rb')
def fetch(c, s, e):
    L, off, lb, lB = fai[c]; s = max(0, s); e = min(L, e); fh_.seek(off + (s // lb) * lB + s % lb)
    return fh_.read(int((e - s) * lB / lb) + 2 * lB).replace(b'\n', b'').replace(b'\r', b'')[:e - s]
MASK = (1 << (2 * K)) - 1
def kcodes(b):
    c = CODE[np.frombuffer(b, dtype=np.uint8)]; n = len(c) - K + 1
    if n <= 0: return np.empty(0, np.int64)
    f = np.zeros(n, np.int64); r = np.zeros(n, np.int64); bad = np.zeros(n, bool)
    for j in range(K):
        x = c[j:j + n]; bad |= x >= 4; x = np.minimum(x, 3)
        f = (f << 2) | x; r |= (3 - x) << (2 * j)
    can = np.minimum(f, r); return can[~bad]
def merge(iv):
    """union of [s, e) intervals"""
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return out
# ---- groups from the catalogue
runs = list(csv.DictReader(open(f'{C}/59_catalogue_{G}.tsv'), delimiter='\t'))
fam_runs = collections.defaultdict(list); fam_bp = collections.Counter(); fam_kind = {}
for r in runs:
    key = r['name_v46'] if r['element_class'] == 'satellite' else r['element_class'] + '|' + r['family_v46']
    fam_runs[key].append((r['chrom'], int(r['start0']), int(r['end']))); fam_bp[key] += int(r['bp'])
    fam_kind[key] = 'satellite' if r['element_class'] == 'satellite' else 'gene_array'
groups = [k for k, b in fam_bp.most_common() if b >= MINBP and (fam_kind[k] == 'satellite' or re.search(r'5S|U1|U2|U6|tRNA', k))]
q = {}            # kmer -> group
inside = collections.Counter()
for gname in groups:
    rr = fam_runs[gname][:]; rng.shuffle(rr); got = 0; cnt = collections.Counter()
    for c, s, e in rr:
        if got >= MAXSEQ: break
        b = fetch(c, s, min(e, s + MAXSEQ - got)); got += len(b)
        for x in kcodes(b).tolist(): cnt[x] += 1
    keep = [x for x, n in cnt.items() if n >= 2 and x not in q]
    rng.shuffle(keep)
    for x in keep[:N_PER]: q[x] = gname
    print(gname, fam_kind[gname], fam_bp[gname], 'kmers kept', min(len(keep), N_PER), flush=True)
# exact in-family counts of the selected k-mers over the UNION of the family's runs (D15)
sel = collections.defaultdict(set)
for x, gname in q.items(): sel[gname].add(x)
for gname in groups:
    qs = np.array(sorted(sel[gname]), np.int64); byc = collections.defaultdict(list)
    for c, s, e in fam_runs[gname]: byc[c].append((s, e))
    for c, iv in byc.items():
        for s, e in merge(iv):
            v = kcodes(fetch(c, s, e)); i = np.searchsorted(qs, v); ok = (i < len(qs)) & (qs[np.minimum(i, len(qs) - 1)] == v)
            for x in v[ok].tolist(): inside[x] += 1
# controls: 5S genes, actb exons, intergenic (paths via satpaths, D08)
ctrl = {}; annot = collections.defaultdict(list)          # annot: chrom -> intervals of every annotation (background filter)
for l in gzip.open(func_path(FUNC, G, '_rRNA_barrnap', f'fClaHyb_{G}.rRNA.sorted.gff3.gz'), 'rt'):
    if l.startswith('#') or '5S' not in l: continue
    f = l.split('\t')
    for x in kcodes(fetch(f[0], int(f[3]) - 1, int(f[4]))).tolist(): ctrl[x] = 'control_5S'
actb = 0
for l in gzip.open(func_path(FUNC, G, '_genes_lifton', f'fClaHyb_{G}.lifton.sorted.gff3.gz'), 'rt'):
    if l.startswith('#'): continue
    f = l.split('\t')
    if f[2] == 'gene': annot[f[0]].append((int(f[3]) - 1, int(f[4])))
    if f[2] == 'exon' and re.search(r'(gene|Name|gene_name)=actb[a-z0-9]*[;,\n]', f[8]):
        actb += 1
        for x in kcodes(fetch(f[0], int(f[3]) - 1, int(f[4]))).tolist(): ctrl.setdefault(x, 'control_actb')
print('actb exons', actb)
te_gff = func_path(FUNC, G, '_repeats_EDTA', f'fClaHyb_{G}.TEanno.sorted.gff3.gz', required=False)
if te_gff:
    for l in gzip.open(te_gff, 'rt'):
        if l.startswith('#'): continue
        f = l.split('\t'); annot[f[0]].append((int(f[3]) - 1, int(f[4])))
else: print('no EDTA TEanno GFF: background not filtered against TEs')
for r in runs: annot[r['chrom']].append((int(r['start0']), int(r['end'])))
if A.trf:
    for l in open(A.trf):
        f = l.split('\t')
        if len(f) >= 3 and not l.startswith(('#', 'track')): annot[f[0]].append((int(f[1]), int(f[2])))
annot = {c: merge(iv) for c, iv in annot.items()}; astart = {c: [s for s, _ in iv] for c, iv in annot.items()}
def annotated(c, s, e):
    iv = annot.get(c, [])
    if not iv: return False
    i = bisect.bisect_right(astart[c], e) - 1
    return i >= 0 and iv[i][1] > s
# TE library + 5S k-mers to remove from satellite sets
excl = set(x for x, g in ctrl.items() if g == 'control_5S')
for p in (func_path(FUNC, G, '_repeats_EDTA', f'fClaHyb_{G}.TElib.fa', required=False), func_path(OTHER, G2, '_repeats_EDTA', f'fClaHyb_{G2}.TElib.fa', required=False)):
    if not p: print('missing TE library (one subgenome)'); continue
    seq = []
    for l in open(p, 'rb'):
        if l.startswith(b'>'):
            if seq: excl.update(kcodes(b''.join(seq)).tolist())
            seq = []
        else: seq.append(l.strip())
    if seq: excl.update(kcodes(b''.join(seq)).tolist())
nrm = sum(1 for x in q if x in excl); q = {x: g for x, g in q.items() if x not in excl}
print('removed (TE library / 5S):', nrm)
# intergenic background: random annotation-free 31-mers, checked later for genome count 1
chroms = list(fai); lens = [fai[c][0] for c in chroms]; bg = {}; tried = 0
while len(bg) < 20000 and tried < 2_000_000:
    tried += 1; c = rng.choices(chroms, weights=lens)[0]; s = rng.randrange(0, fai[c][0] - K)
    if annotated(c, s, s + K): continue
    v = kcodes(fetch(c, s, s + K))
    if len(v): bg[int(v[0])] = (c, s)
print(f'background positions: {len(bg)} annotation-free of {tried} drawn')
# genome-wide exact counts of every query k-mer
allq = np.array(sorted(set(q) | set(ctrl) | set(bg)), np.int64); tot = np.zeros(len(allq), np.int64)
for c in chroms:
    L = fai[c][0]
    for s in range(0, L, 20_000_000):
        v = kcodes(fetch(c, s, min(L, s + 20_000_000 + K - 1))); i = np.searchsorted(allq, v); ok = (i < len(allq)) & (allq[np.minimum(i, len(allq) - 1)] == v)
        np.add.at(tot, i[ok], 1)
    print('counted', c, flush=True)
T = dict(zip(allq.tolist(), tot.tolist()))
with open(f'{DEST}/tx_kmer_table_{G}.tsv', 'w') as tf, open(f'{DEST}/tx_kmers_{G}.fa', 'w') as ff:
    tf.write('kmer_code\tgroup\tkind\tgenome_count\tinside_group\ttier\n'); n = collections.Counter()
    def dec(x): return ''.join('ACGT'[(x >> (2 * (K - 1 - j))) & 3] for j in range(K))
    for x, gname in q.items():
        t = T[x]; ins = inside[x]; out = t - ins
        tier = 'strict' if out <= 0 else 'relaxed' if out <= 0.1 * t else 'not_specific'
        tf.write(f'{x}\t{gname}\t{fam_kind[gname]}\t{t}\t{ins}\t{tier}\n')
        if tier != 'not_specific': ff.write(f'>{gname}|{tier}|{x}\n{dec(x)}\n'); n[(gname, tier)] += 1
    for x, gname in ctrl.items():
        tf.write(f'{x}\t{gname}\tcontrol\t{T[x]}\tNA\tcontrol\n'); ff.write(f'>{gname}|control|{x}\n{dec(x)}\n')
    kept_bg = 0
    for x, (c, s) in bg.items():
        if T[x] == 1 and kept_bg < 5000:
            kept_bg += 1; tf.write(f'{x}\tcontrol_intergenic\tcontrol\t1\tNA\tcontrol\n'); ff.write(f'>control_intergenic|control|{x}\n{dec(x)}\n')
print('intergenic background k-mers (genome count 1, annotation-free):', kept_bg)
for gname in groups: print(gname, 'strict', n[(gname, 'strict')], 'relaxed', n[(gname, 'relaxed')])
