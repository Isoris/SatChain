#!/usr/bin/env python3
"""
59_final_catalogue.py  (v4.7, 2026-11-01; the validated satellite catalogue of one subgenome)

Runs = the copy-chain HMM runs that passed the per-run sequence check (56 -> val56, q0.9 set), existing step-02 arrays and new
scan-57 arrays alike (56 now fits one model per haplotype on both tables; a legacy separate bed for the new arrays is still accepted).
Per run, from its step-03 copies (exact duplicate tiles dropped, audit CB01):
  n_conf (complete, identity >= 80), n_weak (65-80), n_junk (< 65); mean identity of confident copies; dominant class
  G1 dominant-class consensus GC within 5 points of the run's own GC; G2 mean identity of confident copies >= 85 (2026-11-01: the old
  75/85 split could not fail for units > 60 bp because confident copies are >= 80 by definition, audit C06); G3 strand coherent
  (reverse fraction outside 0.35-0.65). The gate cut-offs are swept into 59_gates_sweep_<G>.tsv (G1 3/5/8, G2 80/85/90, G3 0.3-0.7).
  family = majority family of confident copies (existing arrays) or the family matched by the new array's class consensus (famnew;
  'new:<array>' otherwise)
  unit class by the 100-bp rule of the runs of record (<= 6 micro, 7-99 mini, >= 100 satellite; audit C11)
  EDTA overlap fraction on the clipped union of features (<= 1; NOT TE embedding: EDTA labels satellites as TE)
Overlapping runs (existing vs new source) are resolved before anything is counted (audit C08): same family -> one merged run, different
families -> the run with fewer confident copies is trimmed to the other's edge; column overlap_note records it.
Gene arrays are called HERE, before tiers and the centromere flag (audit C04): 59b_gene_arrays.classify() gives element_class =
gene_array:<class> (observed/expected bp ratio >= 0.7 AND >= 0.5 of the units carry a hit), gene_related:<class> (ratio >= 0.2) or
satellite; gene_flag = <class> for both gene_array and gene_related families (audit C25).
Tiers by cumulative share of the SATELLITE bp (gene arrays excluded; tier 'gene_array'), families ranked by bp: core (first 50 %),
major (to 90 %), minor (to 99 %), rare. 59d recomputes them after the NewSat grouping.
Sequence cluster (was 'lineage'; 2026-11-01: the k-mer clustering of step 55 is a composition cluster of monomers, not an evolutionary
lineage): column seq_cluster from the 55d table when given ('none' otherwise); 59f_apply_clusters.py can add it later.
Centromere location (L1 only): putative_centromeric_array_L1 = satellite run overlapping the step-31 candidate interval (as before;
'na' for gene arrays); cen_candidate_ror = the run belongs to the largest cluster of satellite runs (<= 100 kb apart) of its chromosome,
derived from the runs of record themselves (audit C07).
Runs without a confident copy are logged in 59_dropped_runs_<G>.tsv (audit C22); an array missing from both 02 tables raises.
Outputs: 59_catalogue_<G>.tsv, 59_catalogue_<G>.bed (BED9 by seq_cluster), 59_family_tiers_<G>.tsv, 59_summary_<G>.tsv,
59b_gene_arrays_<G>.tsv (+ _sweep), 59_gates_sweep_<G>.tsv, 59_dropped_runs_<G>.tsv
usage: 59_final_catalogue.py <G> <pipeline out> <hmm56 dir> <new03 dir|-> <hmm56new dir|-> <famnew tsv|-> <cluster tsv|-> <func folder> <genome fasta> <dest>
       [--overlays <annotation_overlays/G>] [--carry-new <draft catalogue>] [--validated-name 56_runs_validated_q0.9.bed]
--carry-new: rows of source 'new' are taken from an earlier catalogue when the new-array 03 table is lost (Mac, 2026-10-31); they are
flagged carried_from and get no copy-level recomputation.
"""
import sys, os, csv, collections, bisect, gzip, argparse, importlib.util
csv.field_size_limit(1 << 30)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from satpaths import func_path
ap = argparse.ArgumentParser()
for k in ('G', 'OUT', 'H', 'N03', 'HN', 'FAMNEW', 'CLUS', 'FUNC', 'FA', 'DEST'): ap.add_argument(k)
ap.add_argument('--overlays', default=''); ap.add_argument('--carry-new', default=''); ap.add_argument('--validated-name', default='56_runs_validated_q0.9.bed')
ap.add_argument('--gene-cut', type=float, default=0.7); ap.add_argument('--unit-cut', type=float, default=0.5)
a = ap.parse_args(); G, OUT, H, N03, HN, FAMNEW, CLUS, FUNC, FA, DEST = [getattr(a, k) for k in ('G', 'OUT', 'H', 'N03', 'HN', 'FAMNEW', 'CLUS', 'FUNC', 'FA', 'DEST')]
os.makedirs(DEST, exist_ok=True); P = os.path.join(OUT, f'fClaHyb_{G}_satellites'); N03 = '' if N03 == '-' else N03; HN = '' if HN == '-' else HN
OVL = a.overlays or os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(FUNC))), 'satellite_FISH', 'annotation_overlays', G)
clus = {}
if CLUS and CLUS != '-':
    for r in csv.DictReader(open(CLUS), delimiter='\t'):
        if r['subgenome'] == G: clus[r['family']] = r.get('seq_cluster', r.get('lineage', 'none'))
famnew = {}
if FAMNEW and FAMNEW != '-' and os.path.exists(FAMNEW):
    for r in csv.DictReader(open(FAMNEW), delimiter='\t'): famnew[(r['array_id'], r['monomer_class'])] = r['family']


def load_copies(path, keep):
    cop = collections.defaultdict(list); ndup = 0
    with open(path, newline='') as fh:
        rd = csv.reader(fh, delimiter='\t'); h = next(rd); ix = {k: i for i, k in enumerate(h)}
        for row in rd:
            if row[ix['array_id']] not in keep: continue
            try: idt = float(row[ix['best_alignment_identity']])
            except ValueError: idt = 0.0
            cop[row[ix['array_id']]].append((int(row[ix['genomic_start_0based']]), int(row[ix['genomic_end']]), idt,
                                             row[ix['complete_or_truncated']] == 'complete', row[ix['strand']], row[ix['family_id']], row[ix['monomer_class']]))
    for k, v in cop.items():
        v.sort(); u = []; prev = None
        for x in v:
            if prev is not None and x[:2] == prev: ndup += 1; continue                      # CB01: duplicate tile
            u.append(x); prev = x[:2]
        cop[k] = u
    print(f'{path}: {sum(len(v) for v in cop.values())} copies of {len(cop)} arrays, {ndup} duplicate tiles dropped', flush=True)
    return cop


def runs(bed):
    out = []
    for l in open(bed):
        if l.startswith(('#', 'track')): continue
        c, s, e, nm = l.split('\t')[:4]; seg = [x[1:] for x in nm.split('|') if x.startswith('g') and x[1:].isdigit()]
        out.append((c, int(s), int(e), nm.split('|')[0], int(seg[0]) if seg else 1))
    return out


per = {}; per_old = {}; per_new = {}
for path, d in ((os.path.join(P, '02_satellite_arrays.tsv'), per_old), (os.path.join(N03, '02_satellite_arrays.tsv') if N03 else '', per_new)):
    if path and os.path.exists(path):
        for r in csv.DictReader(open(path), delimiter='\t'): d[r['array_id']] = float(r['period_mode_bp'] or 0)
per.update(per_old); per.update(per_new)
R_all = runs(os.path.join(H, a.validated_name))
R_new_legacy = runs(os.path.join(HN, '56_runs_post0.95.bed')) if HN and os.path.exists(os.path.join(HN, '56_runs_post0.95.bed')) else []
R_old = [r for r in R_all if r[3] in per_old]; R_new = [r for r in R_all if r[3] in per_new] + R_new_legacy
unknown = sorted({r[3] for r in R_all + R_new_legacy} - set(per))
if unknown: raise SystemExit(f'{len(unknown)} run arrays missing from both 02 tables, e.g. {unknown[:5]}')
C_old = load_copies(os.path.join(P, '03_monomer_copies.tsv'), {r[3] for r in R_old})
C_new = load_copies(os.path.join(N03, '03_monomer_copies.tsv'), {r[3] for r in R_new}) if R_new and N03 else {}
NSEG = collections.Counter(); [NSEG.update([(r[3], r[4])]) for r in R_old + R_new]
NSEG = collections.Counter(aid for aid, _ in NSEG)         # arrays that 56 split at assembly gaps get array_id <array>~g<segment>
cgc = {}
for path in (os.path.join(P, '03_array_consensus.tsv'), os.path.join(N03, '03_array_consensus.tsv') if N03 else ''):
    if path and os.path.exists(path):
        for r in csv.DictReader(open(path), delimiter='\t'): cgc[(r['array_id'], r['monomer_class'])] = float(r['consensus_GC_percent'])
te = collections.defaultdict(list)
for l in gzip.open(func_path(FUNC, G, '_repeats_EDTA', f'fClaHyb_{G}.TEanno.sorted.gff3.gz'), 'rt'):
    if l.startswith('#'): continue
    f = l.split('\t', 5)
    if len(f) > 4 and f[2] not in ('repeat_region', 'target_site_duplication', 'long_terminal_repeat'): te[f[0]].append((int(f[3]) - 1, int(f[4])))
for c in te:                                                # clipped union (C11): nested EDTA features counted once
    v = sorted(te[c]); m = []
    for s, e in v:
        if m and s <= m[-1][1]: m[-1] = (m[-1][0], max(m[-1][1], e))
        else: m.append((s, e))
    te[c] = m
cen = collections.defaultdict(list)
for r in csv.DictReader(open(os.path.join(P, '31_centromere_summary.tsv')), delimiter='\t'): cen[r['chrom']].append((int(r['cen_start_0based']), int(r['cen_end'])))


def ov(d, c, s, e):
    v = d.get(c, []); i = bisect.bisect_left(v, (s - 2000000, 0)); t = 0
    for x in v[i:]:
        if x[0] >= e: break
        t += max(0, min(e, x[1]) - max(s, x[0]))
    return t


fai = {l.split()[0]: list(map(int, l.split()[1:5])) for l in open(FA + '.fai')}; fh_ = open(FA, 'rb')


def gc_of(c, s, e):
    L, off, lb, lB = fai[c]; fh_.seek(off + (s // lb) * lB + s % lb); q = fh_.read(int((e - s) * lB / lb) + 2 * lB).replace(b'\n', b'').replace(b'\r', b'')[:e - s].upper()
    n = len(q) - q.count(b'N'); return 100 * (q.count(b'G') + q.count(b'C')) / max(1, n)


def stats(c, s, e, aid, seg, src, C):
    """one run -> dict (None when no confident copy)"""
    cs = [x for x in C.get(aid, []) if x[0] >= s and x[1] <= e]
    conf = [x for x in cs if x[3] and x[2] >= 80]
    if not conf: return None
    u = per.get(aid, 0); ucls = 'micro' if u <= 6 else 'mini' if u < 100 else 'satellite'
    mid = sum(x[2] for x in conf) / len(conf); rev = sum(x[4] == '-' for x in conf) / len(conf)
    domcls = collections.Counter(x[6] for x in conf).most_common(1)[0][0]
    fam = collections.Counter(x[5] for x in conf).most_common(1)[0][0] if src == 'existing' else famnew.get((aid, domcls), f'new:{aid}')
    rgc = gc_of(c, s, e); cg = cgc.get((aid, domcls))
    return dict(chrom=c, start0=s, end=e, bp=e - s, source=src, array_id=f'{aid}~g{seg}' if NSEG[aid] > 1 else aid, array_id_v46=aid, gap_segment=seg,
                unit_bp=int(u), unit_class=ucls, n_conf=len(conf), n_weak=sum(1 for x in cs if 65 <= x[2] < 80 or (not x[3] and x[2] >= 80)),
                n_junk=sum(1 for x in cs if x[2] < 65), mean_identity_conf=round(mid, 2), reverse_fraction=round(rev, 3), run_GC=round(rgc, 2),
                consensus_GC=cg if cg is not None else 'NA', dominant_class=domcls, family=fam, _gc_diff=abs(rgc - cg) if cg is not None else None, _mid=mid, _rev=rev)


rows = []; dropped = []
for src, R, C in (('existing', R_old, C_old), ('new', R_new, C_new)):
    for c, s, e, aid, seg in R:
        d = stats(c, s, e, aid, seg, src, C)
        if d is None: dropped.append((c, s, e, aid, src, 'no confident copy')); continue
        rows.append(d)
if a.carry_new:                                              # Mac: the new-array 03 table is lost; take the draft rows as they are
    for r in csv.DictReader(open(a.carry_new), delimiter='\t'):
        if r['source'] != 'new': continue
        d = dict(chrom=r['chrom'], start0=int(r['start0']), end=int(r['end']), bp=int(r['bp']), source='new', array_id=r['array_id'], array_id_v46=r.get('array_id_v46', r['array_id']),
                 gap_segment=int(r.get('gap_segment', 1)), unit_bp=int(r['unit_bp']), unit_class='micro' if int(r['unit_bp']) <= 6 else 'mini' if int(r['unit_bp']) < 100 else 'satellite',
                 n_conf=int(r['n_conf']), n_weak=int(r['n_weak']), n_junk=int(r['n_junk']), mean_identity_conf=float(r['mean_identity_conf']), reverse_fraction=float(r['reverse_fraction']),
                 run_GC=float(r['run_GC']), consensus_GC=r['consensus_GC'], dominant_class=r.get('dominant_class', 'NA'), family=r.get('family_v46', r['family']),
                 _gc_diff=abs(float(r['run_GC']) - float(r['consensus_GC'])) if r['consensus_GC'] != 'NA' else None, _mid=float(r['mean_identity_conf']), _rev=float(r['reverse_fraction']))
        d['carried_from'] = os.path.basename(a.carry_new); rows.append(d)
    print(f'carried {sum(1 for r in rows if r.get("carried_from"))} new runs from {a.carry_new}')
# ---- overlap resolution (C08)
rows.sort(key=lambda r: (r['chrom'], r['start0'], r['end'])); res = []
for r in rows:
    r.setdefault('overlap_note', 'none')
    if res and res[-1]['chrom'] == r['chrom'] and r['start0'] < res[-1]['end']:
        p = res[-1]
        if p['family'] == r['family'] and p['source'] != r['source'] and not (p.get('carried_from') or r.get('carried_from')):
            C = C_old if p['source'] == 'existing' else C_new; aid = p['array_id_v46']
            m = stats(p['chrom'], min(p['start0'], r['start0']), max(p['end'], r['end']), aid, p['gap_segment'], p['source'], C) or dict(p)
            m['overlap_note'] = f"merged_with:{r['source']}:{r['array_id']}"; res[-1] = m; continue
        loser, keeper = (r, p) if r['n_conf'] <= p['n_conf'] else (p, r)
        if loser is r:
            if r['end'] <= p['end']: r['overlap_note'] = f"contained_in:{p['array_id']}:dropped"; dropped.append((r['chrom'], r['start0'], r['end'], r['array_id_v46'], r['source'], 'contained in ' + p['array_id'])); continue
            C = C_old if r['source'] == 'existing' else C_new
            t = stats(r['chrom'], p['end'], r['end'], r['array_id_v46'], r['gap_segment'], r['source'], C) if not r.get('carried_from') else None
            if t is None:
                r['overlap_note'] = f"overlaps:{p['array_id']}:dropped"; dropped.append((r['chrom'], r['start0'], r['end'], r['array_id_v46'], r['source'], 'trimmed to no copy by ' + p['array_id'])); continue
            t['overlap_note'] = f"trimmed_start_by:{p['array_id']}"; r = t
        else:
            C = C_old if p['source'] == 'existing' else C_new
            t = stats(p['chrom'], p['start0'], r['start0'], p['array_id_v46'], p['gap_segment'], p['source'], C) if not p.get('carried_from') else None
            if t is None:
                p['overlap_note'] = f"overlaps:{r['array_id']}:dropped"; dropped.append((p['chrom'], p['start0'], p['end'], p['array_id_v46'], p['source'], 'trimmed to no copy by ' + r['array_id'])); res.pop()
            else: t['overlap_note'] = f"trimmed_end_by:{r['array_id']}"; res[-1] = t
    res.append(r)
rows = res
# ---- run ids, gates
k = collections.Counter()
for r in rows:
    k[r['source']] += 1; r['run_id'] = f'{G}_{"E" if r["source"] == "existing" else "N"}{k[r["source"]]:05d}'; r['subgenome'] = G
    r['G1'] = 'PASS' if r['_gc_diff'] is not None and r['_gc_diff'] <= 5 else 'FAIL'
    r['G2'] = 'PASS' if r['_mid'] >= 85 else 'FAIL'
    r['G3'] = 'PASS' if not (0.35 < r['_rev'] < 0.65) else 'FAIL'
    r['gates_passed'] = (r['G1'], r['G2'], r['G3']).count('PASS')
with open(os.path.join(DEST, f'59_gates_sweep_{G}.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n'); w.writerow(['G1_gc_points', 'G2_identity', 'G3_band', 'n_runs', 'n_pass_all', 'Mb_pass_all', 'fail_G1', 'fail_G2', 'fail_G3'])
    for g1 in (3, 5, 8):
        for g2 in (80, 85, 90):
            for lo, hi in ((0.3, 0.7), (0.35, 0.65), (0.4, 0.6)):
                f1 = sum(1 for r in rows if r['_gc_diff'] is None or r['_gc_diff'] > g1); f2 = sum(1 for r in rows if r['_mid'] < g2); f3 = sum(1 for r in rows if lo < r['_rev'] < hi)
                ok = [r for r in rows if r['_gc_diff'] is not None and r['_gc_diff'] <= g1 and r['_mid'] >= g2 and not (lo < r['_rev'] < hi)]
                w.writerow([g1, g2, f'{lo}-{hi}', len(rows), len(ok), round(sum(r['bp'] for r in ok) / 1e6, 3), f1, f2, f3])
# ---- gene arrays (59b as a module), element_class before tiers
spec = importlib.util.spec_from_file_location('ga', os.path.join(os.path.dirname(os.path.abspath(__file__)), '59b_gene_arrays.py')); ga = importlib.util.module_from_spec(spec); spec.loader.exec_module(ga)
calls, sweep = ga.classify(rows, G, OVL, cut=a.gene_cut, unit_cut=a.unit_cut)
ga.write(calls, sweep, os.path.join(DEST, f'59b_gene_arrays_{G}.tsv'))
for r in rows:
    cl = calls.get(r['family']); r['element_class'] = cl['element_class'] if cl else 'satellite'; r['gene_flag'] = cl['gene_flag'] if cl else 'none'
    r['seq_cluster'] = clus.get(r['family'], 'none') if r['element_class'] == 'satellite' else 'na'
    r['EDTA_overlap_fraction'] = round(min(1.0, ov(te, r['chrom'], r['start0'], r['end']) / r['bp']), 3)
    r['putative_centromeric_array_L1'] = ('yes' if ov(cen, r['chrom'], r['start0'], r['end']) > 0 else 'no') if r['element_class'] == 'satellite' else 'na'
    r['evidence'] = 'L1'
# candidate centromeric array from the runs of record: largest cluster (<= 100 kb apart) of satellite-class runs per chromosome
byc = collections.defaultdict(list)
for i, r in enumerate(rows):
    if r['element_class'] == 'satellite' and r['unit_class'] == 'satellite': byc[r['chrom']].append((r['start0'], r['end'], i))
for c, v in byc.items():
    v.sort(); groups = []; cur = []
    for s, e, i in v:
        if cur and s - cur[-1][1] > 100_000: groups.append(cur); cur = []
        cur.append((s, e, i))
    if cur: groups.append(cur)
    best = max(groups, key=lambda g: sum(e - s for s, e, _ in g)); bi = {i for _, _, i in best}
    for s, e, i in v: rows[i]['cen_candidate_ror'] = 'yes' if i in bi else 'no'
for r in rows: r.setdefault('cen_candidate_ror', 'na')
# ---- tiers on satellite families only (C04)
fbp = collections.Counter()
for r in rows:
    if r['element_class'] == 'satellite': fbp[r['family']] += r['bp']
tot = sum(fbp.values()); cum = 0; tier = {}
for f, b in fbp.most_common():
    cum += b; frac_before = (cum - b) / tot
    tier[f] = 'core' if frac_before < 0.5 else 'major' if frac_before < 0.9 else 'minor' if frac_before < 0.99 else 'rare'
for r in rows: r['tier'] = tier.get(r['family'], 'gene_array')
cols = ['run_id', 'subgenome', 'chrom', 'start0', 'end', 'bp', 'source', 'array_id', 'array_id_v46', 'gap_segment', 'unit_bp', 'unit_class', 'n_conf', 'n_weak', 'n_junk',
        'mean_identity_conf', 'reverse_fraction', 'run_GC', 'consensus_GC', 'G1', 'G2', 'G3', 'gates_passed', 'dominant_class', 'family', 'element_class', 'gene_flag', 'seq_cluster',
        'EDTA_overlap_fraction', 'putative_centromeric_array_L1', 'cen_candidate_ror', 'tier', 'overlap_note', 'carried_from', 'evidence']
for r in rows: r.setdefault('carried_from', 'none')
with open(os.path.join(DEST, f'59_catalogue_{G}.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', lineterminator='\n', extrasaction='ignore'); w.writeheader(); w.writerows(rows)
PAL = ['31,119,180', '214,39,40', '148,103,189', '44,160,44', '255,127,14', '140,86,75', '227,119,194', '23,190,207']
cl_names = sorted({r['seq_cluster'] for r in rows if r['seq_cluster'] not in ('none', 'na')}); COL = {c: PAL[i % len(PAL)] for i, c in enumerate(cl_names)}
with open(os.path.join(DEST, f'59_catalogue_{G}.bed'), 'w') as fh:
    for r in rows:
        fh.write(f"{r['chrom']}\t{r['start0']}\t{r['end']}\t{r['run_id']}|{r['family']}|{r['element_class']}|{r['seq_cluster']}|{r['unit_bp']}bp|{r['tier']}\t0\t.\t{r['start0']}\t{r['end']}\t{COL.get(r['seq_cluster'], '150,150,150')}\n")
with open(os.path.join(DEST, f'59_family_tiers_{G}.tsv'), 'w') as fh:
    fh.write('family\tbp\tshare\ttier\tseq_cluster\telement_class\n')
    for f, b in fbp.most_common(): fh.write(f'{f}\t{b}\t{b / tot:.4f}\t{tier[f]}\t{clus.get(f, "none")}\tsatellite\n')
    gb = collections.Counter(); ge = {}
    for r in rows:
        if r['element_class'] != 'satellite': gb[r['family']] += r['bp']; ge[r['family']] = r['element_class']
    for f, b in gb.most_common(): fh.write(f'{f}\t{b}\tNA\t{ge[f].split(":")[0]}\tna\t{ge[f]}\n')
with open(os.path.join(DEST, f'59_dropped_runs_{G}.tsv'), 'w') as fh:
    fh.write('chrom\tstart0\tend\tarray_id\tsource\treason\n'); [fh.write('\t'.join(map(str, d)) + '\n') for d in dropped]
S = collections.Counter(); B = collections.Counter()
for r in rows:
    for kk in (('source', r['source']), ('unit_class', r['unit_class']), ('element_class', r['element_class'].split(':')[0]), ('seq_cluster', r['seq_cluster']), ('tier', r['tier']),
               ('gates_passed', r['gates_passed']), ('putative_centromeric_array_L1', r['putative_centromeric_array_L1']), ('cen_candidate_ror', r['cen_candidate_ror']),
               ('overlap_note', r['overlap_note'].split(':')[0])):
        S[kk] += 1; B[kk] += r['bp']
with open(os.path.join(DEST, f'59_summary_{G}.tsv'), 'w') as fh:
    fh.write('kind\tvalue\tn_runs\tMb\n')
    for kk in sorted(S, key=str): fh.write(f'{kk[0]}\t{kk[1]}\t{S[kk]}\t{B[kk] / 1e6:.3f}\n')
print(G, len(rows), 'runs,', round(sum(r['bp'] for r in rows) / 1e6, 3), 'Mb (satellite families', round(tot / 1e6, 3), 'Mb);', len(dropped), 'runs dropped'); print(open(os.path.join(DEST, f'59_summary_{G}.tsv')).read())
