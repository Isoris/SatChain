#!/usr/bin/env python3
"""
59c_te_layer.py  (v4.6; TE content of the catalogue runs without counting EDTA calls on the satellite itself)

Per run, every EDTA feature overlapping it (repeat_region / TSD / LTR-terminal records excluded), grouped by EDTA label:
  carried share = fraction of the run's confident copies that one feature of the label covers >= 50 % of (min of copy / feature length)
  intrinsic     carried share > --intrinsic (swept 0.2 / 0.3 / 0.5) -> the label describes the satellite (EDTA calling it TE): kind
                'uniform' (>= 0.95) or 'partial'; its bp are NOT TE insertions
  insertion     all other features of the run: TE inserted in the array (bp, n, superfamilies)
Plus the step-54 insertions of known-TE families (EDTA missed them inside arrays): members of 54c 'known_TE' families inside the run.
Outputs: 59c_te_<G>.tsv (per run), summary printed (per intrinsic cut). usage: 59c_te_layer.py <G> <catalogue> <pipeline out> <new03> <func> <ins54 dir>
"""
import sys, csv, gzip, collections, bisect, os
csv.field_size_limit(1 << 30)
G, CAT, OUT, N03, FUNC, INS = sys.argv[1:7]
rows = list(csv.DictReader(open(CAT), delimiter='\t'))
need = collections.defaultdict(set)
for r in rows: need['existing' if r['source'] == 'existing' else 'new'].add(r['array_id'])
cop = collections.defaultdict(list)
for src, path in (('existing', os.path.join(OUT, f'fClaHyb_{G}_satellites', '03_monomer_copies.tsv')), ('new', os.path.join(N03, '03_monomer_copies.tsv'))):
    with open(path, newline='') as fh:
        rd = csv.reader(fh, delimiter='\t'); h = next(rd); ix = {k: i for i, k in enumerate(h)}
        for row in rd:
            if row[ix['array_id']] not in need[src]: continue
            try: idt = float(row[ix['best_alignment_identity']])
            except ValueError: continue
            if row[ix['complete_or_truncated']] == 'complete' and idt >= 80: cop[row[ix['array_id']]].append((int(row[ix['genomic_start_0based']]), int(row[ix['genomic_end']])))
for v in cop.values(): v.sort()
te = collections.defaultdict(list)
for l in gzip.open(os.path.join(FUNC, '_repeats_EDTA', f'fClaHyb_{G}.TEanno.sorted.gff3.gz'), 'rt'):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] in ('repeat_region', 'target_site_duplication', 'long_terminal_repeat'): continue
    cl = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x).get('classification', f[2])
    te[f[0]].append((int(f[3]) - 1, int(f[4]), cl))
for v in te.values(): v.sort()
tst = {c: [x[0] for x in v] for c, v in te.items()}; tml = {c: max(x[1] - x[0] for x in v) for c, v in te.items()}
known = set()
cls_path = os.path.join(INS, '54c_ins_families_classified.tsv')
for r in csv.DictReader(open(cls_path), delimiter='\t'):
    if r['class'] == 'known_TE': known.add(r['family'])
kte = collections.defaultdict(list)
for r in csv.DictReader(open(os.path.join(INS, '54b_ins_members.tsv')), delimiter='\t'):
    if r['family'] in known and r['ins_id'].startswith(G + '|'):
        loc = r['ins_id'].split('|')[2]; c, se = loc.split(':'); s, e = map(int, se.split('-')); kte[c].append((s, e, r['family']))
for v in kte.values(): v.sort()
CUTS = (0.2, 0.3, 0.5)
out = []; summ = {c: collections.Counter() for c in CUTS}
for r in rows:
    c, s, e = r['chrom'], int(r['start0']), int(r['end']); cs = [x for x in cop.get(r['array_id'], []) if x[0] >= s and x[1] <= e]
    i0 = bisect.bisect_left(tst.get(c, []), s - tml.get(c, 0)); i1 = bisect.bisect_left(tst.get(c, []), e)
    F = [x for x in te.get(c, [])[i0:i1] if x[1] > s]
    by = collections.defaultdict(list)
    for x in F: by[x[2]].append(x)
    carried = {}
    for lab, v in by.items():
        vs = sorted(v); vst = [x[0] for x in vs]; ml = max(x[1] - x[0] for x in vs); hit = 0
        for a, b in cs:
            j0, j1 = bisect.bisect_left(vst, a - ml), bisect.bisect_left(vst, b)
            hit += any(min(b, vs[k][1]) - max(a, vs[k][0]) >= 0.5 * min(b - a, vs[k][1] - vs[k][0]) for k in range(j0, j1))
        carried[lab] = hit / max(1, len(cs))
    def clipped(xs):
        iv = sorted((max(s, x[0]), min(e, x[1])) for x in xs); t = 0; cur = None
        for a, b in iv:
            if cur and a <= cur[1]: cur[1] = max(cur[1], b)
            else:
                if cur: t += cur[1] - cur[0]
                cur = [a, b]
        return t + (cur[1] - cur[0] if cur else 0)
    rec = dict(run_id=r['run_id'], chrom=c, start0=s, end=e, bp=e - s, family=r['family'], lineage=r['lineage'], element_class=r['element_class'],
               EDTA_overlap_bp=clipped(F), labels_carried=';'.join(f'{k}:{v:.2f}' for k, v in sorted(carried.items(), key=lambda kv: -kv[1])[:4]))
    for cut in CUTS:
        intr = {k for k, v in carried.items() if v > cut}
        ins = [x for x in F if x[2] not in intr]
        rec[f'intrinsic_labels_{cut}'] = ','.join(sorted(intr)) or 'none'
        rec[f'TE_insert_bp_{cut}'] = clipped(ins); rec[f'TE_insert_n_{cut}'] = len(ins)
        summ[cut]['runs'] += 1; summ[cut]['bp'] += e - s; summ[cut]['edta_bp'] += rec['EDTA_overlap_bp']; summ[cut]['ins_bp'] += rec[f'TE_insert_bp_{cut}']
        summ[cut]['runs_with_ins'] += rec[f'TE_insert_n_{cut}'] > 0; summ[cut]['runs_intrinsic'] += bool(intr)
    v = kte.get(c, []); j = bisect.bisect_left(v, (s, 0, '')); k = [x for x in v[j:] if x[0] < e]
    rec['missed_TE_insertions_54'] = len(k); rec['missed_TE_bp_54'] = sum(min(e, x[1]) - x[0] for x in k)
    out.append(rec)
with open(f'59c_te_{G}.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(out)
m54 = sum(x['missed_TE_bp_54'] for x in out); n54 = sum(x['missed_TE_insertions_54'] > 0 for x in out)
for cut in CUTS:
    S = summ[cut]
    print(f"{G} intrinsic>{cut}: runs {S['runs']}, satellite {S['bp'] / 1e6:.2f} Mb, EDTA overlap {S['edta_bp'] / 1e6:.2f} Mb -> TE insertions {S['ins_bp'] / 1e6:.2f} Mb "
          f"({S['ins_bp'] / S['bp']:.3f} of run bp) in {S['runs_with_ins']} runs; runs where EDTA labels the satellite: {S['runs_intrinsic']}")
print(f'{G}: + missed known-TE insertions (54) in {n54} runs, {m54 / 1e6:.2f} Mb')
