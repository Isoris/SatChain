#!/usr/bin/env python3
"""59d: families for the new arrays and one satellitome name for every catalogue family.
New-array class consensuses with no existing-family match (famnew 'new:') are grouped by single linkage on the all-vs-all BLAST
(identity >= 80 % over >= 50 % of the shorter consensus) -> 'NewSat<G><n>'. The family of a new run = that of its dominant class.
All satellite families (gene arrays excluded) are then ranked by catalogue bp: name = Cga|Cma + 'Sat' + rank (3 digits) + '-' + median unit.
Gene arrays keep 'gene_array:<class>' with the family id. Outputs: 59_catalogue_<G>.tsv (family_v46, name_v46 added), 59d_names_<G>.tsv
usage: 59d_family_names.py <G> <cat59 dir> <new03 dir>"""
import sys, csv, collections
csv.field_size_limit(1 << 30)
G, C, N = sys.argv[1:4]; ab = 'Cga' if G == 'Gar' else 'Cma'
fam = {(r['array_id'], r['monomer_class']): r['family'] for r in csv.DictReader(open(f'{N}/famnew.tsv'), delimiter='\t')}
parent = {}
def find(x):
    parent.setdefault(x, x)
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
for l in open(f'{N}/new_allvsall.tsv'):
    q, s, pid, ln, *_ , ql, sl = l.rstrip('\n').split('\t')
    if q == s: continue
    if float(pid) >= 80 and int(ln) >= 0.5 * min(int(ql), int(sl)):
        a, b = tuple(q.split('|')), tuple(s.split('|'))
        if fam.get(a, '').startswith('new:') and fam.get(b, '').startswith('new:'): parent[find(b)] = find(a)
comp = collections.defaultdict(list)
for k, f in fam.items():
    if f.startswith('new:'): comp[find(k)].append(k)
newname = {}
for i, (root, ks) in enumerate(sorted(comp.items(), key=lambda kv: -len(kv[1]))):
    for k in ks: newname[k] = f'NewSat{G}{i + 1:05d}'
# dominant class of each new run
cls = collections.defaultdict(collections.Counter)
runs = list(csv.DictReader(open(f'{C}/59_catalogue_{G}.tsv'), delimiter='\t'))
need = {r['array_id'] for r in runs if r['source'] == 'new'}
with open(f'{N}/03_monomer_copies.tsv', newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if r['array_id'] in need and r['complete_or_truncated'] == 'complete':
            try:
                if float(r['best_alignment_identity']) >= 80: cls[r['array_id']][r['monomer_class']] += 1
            except ValueError: pass
for r in runs:
    if r['source'] == 'new':
        d = cls[r['array_id']].most_common(1)[0][0] if cls[r['array_id']] else 'C01'
        f = fam.get((r['array_id'], d), r['family']); r['family_v46'] = newname.get((r['array_id'], d), f)
    else: r['family_v46'] = r['family']
bp = collections.Counter(); units = collections.defaultdict(list)
for r in runs:
    if r['element_class'] == 'satellite': bp[r['family_v46']] += int(r['bp']); units[r['family_v46']].append(int(r['unit_bp']))
name = {}
for i, (f, b) in enumerate(bp.most_common()):
    u = sorted(units[f])[len(units[f]) // 2]; name[f] = f'{ab}Sat{i + 1:03d}-{u}' if i < 999 else f'{ab}Sat{i + 1}-{u}'
for r in runs: r['name_v46'] = name.get(r['family_v46'], f"{r['element_class']}|{r['family_v46']}")
with open(f'{C}/59_catalogue_{G}.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(runs[0].keys()), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(runs)
lin = {r['family_v46']: r['lineage'] for r in runs}
with open(f'{C}/59d_names_{G}.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n'); w.writerow(['name_v46', 'family_v46', 'bp', 'share', 'unit_bp_median', 'lineage', 'n_runs'])
    nr = collections.Counter(r['family_v46'] for r in runs if r['element_class'] == 'satellite'); tot = sum(bp.values())
    for f, b in bp.most_common(): w.writerow([name[f], f, b, round(b / tot, 4), sorted(units[f])[len(units[f]) // 2], lin.get(f, 'none'), nr[f]])
print(G, 'new families', len(comp), '| satellite families named', len(name), '| top:', [(name[f], f) for f, _ in bp.most_common(6)])
