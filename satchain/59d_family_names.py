#!/usr/bin/env python3
"""59d: families for the new arrays and one satellitome name for every catalogue family.
New-array class consensuses with no existing-family match (famnew 'new:') are grouped by single linkage on the all-vs-all BLAST
(identity >= 80 % over >= 50 % of the shorter consensus) -> 'NewSat<G><n>'. The family of a new run = that of its DOMINANT CLASS AS
COMPUTED BY 59 (column dominant_class, per run; 2026-11-01, audit C14: 59d used to recompute an array-level dominant class and disagreed
with 59 for 17 / 1 runs). Copies are looked up by array_id_v46 (audit C10).
All satellite families (gene arrays excluded) are then ranked by catalogue bp: name = Cga|Cma + 'Sat' + rank (3 digits) + '-' + median unit.
Tiers are RECOMPUTED here on family_v46 and satellite bp (2026-11-01, audit C15: before, a NewSat family carried the tier of its largest
single array). Gene arrays keep 'gene_array:<class>' with the family id. Runs without a 03 table (carried_from != none) keep family_v46 = family.
Outputs: 59_catalogue_<G>.tsv (family_v46, name_v46, tier updated), 59d_names_<G>.tsv
usage: 59d_family_names.py <G> <cat59 dir> <new03 dir|->"""
import sys, csv, collections, os
csv.field_size_limit(1 << 30)
G, C, N = sys.argv[1:4]; ab = 'Cga' if G == 'Gar' else 'Cma'; N = '' if N == '-' else N
runs = list(csv.DictReader(open(f'{C}/59_catalogue_{G}.tsv'), delimiter='\t'))
fam = {}; newname = {}; comp = {}
if N and os.path.exists(f'{N}/famnew.tsv'):
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
    for i, (root, ks) in enumerate(sorted(comp.items(), key=lambda kv: -len(kv[1]))):
        for k in ks: newname[k] = f'NewSat{G}{i + 1:05d}'
for r in runs:
    aid = r.get('array_id_v46') or r['array_id'].split('~g')[0]
    if r['source'] == 'new' and r.get('carried_from', 'none') == 'none' and fam:
        d = r.get('dominant_class', 'C01'); f = fam.get((aid, d), r['family']); r['family_v46'] = newname.get((aid, d), f)
    else: r['family_v46'] = r['family']
bp = collections.Counter(); units = collections.defaultdict(list)
for r in runs:
    if r['element_class'] == 'satellite': bp[r['family_v46']] += int(r['bp']); units[r['family_v46']].append(int(r['unit_bp']))
tot = sum(bp.values()); cum = 0; tier = {}
for f, b in bp.most_common():
    cum += b; fb = (cum - b) / tot
    tier[f] = 'core' if fb < 0.5 else 'major' if fb < 0.9 else 'minor' if fb < 0.99 else 'rare'
name = {}
for i, (f, b) in enumerate(bp.most_common()):
    u = sorted(units[f])[len(units[f]) // 2]; name[f] = f'{ab}Sat{i + 1:03d}-{u}' if i < 999 else f'{ab}Sat{i + 1}-{u}'
for r in runs:
    r['name_v46'] = name.get(r['family_v46'], f"{r['element_class']}|{r['family_v46']}")
    r['tier'] = tier.get(r['family_v46'], 'gene_array')
with open(f'{C}/59_catalogue_{G}.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(runs[0].keys()), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(runs)
cl = {r['family_v46']: r.get('seq_cluster', 'none') for r in runs}
with open(f'{C}/59d_names_{G}.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n'); w.writerow(['name_v46', 'family_v46', 'bp', 'share', 'unit_bp_median', 'tier', 'seq_cluster', 'n_runs'])
    nr = collections.Counter(r['family_v46'] for r in runs if r['element_class'] == 'satellite')
    for f, b in bp.most_common(): w.writerow([name[f], f, b, round(b / tot, 4), sorted(units[f])[len(units[f]) // 2], tier[f], cl.get(f, 'none'), nr[f]])
with open(f'{C}/59_family_tiers_{G}.tsv', 'w') as fh:                       # tiers table on family_v46
    fh.write('family\tbp\tshare\ttier\tseq_cluster\telement_class\n')
    for f, b in bp.most_common(): fh.write(f'{f}\t{b}\t{b / tot:.4f}\t{tier[f]}\t{cl.get(f, "none")}\tsatellite\n')
    gb = collections.Counter(); ge = {}
    for r in runs:
        if r['element_class'] != 'satellite': gb[r['family_v46']] += int(r['bp']); ge[r['family_v46']] = r['element_class']
    for f, b in gb.most_common(): fh.write(f'{f}\t{b}\tNA\t{ge[f].split(":")[0]}\tna\t{ge[f]}\n')
print(G, 'new families', len(comp), '| satellite families named', len(name), '| tiers', dict(collections.Counter(tier.values())), '| top:', [(name[f], f) for f, _ in bp.most_common(6)])
