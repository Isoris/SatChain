#!/usr/bin/env python3
"""Family of each new array class: BLAST of the new class consensuses (new03 03_array_consensus.tsv) against all family consensuses of both
subgenomes; best hit with identity >= 80 % over >= 50 % of the consensus; a same-subgenome family wins over a cross-subgenome one at
equal score; otherwise 'new:<array>'. usage: famnew.py export <G> <new03 dir> <fa> | assign <G> <new03 dir> <blast tsv> <out tsv>"""
import sys, csv, collections
csv.field_size_limit(1 << 30)
mode, G, D = sys.argv[1:4]
if mode == 'export':
    with open(sys.argv[4], 'w') as fh:
        for r in csv.DictReader(open(f'{D}/03_array_consensus.tsv'), delimiter='\t'): fh.write(f">{r['array_id']}|{r['monomer_class']}\n{r['consensus_sequence']}\n")
else:
    best = {}
    for l in open(sys.argv[4]):
        q, s, pid, ln, *_ , ql, sl = l.rstrip('\n').split('\t'); pid = float(pid); cov = int(ln) / int(ql)
        if pid < 80 or cov < 0.5: continue
        g, fam = s.split('|')[0], s.split('|')[1]; sc = pid * cov + (1 if g == G else 0)
        if q not in best or sc > best[q][0]: best[q] = (sc, fam, g, pid, cov)
    with open(sys.argv[5], 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n'); w.writerow(['array_id', 'monomer_class', 'family', 'hit_subgenome', 'identity', 'coverage'])
        n = collections.Counter()
        for r in csv.DictReader(open(f'{D}/03_array_consensus.tsv'), delimiter='\t'):
            q = f"{r['array_id']}|{r['monomer_class']}"; b = best.get(q)
            if b: w.writerow([r['array_id'], r['monomer_class'], b[1], b[2], b[3], round(b[4], 3)]); n['same' if b[2] == G else 'cross'] += 1
            else: w.writerow([r['array_id'], r['monomer_class'], f"new:{r['array_id']}", 'NA', 'NA', 'NA']); n['new'] += 1
    print(G, dict(n))
