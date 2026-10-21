#!/usr/bin/env python3
"""55a: sample confident monomer copies (identity >= 80 %, complete) from every step-53 segment of the centromeric arrays of both
subgenomes, up to --per-seg copies per segment (seeded). Writes 55_sample.fa and 55_sample.tsv (copy metadata incl. the OLD family
lineage label, kept only for comparison). usage: 55a_lineage_sample.py <dest> <per_seg> <G>=<pipeline out>|<53 run dir> ..."""
import csv, collections, random, sys, os
csv.field_size_limit(1 << 30)
DEST, PER = sys.argv[1], int(sys.argv[2]); os.makedirs(DEST, exist_ok=True)
OLD = {}
# optional: a previous family-level lineage table (assembly, lineage, families), kept only for comparison: set SATCHAIN_OLD_LINEAGE
if os.environ.get('SATCHAIN_OLD_LINEAGE'):
    for r in csv.DictReader(open(os.environ['SATCHAIN_OLD_LINEAGE']), delimiter='\t'):
        for f in r['families'].split(','): OLD[(r['assembly'], f)] = r['lineage']
rng = random.Random(55)
fa = open(os.path.join(DEST, '55_sample.fa'), 'w'); tw = csv.writer(open(os.path.join(DEST, '55_sample.tsv'), 'w', newline=''), delimiter='\t', lineterminator='\n')
tw.writerow(['copy_id', 'subgenome', 'array_id', 'segment', 'segment_profile', 'chrom', 'start0', 'end', 'length', 'family', 'class', 'strand', 'identity', 'old_lineage'])
for spec in sys.argv[3:]:
    G, rest = spec.split('=', 1); OUT, RUN = rest.split('|', 1)
    segs = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(RUN, '53_segments.tsv')), delimiter='\t'): segs[r['array_id']].append(r)
    pool = collections.defaultdict(list)
    with open(os.path.join(OUT, '03_monomer_copies.tsv'), newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            aid = r['array_id']
            if aid not in segs or r['complete_or_truncated'] != 'complete': continue
            try: idt = float(r['best_alignment_identity'])
            except ValueError: continue
            if idt < 80: continue
            s = int(r['genomic_start_0based'])
            for sg in segs[aid]:
                if int(sg['start0']) <= s < int(sg['end']): pool[(aid, sg['segment'], sg['profile'])].append((r, idt)); break
    n = 0
    for (aid, sg, prof), cs in pool.items():
        for r, idt in rng.sample(cs, min(PER, len(cs))):
            cid = f"{G}|{r['monomer_id']}"; seq = r['sequence'].upper()
            fa.write(f'>{cid}\n{seq}\n'); n += 1
            tw.writerow([cid, G, aid, sg, prof, r['chrom'], r['genomic_start_0based'], r['genomic_end'], len(seq), r['family_id'], r['monomer_class'],
                         r['strand'], idt, OLD.get((G, r['family_id']), 'none')])
    print(G, n, 'copies from', len(pool), 'segments')
