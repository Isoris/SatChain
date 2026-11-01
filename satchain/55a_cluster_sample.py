#!/usr/bin/env python3
"""55a (2026-11-01; replaces 55a_lineage_sample.py): sample confident monomer copies (identity >= 80 %, complete) from the catalogue runs of
record of one or more haplotypes, for the k-mer composition clustering of 55c. Pooled design: give every haplotype in ONE call (the
published clustering pooled Gar and Mac; per-haplotype partitions are computed by 55c on the same sample, audit C18).
Only satellite runs enter (element_class == 'satellite'; gene arrays and gene-related families are excluded, audit C04); by default the
candidate centromeric runs (putative_centromeric_array_L1 == yes or cen_candidate_ror == yes), --all-satellite takes every satellite run.
Up to --per-run copies per run, drawn with a generator seeded on the run id (adding or removing runs does not reshuffle the others).
No step-53 segments are needed any more (the old 55a sampled per 53 segment, which made 59 -> 53 -> 55 -> 59 circular, audit C09).
Outputs: <dest>/55_sample.fa, 55_sample.tsv (copy_id, subgenome, run_id, array_id, family, chrom, start0, end, length, class, strand, identity, cen flags)
usage: 55a_cluster_sample.py <dest> <per_run> [--all-satellite] <G>=<pipeline out>|<catalogue tsv>[|<new03 dir>] ..."""
import csv, collections, random, sys, os, zlib
csv.field_size_limit(1 << 30)
args = [x for x in sys.argv[1:] if x != '--all-satellite']; ALL = '--all-satellite' in sys.argv
DEST, PER = args[0], int(args[1]); os.makedirs(DEST, exist_ok=True)
fa = open(os.path.join(DEST, '55_sample.fa'), 'w'); tw = csv.writer(open(os.path.join(DEST, '55_sample.tsv'), 'w', newline=''), delimiter='\t', lineterminator='\n')
tw.writerow(['copy_id', 'subgenome', 'run_id', 'array_id', 'family', 'chrom', 'start0', 'end', 'length', 'monomer_class', 'strand', 'identity', 'cen_L1', 'cen_ror'])
for spec in args[2:]:
    G, rest = spec.split('=', 1); parts = rest.split('|'); OUT, CAT = parts[0], parts[1]; N03 = parts[2] if len(parts) > 2 else ''
    runs = collections.defaultdict(list); nrun = 0
    for r in csv.DictReader(open(CAT), delimiter='\t'):
        if r['subgenome'] != G or r['element_class'] != 'satellite' or r.get('carried_from', 'none') != 'none': continue
        if not ALL and not (r.get('putative_centromeric_array_L1') == 'yes' or r.get('cen_candidate_ror') == 'yes'): continue
        runs[r['array_id_v46']].append(r); nrun += 1
    pool = collections.defaultdict(list)
    tables = [os.path.join(OUT, f'fClaHyb_{G}_satellites', '03_monomer_copies.tsv')] + ([os.path.join(N03, '03_monomer_copies.tsv')] if N03 else [])
    for path in tables:
        if not os.path.exists(path): continue
        with open(path, newline='') as fh:
            prev = None
            for r in csv.DictReader(fh, delimiter='\t'):
                aid = r['array_id']
                if aid not in runs or r['complete_or_truncated'] != 'complete': continue
                try: idt = float(r['best_alignment_identity'])
                except ValueError: continue
                if idt < 80: continue
                s, e = int(r['genomic_start_0based']), int(r['genomic_end'])
                if prev == (aid, s, e): continue                                  # CB01 duplicate tile
                prev = (aid, s, e)
                for rr in runs[aid]:
                    if int(rr['start0']) <= s and e <= int(rr['end']): pool[rr['run_id']].append((r, idt, rr)); break
    n = 0
    for rid, cs in sorted(pool.items()):
        rng = random.Random(zlib.crc32(rid.encode()))
        for r, idt, rr in rng.sample(cs, min(PER, len(cs))):
            cid = f"{G}|{r['monomer_id']}"; seq = r['sequence'].upper()
            fa.write(f'>{cid}\n{seq}\n'); n += 1
            tw.writerow([cid, G, rid, rr['array_id'], rr['family'], r['chrom'], r['genomic_start_0based'], r['genomic_end'], len(seq), r['monomer_class'], r['strand'], idt,
                         rr.get('putative_centromeric_array_L1', 'na'), rr.get('cen_candidate_ror', 'na')])
    print(G, n, 'copies from', len(pool), 'of', nrun, 'eligible runs', '(all satellite runs)' if ALL else '(candidate centromeric runs)')
