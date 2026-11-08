#!/usr/bin/env python3
"""igv_tracks_v47.py  (SatChain v4.7, 2026-11-07): IGV tracks of the final catalogue, plain BED / GFF3 / bedGraph (no bgzip needed at
this size), one folder per haplotype with an IGV session file (relative paths) that loads everything in one click.

Purpose (declared): look at the annotation. Nothing here is a new measurement; every column comes from the 59 catalogue, 59g, 59h,
the chromosome loci table, centromere_check.tsv and the v4.5 telomere blocks.

Per haplotype <dest>/<G>/:
  satDNA_v47.runs.gff3                     one feature per run of record, all attributes (family, name, unit, class, 59g decision, copies, identity, tier, cluster, centromere flags)
  satDNA_v47.runs_by_family.bed            BED9: colour = family (12 largest satellite families distinct, others grey, gene arrays black); name = family name
  satDNA_v47.runs_by_decision.bed          BED9: colour = 59g decision (satellite / gene_array / gene_derived_satellite / unresolved)
  satDNA_v47.gene_arrays.bed               BED6 subset: runs of gene-array families, name = family:gene class
  satDNA_v47.largest_satellite_cluster.bed BED6: runs in the largest satellite cluster of their chromosome (cen_candidate_ror)
  satDNA_v47.loci_top3.bed                 BED9: the three largest run clusters per chromosome, colour = locus class, name = rank:class:kb
  candidate_centromere_intervals.bed       BED6: the June synteny-reconciled candidate BED (identical to the step-31 interval), name = support / verdict / dominant family
  telomere_blocks.bed                      BED6: v4.5 step-11 telomere blocks >= 500 bp
  satDNA_v47.density_100kb.bedgraph        satellite bp (satellite-class runs) per 100-kb window
  family_metrics.tsv, unit_content.tsv     copies of the 59h and 59g tables for lookup
  igv_session.xml, README.md
usage: igv_tracks_v47.py <final dir> <annotation root> <v4.5 out root> <loci tsv> <dest> [--tags Gar,Mac]
"""
import sys, os, csv, argparse, collections, shutil, colorsys
ap = argparse.ArgumentParser(); ap.add_argument('final'); ap.add_argument('ann'); ap.add_argument('out45'); ap.add_argument('loci'); ap.add_argument('dest'); ap.add_argument('--tags', default='Gar,Mac')
a = ap.parse_args(); tags = a.tags.split(',')
DEC = {'satellite': (138, 28, 43), 'gene_array': (17, 17, 17), 'gene_derived_satellite': (217, 119, 6), 'unresolved': (111, 66, 193)}
LOC = {'satellite array': (138, 28, 43), 'tandem gene array': (17, 17, 17), 'gene-related satellite (tRNA-derived)': (154, 107, 0), 'minisatellite cluster': (74, 93, 138)}
def rgb(t): return ','.join(map(str, t))
def pal(n):
    """n distinct colours, evenly spaced hues, two lightness levels so neighbours differ."""
    out = []
    for i in range(n):
        h = i / n; l = 0.42 if i % 2 == 0 else 0.58; r, g, b = colorsys.hls_to_rgb(h, l, 0.75); out.append((int(r * 255), int(g * 255), int(b * 255)))
    return out
for G in tags:
    D = os.path.join(a.dest, G); os.makedirs(D, exist_ok=True)
    fa = f'{a.ann}/fClaHyb_{G}/fClaHyb_{G}_genome/fClaHyb_{G}_LG.fa'; sizes = {l.split()[0]: int(l.split()[1]) for l in open(fa + '.fai')}
    cat = list(csv.DictReader(open(f'{a.final}/cat59/59_catalogue_{G}.tsv'), delimiter='\t'))
    g59 = {r['family']: r for r in csv.DictReader(open(f'{a.final}/unit_content_2026-11-07/59g_unit_content_{G}.tsv'), delimiter='\t')}
    met = {r['family']: r for r in csv.DictReader(open(f'{a.final}/family_metrics_2026-11-07/59h_family_metrics_{G}.tsv'), delimiter='\t')}
    fam = lambda r: r['family_v46'] or r['family']
    fbp = collections.Counter()
    for r in cat:
        if r['element_class'] == 'satellite': fbp[fam(r)] += int(r['bp'])
    top = [f for f, _ in fbp.most_common(12)]; col = dict(zip(top, pal(12)))
    cat.sort(key=lambda r: (r['chrom'], int(r['start0'])))
    with open(f'{D}/satDNA_v47.runs.gff3', 'w') as gff, open(f'{D}/satDNA_v47.runs_by_family.bed', 'w') as bf, open(f'{D}/satDNA_v47.runs_by_decision.bed', 'w') as bd, \
         open(f'{D}/satDNA_v47.gene_arrays.bed', 'w') as bg, open(f'{D}/satDNA_v47.largest_satellite_cluster.bed', 'w') as bc:
        gff.write('##gff-version 3\n'); bf.write(f'track name="satDNA v4.7 runs by family ({G})" description="SatChain v4.7 runs of record; colour = family (12 largest satellite families), black = gene array" itemRgb="On" visibility=2\n')
        bd.write(f'track name="satDNA v4.7 runs by 59g decision ({G})" description="dark red satellite / black gene array / orange gene-derived satellite / purple unresolved" itemRgb="On" visibility=2\n')
        bg.write(f'track name="tandem gene arrays v4.7 ({G})" description="runs of families called gene_array by 59b" visibility=2\n'); bc.write(f'track name="largest satellite cluster per chromosome ({G})" description="runs with cen_candidate_ror = yes" visibility=1\n')
        dens = collections.Counter()
        for i, r in enumerate(cat, 1):
            f = fam(r); g = g59.get(f, {}); m = met.get(f, {}); c, s, e = r['chrom'], int(r['start0']), int(r['end']); name = r.get('name_v46') or f
            dec = g.get('unit_content_rule', 'satellite' if r['element_class'] == 'satellite' else 'gene_array')
            ga = r['element_class'].startswith('gene_array')
            attrs = [('ID', f'{G}_R{i:05d}'), ('Name', name), ('family', f), ('element_class', r['element_class']), ('gene_flag', r['gene_flag']), ('unit_content_rule', dec),
                     ('sequence_axis', g.get('sequence_axis', '')), ('unit_organization', g.get('unit_organization', '')), ('gene_class', g.get('gene_best_class', '')), ('gene_identity', g.get('gene_best_ident', '')),
                     ('te_frac_of_unit', g.get('te_frac', '')), ('unit_bp', r['unit_bp']), ('unit_class', r['unit_class']), ('dominant_class', r['dominant_class']), ('n_confident_copies', r['n_conf']),
                     ('mean_identity', r['mean_identity_conf']), ('reverse_fraction', r['reverse_fraction']), ('gates_passed', r['gates_passed']), ('tier', r['tier']), ('seq_cluster', r['seq_cluster']),
                     ('putative_centromeric_array_L1', r['putative_centromeric_array_L1']), ('cen_candidate_ror', r['cen_candidate_ror']), ('EDTA_overlap_fraction', r['EDTA_overlap_fraction']),
                     ('family_total_bp', m.get('total_bp', '')), ('family_genome_pct', m.get('genome_pct', '')), ('family_n_loci', m.get('n_loci', '')), ('family_divergence_pct', m.get('divergence_pct', '')),
                     ('homolog_other_haplotype', m.get('homolog_other_haplotype', '')), ('source', r['source']), ('array_id', r['array_id_v46'] or r['array_id'])]
            gff.write(f"{c}\tSatChain_v4.7\t{'tandem_gene_array' if ga else 'satellite_DNA'}\t{s + 1}\t{e}\t.\t.\t.\t" + ';'.join(f'{k}={v}' for k, v in attrs if v not in ('', None)) + '\n')
            rc = (17, 17, 17) if ga else col.get(f, (160, 160, 160))
            bf.write(f'{c}\t{s}\t{e}\t{name}\t{min(1000, int(r["n_conf"]))}\t.\t{s}\t{e}\t{rgb(rc)}\n')
            bd.write(f'{c}\t{s}\t{e}\t{name}|{dec}\t{min(1000, int(r["n_conf"]))}\t.\t{s}\t{e}\t{rgb(DEC.get(dec, (160, 160, 160)))}\n')
            if ga: bg.write(f"{c}\t{s}\t{e}\t{f}:{r['element_class'].split(':', 1)[-1]}\t{r['n_conf']}\t.\n")
            if r['cen_candidate_ror'] == 'yes': bc.write(f'{c}\t{s}\t{e}\t{name}\t{r["n_conf"]}\t.\n')
            if r['element_class'] == 'satellite':
                for w in range(s // 100000, (e - 1) // 100000 + 1):
                    ws, we = w * 100000, (w + 1) * 100000; dens[(c, w)] += min(e, we) - max(s, ws)
    with open(f'{D}/satDNA_v47.density_100kb.bedgraph', 'w') as fh:
        fh.write(f'track type=bedGraph name="satellite bp per 100 kb ({G})" visibility=full color=138,28,43 autoScale=on\n')
        for (c, w), v in sorted(dens.items(), key=lambda kv: (kv[0][0], kv[0][1])): fh.write(f'{c}\t{w * 100000}\t{min(sizes[c], (w + 1) * 100000)}\t{v}\n')
    # loci computed here at base-pair precision (the loci TSV in docs/ is rounded to 10 kb): runs <= 100 kb apart, any class; class by bp majority
    with open(f'{D}/satDNA_v47.loci_top3.bed', 'w') as fh:
        fh.write(f'track name="three largest run clusters per chromosome ({G})" description="rank:class:kb; black = tandem gene array" itemRgb="On" visibility=2\n')
        byc = collections.defaultdict(list)
        for r in cat: byc[r['chrom']].append(r)
        for c in sorted(byc):
            rs = sorted(byc[c], key=lambda r: int(r['start0'])); cl = []; cur = [rs[0]]
            for p, q in zip(rs, rs[1:]):
                if int(q['start0']) - int(p['end']) <= 100000: cur.append(q)
                else: cl.append(cur); cur = [q]
            cl.append(cur); cl.sort(key=lambda L: -sum(int(r['bp']) for r in L))
            for k, L in enumerate(cl[:3], 1):
                bp = sum(int(r['bp']) for r in L); comp = collections.Counter()
                for r in L:
                    lab = r['element_class'] if r['element_class'] != 'satellite' else ('minisat' if r['unit_class'] != 'satellite' else 'satellite')
                    if r['gene_flag'].startswith('gene_related'): lab = 'gene_related'
                    comp[lab] += int(r['bp'])
                t = comp.most_common(1)[0][0]
                lc = 'tandem gene array' if t.startswith('gene_array') else 'gene-related satellite (tRNA-derived)' if t == 'gene_related' else 'minisatellite cluster' if t == 'minisat' else 'satellite array'
                s, e = int(L[0]['start0']), max(int(r['end']) for r in L)
                fh.write(f"{c}\t{s}\t{e}\t{k}:{lc.split(' (')[0].replace(' ', '_')}:{bp // 1000}kb\t{min(1000, bp // 1000)}\t.\t{s}\t{e}\t{rgb(LOC[lc])}\n")
    with open(f'{D}/candidate_centromere_intervals.bed', 'w') as fh:
        fh.write(f'track name="candidate centromere intervals, June synteny BED ({G})" description="support = satellite share inside; location call only" visibility=2 color=0,0,0\n')
        for r in csv.DictReader(open(f'{a.final}/runs_of_record/{G}/centromere_check.tsv'), delimiter='\t'):
            fh.write(f"{r['chrom']}\t{r['cen_start0']}\t{r['cen_end']}\tsupport={float(r['support']):.2f}|{r['verdict']}|{r['dominant_family']}\t{int(float(r['support']) * 1000)}\t.\n")
    tel = f'{a.out45}/fClaHyb_{G}_satellites/11_telomere_blocks.tsv'
    if os.path.exists(tel):
        with open(f'{D}/telomere_blocks.bed', 'w') as fh:
            fh.write(f'track name="telomere blocks >= 500 bp ({G})" visibility=1 color=30,132,73\n')
            for r in csv.DictReader(open(tel), delimiter='\t'):
                if int(r['block_length_bp']) >= 500: fh.write(f"{r['chrom']}\t{r['start_0based']}\t{r['end']}\ttelomere_{r['block_length_bp']}bp\t{min(1000, int(r['block_length_bp']) // 10)}\t.\n")
    shutil.copy(f'{a.final}/family_metrics_2026-11-07/59h_family_metrics_{G}.tsv', f'{D}/family_metrics.tsv'); shutil.copy(f'{a.final}/unit_content_2026-11-07/59g_unit_content_{G}.tsv', f'{D}/unit_content.tsv')
    tracks = ['satDNA_v47.runs_by_decision.bed', 'satDNA_v47.runs_by_family.bed', 'satDNA_v47.runs.gff3', 'satDNA_v47.gene_arrays.bed', 'satDNA_v47.largest_satellite_cluster.bed', 'satDNA_v47.loci_top3.bed',
              'candidate_centromere_intervals.bed', 'telomere_blocks.bed', 'satDNA_v47.density_100kb.bedgraph']
    tracks = [t for t in tracks if os.path.exists(f'{D}/{t}')]
    with open(f'{D}/igv_session.xml', 'w') as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'); fh.write(f'<Session genome="{fa.replace(chr(92), "/")}" locus="All" version="8">\n  <Resources>\n')
        for t in tracks: fh.write(f'    <Resource name="{t}" path="{t}"/>\n')
        fh.write('  </Resources>\n</Session>\n')
    with open(f'{D}/README.md', 'w') as fh:
        fh.write(f"# IGV tracks, SatChain v4.7 final catalogue, {G} haplotype (2026-11-07)\n\nOpen `igv_session.xml` in IGV (File > Open Session). The genome is the annotation FASTA `{fa}` (its .fai is beside it); the track paths are relative to this folder.\n\n"
                 "| file | what | colour |\n|---|---|---|\n"
                 "| satDNA_v47.runs_by_decision.bed | every run of record, name = family name and the 59g decision | dark red satellite, black gene array, orange gene-derived satellite, purple unresolved |\n"
                 "| satDNA_v47.runs_by_family.bed | the same runs, name = family name | 12 largest satellite families distinct, others grey, gene arrays black |\n"
                 "| satDNA_v47.runs.gff3 | the same runs with every attribute (click a feature): unit, class, copies, identity, tier, cluster, centromere flags, 59g axes, family metrics, homolog | IGV default |\n"
                 "| satDNA_v47.gene_arrays.bed | runs of the families 59b calls gene arrays; name = family:gene class | |\n"
                 "| satDNA_v47.largest_satellite_cluster.bed | runs in the largest satellite cluster of their chromosome (cen_candidate_ror) | |\n"
                 "| satDNA_v47.loci_top3.bed | the three largest run clusters per chromosome; name = rank:class:kb | by locus class; black = tandem gene array |\n"
                 "| candidate_centromere_intervals.bed | the June synteny-reconciled candidate BED (= the step-31 interval); name = satellite share, verdict, dominant family | black |\n"
                 "| telomere_blocks.bed | v4.5 step-11 telomere blocks >= 500 bp | green |\n"
                 "| satDNA_v47.density_100kb.bedgraph | satellite bp per 100-kb window (satellite-class runs only) | dark red |\n\n"
                 "Evidence level L1 (sequence location). The candidate centromere intervals are location calls, not validated centromeres. Gene arrays are kept apart from satellite totals. "
                 "The 59g decision (satellite / gene_array / gene_derived_satellite / unresolved) is per family; `unresolved` marks families where annotation and unit sequence disagree or the 59b ratio is threshold-sensitive.\n")
    print(G, 'runs', len(cat), 'tracks', len(tracks), '->', D)
