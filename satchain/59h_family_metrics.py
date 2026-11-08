#!/usr/bin/env python3
"""59h_family_metrics.py  (SatChain v4.7, 2026-11-07): the standard per-family satellitome table, from tables that already exist.

Purpose (declared): results table for the manuscript and for comparison with other satellitomes. The metrics are the ones the field
reports per satDNA family (Garrido-Ramos, Plohl & Satovic-Vuksic 2025, Int J Mol Sci 26:11291, sections 3-7): monomer length, genome
share and copy number, number and length distribution of arrays (short 'seed' arrays vs long arrays), chromosomal distribution
(clustered vs dispersed, candidate-centromeric vs elsewhere), intra-family divergence, monomer structure (GC, palindromy), gene / TE
association (59g), and the library view (families shared by the two haplotypes vs haplotype-specific, copy-number ratio). Nothing here
is a new measurement on sequence: every column is an aggregate of the 59 catalogue, the 03 consensus tables and the 59g table.

Per family (family_v46, one row per haplotype it occurs in):
  unit_bp (bp-weighted median of the runs), unit_min / unit_max; consensus GC, self-vs-revcomp identity and palindromy call (03)
  total_bp, genome_pct (of the haplotype assembly), copies (sum of confident copies), copies_from_bp (total_bp / unit_bp)
  n_arrays (runs), array_bp min / median / max, bp share in arrays < 1 kb (seeds), 1-10 kb, >= 10 kb (long); longest array
  n_chromosomes, n_loci (runs <= 100 kb apart), largest_locus_share (bp of the largest locus / total_bp: 1 = one cluster, small = dispersed),
  bp share inside the candidate centromere intervals (putative_centromeric_array_L1) and in the largest satellite cluster (cen_candidate_ror)
  divergence: 100 - bp-weighted mean identity of confident copies to their class consensus; identity min / max over runs
  strand coherence (bp-weighted |reverse_fraction - 0.5| x 2), EDTA overlap, tier, seq_cluster, element_class (59b), unit_content_rule
  and the three axes (59g), gene class, te_frac, othersat_frac
  low_copy_flag: fewer than --min-copies confident copies OR total_bp < --min-bp (the review's 'where to set the threshold' point)
Library table (59h_library.tsv): one row per family id seen in either haplotype: bp and copies in Gar and in Mac, shared / Gar-only /
Mac-only, log2 copy ratio, unit in each haplotype.
Summary (59h_summary.tsv): families, total satellite bp and % per haplotype, most abundant family and its %, families below the low-copy
line, families shared between haplotypes, bp share of the five most abundant families.
Figure (59h_satellitome.png/.pdf, 450 dpi): a = abundance rank curve per haplotype; b = array-length distribution (bp share) per tier.
usage: 59h_family_metrics.py <final cat59 dir> <unit_content dir (59g)> <v4.5 out root> <new03 pattern {G}> <annotation root> <dest>
       [--tags Gar,Mac] [--min-copies 10] [--min-bp 1000] [--locus-gap 100000]
"""
import sys, os, csv, math, collections, argparse
import numpy as np
ap = argparse.ArgumentParser()
ap.add_argument('cat59'); ap.add_argument('g59'); ap.add_argument('out45'); ap.add_argument('new03'); ap.add_argument('ann'); ap.add_argument('dest')
ap.add_argument('--tags', default='Gar,Mac'); ap.add_argument('--min-copies', type=int, default=10); ap.add_argument('--min-bp', type=int, default=1000)
ap.add_argument('--locus-gap', type=int, default=100_000); ap.add_argument('--hom-identity', type=float, default=0.80, help='homolog cut; swept 0.70-0.90 in 59h_library_sweep.tsv')
a = ap.parse_args(); os.makedirs(a.dest, exist_ok=True); tags = a.tags.split(','); HOM_ALL = {}; ROWS = {}
def num(v, d=0.0):
    try: return float(v)
    except Exception: return d
def wmedian(vals, w):
    o = np.argsort(vals); v = np.asarray(vals)[o]; ww = np.asarray(w)[o]; c = np.cumsum(ww); return float(v[np.searchsorted(c, c[-1] / 2)])
rows_all = []; lib = collections.defaultdict(dict); summ = []
for G in tags:
    fai = f'{a.ann}/fClaHyb_{G}/fClaHyb_{G}_genome/fClaHyb_{G}_LG.fa.fai'; gsize = sum(int(l.split()[1]) for l in open(fai))
    cat = list(csv.DictReader(open(f'{a.cat59}/59_catalogue_{G}.tsv'), delimiter='\t'))
    g59 = {r['family']: r for r in csv.DictReader(open(f'{a.g59}/59g_unit_content_{G}.tsv'), delimiter='\t')} if os.path.exists(f'{a.g59}/59g_unit_content_{G}.tsv') else {}
    cons = {}
    for p in (f'{a.out45}/fClaHyb_{G}_satellites/03_array_consensus.tsv', f'{a.new03.format(G=G)}/03_array_consensus.tsv'):
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter='\t'): cons[(r['array_id'], r['monomer_class'])] = r
    fam = collections.defaultdict(list)
    for r in cat: fam[r['family_v46'] or r['family']].append(r)
    # homolog in the OTHER haplotype from the 59g alignments (unit vs every family consensus of both haplotypes): family ids are
    # haplotype-specific by naming (GarSatNN / MacSatNN), so 'shared' must come from sequence, not from the id
    other = 'MacSat' if G == 'Gar' else 'GarSat'; HOM_ALL[G] = collections.defaultdict(list)
    hp = f'{a.g59}/59g_hits_{G}.tsv'
    if os.path.exists(hp):
        for h in csv.DictReader(open(hp), delimiter='\t'):
            # 2026-11-08: the alignment must cover >= 50 % of BOTH units (unit_frac = share of this unit, ref_cov = share of the other
            # family's consensus); a shared tRNA / snRNA segment inside a long unit passed the one-sided rule (GarSat4282 vs MacSat1286)
            if h['ref_group'] != 'SAT' or not h['ref'].startswith(other) or num(h['unit_frac']) < 0.5 or num(h['ref_cov']) < 0.5: continue
            HOM_ALL[G][h['family']].append((h['ref'], num(h['identity']), num(h['unit_frac'])))   # decided in the post-pass
    sat_total = sum(int(r['bp']) for r in cat if r['element_class'] == 'satellite')
    rows = []
    for f, rs in fam.items():
        bp = np.array([int(r['bp']) for r in rs]); tot = int(bp.sum()); units = np.array([num(r['unit_bp'], np.nan) for r in rs]); ok = ~np.isnan(units)
        unit = wmedian(units[ok], bp[ok]) if ok.any() else float('nan')
        big = max(rs, key=lambda r: int(r['bp'])); c = cons.get((big['array_id_v46'] or big['array_id'], big['dominant_class'])) or cons.get((big['array_id'], big['dominant_class'])) or {}
        copies = sum(int(r['n_conf']) for r in rs)
        # loci: runs <= locus-gap apart on one chromosome
        byc = collections.defaultdict(list)
        for r in rs: byc[r['chrom']].append((int(r['start0']), int(r['end']), int(r['bp'])))
        loci = []
        for ch, v in byc.items():
            v.sort(); cur = [v[0]]
            for x in v[1:]:
                if x[0] - cur[-1][1] <= a.locus_gap: cur.append(x)
                else: loci.append(sum(y[2] for y in cur)); cur = [x]
            loci.append(sum(y[2] for y in cur))
        ident = np.array([num(r['mean_identity_conf']) for r in rs]); rev = np.array([num(r['reverse_fraction'], 0.5) for r in rs])
        g = g59.get(f, {})
        row = dict(haplotype=G, family=f, name=big.get('name_v46', ''), element_class_59b=rs[0]['element_class'], unit_content_rule=g.get('unit_content_rule', ''), sequence_axis=g.get('sequence_axis', ''),
                   unit_organization=g.get('unit_organization', ''), gene_class=g.get('gene_best_class', ''), gene_ident=g.get('gene_best_ident', ''), te_frac=g.get('te_frac', ''), othersat_frac=g.get('othersat_frac', ''),
                   tier=rs[0]['tier'], seq_cluster=rs[0]['seq_cluster'],
                   unit_bp=round(unit, 1), unit_min=int(np.nanmin(units)) if ok.any() else '', unit_max=int(np.nanmax(units)) if ok.any() else '',
                   consensus_GC=c.get('consensus_GC_percent', ''), self_vs_revcomp_identity=c.get('self_vs_revcomp_identity', ''), palindromy=c.get('palindromy_call', ''),
                   total_bp=tot, genome_pct=round(100 * tot / gsize, 4), satellite_share_pct=round(100 * tot / sat_total, 3) if rs[0]['element_class'] == 'satellite' else '',
                   copies=copies, copies_from_bp=round(tot / unit, 1) if unit and not math.isnan(unit) else '',
                   n_arrays=len(rs), array_bp_min=int(bp.min()), array_bp_median=int(np.median(bp)), array_bp_max=int(bp.max()),
                   bp_share_arrays_lt1kb=round(float(bp[bp < 1000].sum()) / tot, 3), bp_share_arrays_1to10kb=round(float(bp[(bp >= 1000) & (bp < 10000)].sum()) / tot, 3), bp_share_arrays_ge10kb=round(float(bp[bp >= 10000].sum()) / tot, 3),
                   n_chromosomes=len(byc), n_loci=len(loci), largest_locus_share=round(max(loci) / tot, 3), loci_ge10kb=sum(1 for x in loci if x >= 10000),
                   bp_share_in_cen_interval=round(sum(int(r['bp']) for r in rs if r['putative_centromeric_array_L1'] == 'yes') / tot, 3), bp_share_in_cen_ror=round(sum(int(r['bp']) for r in rs if r['cen_candidate_ror'] == 'yes') / tot, 3),
                   divergence_pct=round(100 - float(np.average(ident, weights=bp)), 2), identity_min=round(float(ident.min()), 1), identity_max=round(float(ident.max()), 1),
                   strand_coherence=round(float(np.average(np.abs(rev - 0.5) * 2, weights=bp)), 3), edta_overlap=round(float(np.average([num(r['EDTA_overlap_fraction']) for r in rs], weights=bp)), 3),
                   low_copy_flag='yes' if copies < a.min_copies or tot < a.min_bp else 'no',
                   homolog_other_haplotype='', homolog_identity='', homolog_unit_frac='', homolog_status='')
        rows.append(row); lib[f][G] = (tot, copies, unit)
    rows.sort(key=lambda r: -r['total_bp']); ROWS[G] = rows
    sat = [r for r in rows if r['element_class_59b'] == 'satellite']; top = sat[0]
    summ.append(dict(haplotype=G, assembly_bp=gsize, families_all=len(rows), families_satellite=len(sat), satellite_bp=sat_total, satellite_pct=round(100 * sat_total / gsize, 3),
                     gene_array_families=sum(1 for r in rows if r['element_class_59b'].startswith('gene_array')), most_abundant=top['family'], most_abundant_pct=top['genome_pct'],
                     top5_share_of_satellite_pct=round(sum(r['total_bp'] for r in sat[:5]) * 100 / sat_total, 1), families_low_copy=sum(1 for r in sat if r['low_copy_flag'] == 'yes'),
                     families_ge_0_01pct=sum(1 for r in sat if r['genome_pct'] >= 0.01), families_single_locus=sum(1 for r in sat if r['n_loci'] == 1), families_dispersed_ge5loci=sum(1 for r in sat if r['n_loci'] >= 5),
                     median_divergence_pct=round(float(np.median([r['divergence_pct'] for r in sat])), 2)))
    rows_all += rows
    print(G, 'families', len(rows), 'satellite', len(sat), f'{sat_total / 1e6:.2f} Mb = {100 * sat_total / gsize:.2f} %', 'top', top['family'], top['genome_pct'], '% | low-copy', summ[-1]['families_low_copy'], '| single-locus', summ[-1]['families_single_locus'])
# ---- post-pass: homolog in the other haplotype, decided now that both catalogues are known
# shared = best homolog WITH a run of record at identity >= cut (alignment covers >= 50 % of both units); when the only homologs are v4.5 consensuses with no
# run in the final catalogue -> homolog_without_run_of_record; else haplotype-specific. The identity cut is swept (59h_library_sweep.tsv).
bp_of = {(G, f): v for f, d in lib.items() for G, v in d.items()}
def pick(G, f, cut):
    oth = [t for t in tags if t != G][0]; c = HOM_ALL.get(G, {}).get(f, [])
    with_run = [x for x in c if (oth, x[0]) in bp_of and x[1] >= cut]; anyh = [x for x in c if x[1] >= cut]
    if with_run: return 'shared', max(with_run, key=lambda x: x[1])
    if anyh: return 'homolog_without_run_of_record', max(anyh, key=lambda x: x[1])
    return f'{G}_only', None
HOM = {}
for G in tags:
    for r in ROWS[G]:
        st_, h = pick(G, r['family'], a.hom_identity); HOM[(G, r['family'])] = (st_, h)
        r['homolog_status'] = st_; r['homolog_other_haplotype'] = h[0] if h else ''; r['homolog_identity'] = round(h[1], 3) if h else ''; r['homolog_unit_frac'] = round(h[2], 3) if h else ''
    with open(f'{a.dest}/59h_family_metrics_{G}.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(ROWS[G][0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(ROWS[G])
with open(f'{a.dest}/59h_library.tsv', 'w') as fh:
    fh.write('haplotype\tfamily\tstatus\thomolog\thomolog_identity\treciprocal\tbp\tbp_homolog\tcopies\tcopies_homolog\tunit\tunit_homolog\tlog2_copy_ratio_self_over_homolog\n'); st = collections.Counter()
    for G in tags:
        oth = [t for t in tags if t != G][0]
        for r in ROWS[G]:
            f = r['family']; st_, h = HOM[(G, f)]; tot, copies, unit = lib[f][G]; st[(G, st_)] += 1
            if st_ == 'shared':
                hv = bp_of[(oth, h[0])]; b = HOM.get((oth, h[0]), ('', None))[1]; back = 'yes' if b and b[0] == f else 'no'
                fh.write(f"{G}\t{f}\t{st_}\t{h[0]}\t{h[1]:.3f}\t{back}\t{tot}\t{hv[0]}\t{copies}\t{hv[1]}\t{unit:.1f}\t{hv[2]:.1f}\t{math.log2((copies + 1) / (hv[1] + 1)):.2f}\n")
            elif st_ == 'homolog_without_run_of_record':
                fh.write(f"{G}\t{f}\t{st_}\t{h[0]}\t{h[1]:.3f}\t\t{tot}\t\t{copies}\t\t{unit:.1f}\t\t\n")
            else: fh.write(f"{G}\t{f}\t{st_}\t\t\t\t{tot}\t\t{copies}\t\t{unit:.1f}\t\t\n")
for s_ in summ:
    G = s_['haplotype']; tot = sum(v[0] for (g, f), v in bp_of.items() if g == G)
    s_['families_with_homolog_in_other_haplotype'] = st[(G, 'shared')]; s_['families_homolog_without_run_of_record'] = st[(G, 'homolog_without_run_of_record')]; s_['families_haplotype_specific'] = st[(G, f'{G}_only')]
    s_['satellite_bp_share_shared'] = round(sum(lib[f][G][0] for (g, f), (st_, h) in HOM.items() if g == G and st_ == 'shared') * 100 / tot, 1)
with open(f'{a.dest}/59h_summary.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(summ[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(summ)
print('library', dict(st))
with open(f'{a.dest}/59h_library_sweep.tsv', 'w') as fh:
    fh.write('haplotype\tidentity_cut\tfamilies_shared\tsatellite_bp_share_shared_pct\tGarSat04\tGarSat09\tGarSat03\tMacSat03\tMacSat06\n')
    for G in tags:
        tot = sum(v[0] for (g, f), v in bp_of.items() if g == G)
        for cut in (0.70, 0.75, 0.80, 0.85, 0.90):
            res = {r['family']: pick(G, r['family'], cut) for r in ROWS[G]}
            n = sum(1 for v in res.values() if v[0] == 'shared'); share = sum(lib[f][G][0] for f, v in res.items() if v[0] == 'shared') * 100 / tot
            def ex(f):
                if f not in res: return '-'
                return f"{res[f][1][0]} {res[f][1][1]:.2f}" if res[f][0] == 'shared' else res[f][0]
            fh.write(f"{G}\t{cut}\t{n}\t{share:.1f}\t{ex('GarSat04')}\t{ex('GarSat09')}\t{ex('GarSat03')}\t{ex('MacSat03')}\t{ex('MacSat06')}\n")
# ---- figure
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, 2, figsize=(10, 4)); fig.patch.set_facecolor('white')
for G, col in zip(tags, ('#8a1c2b', '#1f2d4d')):
    sat = sorted([r for r in rows_all if r['haplotype'] == G and r['element_class_59b'] == 'satellite'], key=lambda r: -r['total_bp'])
    axs[0].plot(range(1, len(sat) + 1), [r['genome_pct'] for r in sat], color=col, lw=1.2, label=f'{G} ({len(sat)} families)')
axs[0].set_yscale('log'); axs[0].set_xscale('log'); axs[0].set_xlabel('family rank by abundance'); axs[0].set_ylabel('% of the haplotype assembly'); axs[0].set_title('a  abundance rank curve (satellite families)', fontsize=10, loc='left'); axs[0].legend(frameon=False, fontsize=8)
bins = [0, 500, 1000, 2000, 5000, 10000, 20000, 50000, 1e9]; labels = ['<0.5', '0.5-1', '1-2', '2-5', '5-10', '10-20', '20-50', '>=50 kb']
for i, (G, col) in enumerate(zip(tags, ('#8a1c2b', '#1f2d4d'))):
    cat = list(csv.DictReader(open(f'{a.cat59}/59_catalogue_{G}.tsv'), delimiter='\t')); bp = np.array([int(r['bp']) for r in cat if r['element_class'] == 'satellite'])
    h = [bp[(bp >= lo) & (bp < hi)].sum() / bp.sum() for lo, hi in zip(bins[:-1], bins[1:])]
    axs[1].bar(np.arange(len(h)) + (i - 0.5) * 0.4, h, width=0.4, color=col, label=G)
axs[1].set_xticks(range(len(labels))); axs[1].set_xticklabels(labels, fontsize=8); axs[1].set_xlabel('array (run) length'); axs[1].set_ylabel('share of satellite bp'); axs[1].set_title('b  where the satellite bp sit: array-length classes', fontsize=10, loc='left'); axs[1].legend(frameon=False, fontsize=8)
fig.suptitle('Satellitome metrics of the two fClaHyb haplotypes (SatChain 59h, final v4.7 catalogue)', fontsize=10, x=0.01, ha='left'); fig.tight_layout()
fig.savefig(f'{a.dest}/59h_satellitome.png', dpi=450); fig.savefig(f'{a.dest}/59h_satellitome.pdf'); print('written', a.dest)
