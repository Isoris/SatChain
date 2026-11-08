#!/usr/bin/env python3
"""runs_of_record.py (v4.7, 2026-11-01; was satellite_FISH/runs_of_record_v4.6_2026-10-27/code/runs_of_record_v46.py, now in the repo
because the manuscript figures and tables read this format). From the 59 catalogue of each haplotype (after 59d: family_v46, name_v46,
tier; after 59f: seq_cluster) write the runs-of-record tables the downstream scripts expect:
  <dest>/<G>/runs_of_record.tsv  (haplotype, run_id, chrom, start0, end, run_bp, family, unit_bp, unit_class, copies, tier, verdict, fish_site_ok,
                                  source, array_id, name_v46, seq_cluster, element_class, gene_flag, centromeric_L1, cen_candidate_ror, carried_from)
  <dest>/<G>/runs_of_record.bed, family_totals.tsv, chrom_totals.tsv, centromere_check.tsv; <dest>/satellite_family_nomenclature.tsv
unit_class = catalogue unit_class (100-bp rule) with gene arrays as 'gene_array' (element_class); the optional array-levels table of
benchmarks/array_definition_2026-10-27 (rlvl1/rlvl2) is joined when given (--levels DIR) and kept as columns only.
Satellite bp of record = sum of satellite-class runs (overlaps were resolved in 59, so sum = union).
usage: runs_of_record.py <cat59 dir> <func root (fClaHyb_func_annotation_v1.0)> <dest> [--tags Gar,Mac] [--levels DIR] [--version v4.7]"""
import os, csv, collections, argparse, pandas as pd
ap = argparse.ArgumentParser(); ap.add_argument('CAT'); ap.add_argument('R0'); ap.add_argument('OUT'); ap.add_argument('--tags', default='Gar,Mac'); ap.add_argument('--levels', default=''); ap.add_argument('--version', default='v4.7')
a = ap.parse_args(); nom = []


for h in a.tags.split(','):
    D = os.path.join(a.OUT, h); os.makedirs(D, exist_ok=True)
    C = pd.read_csv(f'{a.CAT}/59_catalogue_{h}.tsv', sep='\t', dtype={'consensus_GC': str})
    fam = C['family_v46'] if 'family_v46' in C else C['family']; name = C['name_v46'] if 'name_v46' in C else fam
    ucls = C.unit_class.where(~C.element_class.str.startswith('gene_array'), 'gene_array')
    R = pd.DataFrame(dict(haplotype=h, run_id=C.run_id, chrom=C.chrom, start0=C.start0, end=C.end, run_bp=C.bp, family=fam, unit_bp=C.unit_bp, unit_class=ucls, copies=C.n_conf,
                          tier='S1', verdict=f'{a.version}_run_of_record', fish_site_ok=(C.bp >= 3000).astype(int), source=f'SatChain_{a.version}', array_id=C.array_id,
                          name_v46=name, seq_cluster=C.get('seq_cluster', 'none'), element_class=C.element_class, gene_flag=C.get('gene_flag', 'none'),
                          centromeric_L1=C.putative_centromeric_array_L1, cen_candidate_ror=C.get('cen_candidate_ror', 'na'), family_tier=C.tier, gates_passed=C.gates_passed,
                          carried_from=C.get('carried_from', 'none')))
    if a.levels and os.path.exists(f'{a.levels}/array_levels_{h}.tsv'):
        L = pd.read_csv(f'{a.levels}/array_levels_{h}.tsv', sep='\t').set_index('array_id'); R['rlvl1'] = C.array_id.map(L.rlvl1); R['rlvl2'] = C.array_id.map(L.rlvl2)
    R.to_csv(os.path.join(D, 'runs_of_record.tsv'), sep='\t', index=False)
    R[['chrom', 'start0', 'end', 'family', 'run_bp']].to_csv(os.path.join(D, 'runs_of_record.bed'), sep='\t', index=False, header=False)
    ft = R.groupby(['family', 'unit_class']).agg(S1_bp=('run_bp', 'sum'), n_S1=('run_bp', 'size')).reset_index(); ft['S2_bp'] = 0; ft['excluded_bp'] = 0; ft['n_S2'] = 0; ft['n_excluded'] = 0
    ft[['family', 'unit_class', 'S1_bp', 'S2_bp', 'excluded_bp', 'n_S1', 'n_S2', 'n_excluded']].sort_values('S1_bp', ascending=False).to_csv(os.path.join(D, 'family_totals.tsv'), sep='\t', index=False)
    Lc = {l.split()[0]: int(l.split()[1]) for l in open(f'{a.R0}/fClaHyb_{h}/fClaHyb_{h}_genome/fClaHyb_{h}_LG.fa.fai')}
    sat = R[R.unit_class == 'satellite']
    with open(os.path.join(D, 'chrom_totals.tsv'), 'w', newline='\n') as o:
        o.write('chrom\tchrom_bp\tsatellite_bp\tsatellite_frac\tn_runs\tmini_bp\tgene_array_bp\n')
        for c, L_ in Lc.items():
            s = sat[sat.chrom == c]; m = R[(R.chrom == c) & (R.unit_class == 'mini')]; g = R[(R.chrom == c) & (R.unit_class == 'gene_array')]
            o.write(f'{c}\t{L_}\t{int(s.run_bp.sum())}\t{s.run_bp.sum() / L_:.5f}\t{len(s)}\t{int(m.run_bp.sum())}\t{int(g.run_bp.sum())}\n')
    with open(os.path.join(D, 'centromere_check.tsv'), 'w', newline='\n') as o:
        o.write('chrom\tcen_start0\tcen_end\tcen_bp\tsatellite_bp_inside\tsupport\tverdict\tdominant_family\n')
        p = f'{a.R0}/fClaHyb_{h}/fClaHyb_{h}_centromeres_satellites/fClaHyb_{h}.centromeres.confidence.bed'
        for l in (open(p) if os.path.exists(p) else []):
            if l.startswith(('track', '#')): continue
            f = l.split('\t'); c, s, e = f[0], int(f[1]), int(f[2]); ins = sat[(sat.chrom == c) & (sat.end > s) & (sat.start0 < e)]
            bp = int(sum(min(e, r.end) - max(s, r.start0) for r in ins.itertuples())); sup = bp / max(1, e - s)
            dom = ins.groupby('family').run_bp.sum().idxmax() if len(ins) else 'none'
            o.write(f'{c}\t{s}\t{e}\t{e - s}\t{bp}\t{sup:.3f}\t{"satellite_supported" if sup >= 0.3 else "weak" if sup >= 0.05 else "unsupported"}\t{dom}\n')
    for f_, g in R.groupby('family'):
        nom.append(dict(haplotype=h, family=f_, name=g.name_v46.iloc[0], element_class=g.element_class.iloc[0], unit_bp_median=int(g.unit_bp.median()), bp=int(g.run_bp.sum()), n_runs=len(g), seq_cluster=g.seq_cluster.iloc[0]))
    print(h, 'runs', len(R), '| satellite (>= 100 bp)', round(sat.run_bp.sum() / 1e6, 3), 'Mb | mini', round(R[R.unit_class == 'mini'].run_bp.sum() / 1e6, 3), '| gene arrays', round(R[R.unit_class == 'gene_array'].run_bp.sum() / 1e6, 3),
          '| centromere check:', pd.read_csv(os.path.join(D, 'centromere_check.tsv'), sep='\t').verdict.value_counts().to_dict())
pd.DataFrame(nom).to_csv(os.path.join(a.OUT, 'satellite_family_nomenclature.tsv'), sep='\t', index=False)
