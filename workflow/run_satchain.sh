#!/usr/bin/env bash
# SatChain, commands in order (experimental). Edit the variables, run step by step; every step is resumable and writes into $WORK.
set -euo pipefail
G=Gar                                    # haplotype tag used in file names (fClaHyb_<G>)
FUNC=/path/to/annotation/fClaHyb_$G      # annotation folder (docs/INPUTS.md)
OVL=/path/to/overlays/$G                 # AnnoSINE + annotation overlays
OUT=/path/to/pipeline/out                # folder holding fClaHyb_<G>_satellites/ (existing arrays: 02, 03, 31 tables)
FA=$FUNC/_genome/fClaHyb_${G}_LG.fa
WORK=work; SC=$(dirname "$0")/../satchain; mkdir -p $WORK
# 57-58: genome-wide annotation-free scan, calibration, new candidates (3-copy chain test, TRF agreement)
python3 $SC/57_genome_tandem_scan.py $FA $G $WORK/scan57/scan57_${G}_job0.tsv all
python3 $SC/58_scan_calibrate.py $G $WORK/scan57 $FUNC $OUT $WORK/hmm56_$G $WORK/cal58_$G
python3 $SC/58b_new_candidates.py $G $WORK/cal58_$G $FUNC $OUT $FA
# 03: decompose new candidates (02-format table written from 58b; see 58c in the docs) into monomer copies
python3 $SC/03_monomer_copies.py --genome $FA --arrays $WORK/new03_$G/02_satellite_arrays.tsv --tag $G --out $WORK/new03_$G --threads 3 --min-array-bp 500
# 56: copy-chain HMM on existing and new arrays, then the per-run sequence check
python3 $SC/56_path_hmm.py --out $OUT/fClaHyb_${G}_satellites --dest $WORK/hmm56_$G
python3 $SC/56_path_hmm.py --out $WORK/new03_$G --dest $WORK/hmm56new_$G
python3 $SC/val56_runfilter.py $G <tandem_signal_1kb.bedgraph> $WORK/hmm56_$G/56_runs_post0.95.bed $OUT/fClaHyb_${G}_satellites/36_characterized_intervals.bed $OUT/fClaHyb_${G}_satellites/02_satellite_arrays.tsv $WORK/hmm56_$G
# 59: catalogue (gates, family, lineage, tiers), gene arrays, TE layer, names
python3 $SC/59_final_catalogue.py $G $OUT $WORK/hmm56_$G $WORK/new03_$G $WORK/hmm56new_$G $WORK/new03_$G/famnew.tsv $WORK/lin55/55d_family_lineage.tsv $FUNC $FA $WORK/cat59
python3 $SC/59b_gene_arrays.py $G $WORK/cat59/59_catalogue_$G.tsv $OVL $WORK/cat59/59b_gene_arrays_$G.tsv
python3 $SC/59d_family_names.py $G $WORK/cat59 $WORK/new03_$G
# 53/53b: segmentation and topology on every catalogue run (runs table: run_id, array_id, chrom, start0, end, unit_bp)
python3 $SC/53_array_segmentation.py --out $OUT/fClaHyb_${G}_satellites --genome $FA --runs $WORK/cat59/runs_${G}_existing.tsv --dest $WORK/cat53/${G}_existing
python3 $SC/53b_topology_compare.py --out $OUT/fClaHyb_${G}_satellites --run $WORK/cat53/${G}_existing --func $FUNC --overlays $OVL --genome-tag $G --runs $WORK/cat59/runs_${G}_existing.tsv
# 55: lineage per copy (sample centromeric segments, cluster k = 5-7, K = 2-6), then family -> lineage
python3 $SC/55a_lineage_sample.py $WORK/lin55 150 "$G=$OUT/fClaHyb_${G}_satellites|$WORK/cat53/${G}_existing"
# 55b: BLAST the sample copies against the 5S genes first -> $WORK/lin55/55_sample_vs_5S.tsv (outfmt 6 + qlen slen)
python3 $SC/55c_lineage_cluster.py $WORK/lin55 && python3 $SC/55d_family_lineage.py $WORK/lin55
