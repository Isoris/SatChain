#!/usr/bin/env bash
# SatChain v4.7 (2026-11-01), commands in order. Edit the variables, run step by step; every step writes into $WORK and is resumable.
# Order fixed on 2026-11-01 (audit C09: the old 59 -> 53 -> 55 -> 59 loop is gone): sweep -> copies -> dedup -> HMM (one model per
# haplotype on existing + new tables) -> per-run check -> catalogue (gene arrays inside, before tiers) -> TE layer -> names -> interruptions
# -> sequence clusters (sampled from the catalogue runs, pooled over haplotypes) written back -> topology on the runs of record -> array definition.
set -euo pipefail
G=Gar                                    # haplotype tag used in file names (fClaHyb_<G>); step 55 needs both haplotypes (see below)
FUNC=/path/to/annotation/fClaHyb_$G      # annotation folder, current layout (docs/INPUTS.md; satpaths.py resolves the old one too)
OVL=/path/to/annotation_overlays/$G      # rRNA / tRNA / miRNA / Rfam / AnnoSINE overlays
OUT=/path/to/pipeline/out                # folder holding fClaHyb_<G>_satellites/ (existing arrays: 02, 03, 31, 36 tables)
INS54=/path/to/step54_insertions         # 54b_ins_members.tsv, 54c_ins_families_classified.tsv (59c, 59e)
FA=$FUNC/fClaHyb_${G}_genome/fClaHyb_${G}_LG.fa
GAPS=$FUNC/fClaHyb_${G}_genome/fClaHyb_${G}_LG.gaps.detgaps.bed   # assembly gaps and contig joins are hard stops (56, 53b, 59e, 61, 61b, 64)
SIG=/path/to/igv/${G}_tandem_signal_1kb.bedgraph                   # step-48 track (val56 --signal; prefer --scan57 when the 57 tables exist)
WORK=work; SC=$(dirname "$0")/../satchain; mkdir -p $WORK; export SATCHAIN_GAPS=$GAPS
# 57-58: genome-wide annotation-free scan, calibration (null on lag >= 7 windows outside arrays and TRF), candidates (3-copy chain test,
#        TRF coverage), 58c = the 02-format table of kept candidates (rule new_array OR (7-39 bp AND TRF cover >= 0.5); restored 2026-11-01)
python3 $SC/57_genome_tandem_scan.py $FA $G $WORK/scan57/scan57_${G}_job0.tsv all
python3 $SC/58_scan_calibrate.py $G $WORK/scan57 $FUNC $OUT $WORK/hmm56_$G $WORK/cal58_$G --null lag7      # hmm56 recall column is NA on a first pass
python3 $SC/58b_new_candidates.py $G $WORK/cal58_$G $FUNC $OUT $FA
python3 $SC/58c_new_arrays_02format.py $G $WORK/cal58_$G/58b_new_candidates.bed $WORK/new03_$G/02_satellite_arrays.tsv --route chain_or_trf --trf-field cover
# 03: decompose the new candidates into monomer copies (StringDecomposer; duplicates of the 5-kb batch overlap removed inside 03 since 2026-11-01)
python3 $SC/03_monomer_copies.py --genome $FA --arrays $WORK/new03_$G/02_satellite_arrays.tsv --tag $G --out $WORK/new03_$G --threads 3 --min-array-bp 500
# 03b: tables made before 2026-11-01 carry StringDecomposer batch duplicates (audit CB01): dedup them once into a fresh folder
#   python3 $SC/03b_dedup_copies.py $OUT/fClaHyb_${G}_satellites/03_monomer_copies.tsv <fresh out>/fClaHyb_${G}_satellites/03_monomer_copies
# 56: copy-chain HMM, ONE model per haplotype fitted on the existing and the new tables together (HB05); runs stop at gaps and joins
python3 $SC/56_path_hmm.py --out $OUT/fClaHyb_${G}_satellites $WORK/new03_$G --dest $WORK/hmm56_$G --gaps $GAPS
# val56: per-run sequence check for every run (existing and new) against a background outside all candidate intervals; writes the benchmark table
python3 $SC/val56_runfilter.py --tag $G --runs $WORK/hmm56_$G/56_runs_post0.95.bed --bed36 $OUT/fClaHyb_${G}_satellites/36_characterized_intervals.bed \
  --a02 $OUT/fClaHyb_${G}_satellites/02_satellite_arrays.tsv --dest $WORK/hmm56_$G --scan57 $WORK/scan57/scan57_${G}_job0.tsv --exclude $WORK/cal58_$G/58b_new_candidates.bed
# famnew: family of each new class consensus (BLAST against the family consensuses of both haplotypes; see famnew.py)
# 59: catalogue (overlaps resolved, gates + sweep, gene arrays called before tiers, cen flags); no cluster table yet ('-')
python3 $SC/59_final_catalogue.py $G $OUT $WORK/hmm56_$G $WORK/new03_$G - $WORK/new03_$G/famnew.tsv - $FUNC $FA $WORK/cat59 --overlays $OVL
python3 $SC/59c_te_layer.py $G $WORK/cat59/59_catalogue_$G.tsv $OUT $WORK/new03_$G $FUNC $INS54 $WORK/cat59       # SATCHAIN_59C_CARRY=min for the v4.6 rule
python3 $SC/59d_family_names.py $G $WORK/cat59 $WORK/new03_$G
python3 $SC/59e_interruptions.py $G $WORK/cat59 $FUNC $INS54 --gaps $GAPS --sweep
# 55: sequence clusters of the satellite monomers, POOLED over haplotypes: run once with every haplotype's catalogue, then write back per haplotype
python3 $SC/55a_cluster_sample.py $WORK/lin55 150 "Gar=$OUT|$WORK/cat59/59_catalogue_Gar.tsv|$WORK/new03_Gar" "Mac=$OUT|$WORK/cat59/59_catalogue_Mac.tsv|$WORK/new03_Mac"
python3 $SC/55b_gene_overlap.py $WORK/lin55 Gar=$(dirname $OVL)/Gar Mac=$(dirname $OVL)/Mac
python3 $SC/55c_cluster_kmer.py $WORK/lin55 && python3 $SC/55d_family_cluster.py $WORK/lin55
python3 $SC/59f_apply_clusters.py $G $WORK/cat59 $WORK/lin55/55d_family_cluster.tsv
# 53/53b: segmentation and topology on the runs of record (run mode is the only mode; 02 extents need --extents-ok)
python3 $SC/53_array_segmentation.py --out $OUT/fClaHyb_${G}_satellites --genome $FA --runs $WORK/cat59/59_catalogue_$G.tsv --dest $WORK/cat53/$G
python3 $SC/53b_topology_compare.py --out $OUT/fClaHyb_${G}_satellites --run $WORK/cat53/$G --func $FUNC --overlays $OVL --genome-tag $G --runs $WORK/cat59/59_catalogue_$G.tsv --gaps $GAPS
python3 $SC/53b_o2_summary.py --topologies $G:existing=$WORK/cat53/$G/53b/53b_topologies.tsv --catalogue $G=$WORK/cat59/59_catalogue_$G.tsv --dest $WORK/cat53
# 61 / 61b: unit boundaries of interrupted arrays (chains built by chains.py on family_v46), composition boundaries (JSD, clipped to the spacer)
python3 $SC/61_unit_boundaries.py --catalogue $WORK/cat59/59_catalogue_$G.tsv \
  --copies $OUT/fClaHyb_${G}_satellites/03_monomer_copies.tsv $WORK/new03_$G/03_monomer_copies.tsv \
  --consensus $OUT/fClaHyb_${G}_satellites/03_array_consensus.tsv $WORK/new03_$G/03_array_consensus.tsv \
  --genome $FA --dest $WORK/units61 --tag $G --gaps $GAPS
python3 $SC/61b_jsd_boundaries.py --catalogue $WORK/cat59/59_catalogue_$G.tsv --genome $FA --gaps $GAPS --dest $WORK/jsd61b --tag $G
# 62-65 (exploratory until their sweeps are on disk): see docs/STEP_METHOD.md section 8 and each script's usage line
# benchmarks: python3 $SC/benchmark_band_test.py $OUT/fClaHyb_${G}_satellites $FA $G $WORK/cat59/59_catalogue_$G.tsv $FUNC/fClaHyb_${G}_TR/TRF/fClaHyb_${G}.TRF.bed.gz benchmarks/threshold_band_test_$G
