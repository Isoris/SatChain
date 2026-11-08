#!/bin/bash
# Mac rerun (2026-10-22): writing the BLAST DB and hits to F: through WSL failed (makeblastdb incomplete, "BLAST failed to write output"); build and search on the WSL disk, copy results to F:.
source ~/mambaforge/etc/profile.d/conda.sh; conda activate assembly
W=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22; T=~/ltr_tmp; mkdir -p $T
makeblastdb -in /mnt/f/data/fClaHyb_func_annotation_v1.0/fClaHyb_Mac/fClaHyb_Mac_genome/fClaHyb_Mac_LG.fa -dbtype nucl -out $T/Mac -parse_seqids > $T/mkdb.log 2>&1 || { echo "makeblastdb failed"; tail -3 $T/mkdb.log; exit 1; }
blastn -task dc-megablast -query $W/seq/TE_00003409_library.fa -db $T/Mac -evalue 1e-10 -num_threads 2 -max_target_seqs 100 -max_hsps 100000 -outfmt "6 sseqid sstart send qstart qend pident length evalue bitscore" > $T/Mac.hits.tsv && echo "Mac HSPs: $(wc -l < $T/Mac.hits.tsv)" && cp $T/Mac.hits.tsv $W/copies/Mac.hits.tsv && python3 $W/code/02b_summary.py
