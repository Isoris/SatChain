#!/bin/bash
source ~/mambaforge/etc/profile.d/conda.sh; conda activate assembly
W=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22
grep ">" $W/domains/orfs.faa
makeblastdb -in $W/seq/TE_00003409_library.fa -dbtype nucl -out $W/domains/element >/dev/null
tblastn -query $W/domains/refdomains.faa -db $W/domains/element -evalue 1e-2 -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen" -num_threads 2 > $W/domains/refdomains_vs_element_tblastn.tsv
sort -t$'\t' -k10,10gr $W/domains/refdomains_vs_element_tblastn.tsv | head -15
