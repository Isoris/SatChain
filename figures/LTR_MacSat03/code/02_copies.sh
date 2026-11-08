#!/bin/bash
# Step 3-4 (2026-10-22): copies of the MacSat03-interleaved RTE-like element in both haplotypes, by sequence (EDTA library IDs differ per genome).
# blastn (dc-megablast) of the 8,164-bp consensus against fClaHyb_Gar_LG and fClaHyb_Mac_LG; 2 threads (shared machine).
# Output: copies/<hap>.hits.tsv (all HSPs, e <= 1e-10), then 02b_summary.py makes copies/<hap>.copies.tsv and copies/summary.txt.
set -eo pipefail
source ~/mambaforge/etc/profile.d/conda.sh; conda activate assembly
W=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22; G=/mnt/f/data/fClaHyb_func_annotation_v1.0
mkdir -p $W/copies/db
for h in Gar Mac; do
  [ -s $W/copies/db/$h.nsq ] || [ -s $W/copies/db/$h.00.nsq ] || makeblastdb -in $G/fClaHyb_$h/_genome/fClaHyb_${h}_LG.fa -dbtype nucl -out $W/copies/db/$h -parse_seqids >/dev/null
  blastn -task dc-megablast -query $W/seq/TE_00003409_library.fa -db $W/copies/db/$h -evalue 1e-10 -num_threads 2 \
    -max_target_seqs 100 -max_hsps 100000 -outfmt '6 sseqid sstart send qstart qend pident length evalue bitscore' > $W/copies/$h.hits.tsv
  echo "$h HSPs: $(wc -l < $W/copies/$h.hits.tsv)"
done
python3 $W/code/02b_summary.py
