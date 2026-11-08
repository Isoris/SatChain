#!/bin/bash
source ~/mambaforge/etc/profile.d/conda.sh; conda activate assembly
W=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22; C=/mnt/f/data/satellite_FISH/pipeline_v4_5_fixed/out/fClaHyb_Mac_satellites/04_satellite_family_consensuses.fa
blastn -query $W/seq/TE_00003409_library.fa -subject $C -evalue 1e-10 -outfmt "6 sseqid pident length qstart qend evalue bitscore" -max_target_seqs 20 > $W/domains/element_vs_Mac_satellite_families.tsv
awk -F"\t" "{c[\$1]+=\$3} END {for (k in c) print c[k], k}" $W/domains/element_vs_Mac_satellite_families.tsv | sort -rn | head -6
awk -F"\t" "{print (\$4<\$5?\$4:\$5), (\$4<\$5?\$5:\$4)}" $W/domains/element_vs_Mac_satellite_families.tsv | sort -n | awk "BEGIN{s=-1;e=-1;t=0} {if(\$1>e){if(e>0)t+=e-s+1; s=\$1;e=\$2} else if(\$2>e) e=\$2} END {t+=e-s+1; print \"bp of the 8,164-bp consensus covered by satellite families:\", t}"
