set -e
S=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22/seq; D=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22/domains
cd ~/ltr_tmp
head -2 $S/Gar_full_orfs_minus.faa > orf.faa
cp $D/refdomains.faa .; makeblastdb -in refdomains.faa -dbtype prot -out refd >/dev/null
echo "== ORF (1,027 aa) vs reference domains (blastp, e<=1e-3)"
blastp -query orf.faa -db refd -evalue 1e-3 -outfmt "6 sseqid qstart qend sstart send pident evalue bitscore" -max_target_seqs 400 | sort -k8,8nr | awk '{split($1,a,"_"); k=$1; sub(/[:|].*/,"",k)} !seen[$2" "$3]++' | head -25
echo "== source consensus vs EDTA consensus (blastn)"
cp $S/Gar_full_consensus.fa src.fa; cp $S/TE_00003409_library.fa edta.fa
blastn -query src.fa -subject edta.fa -task dc-megablast -evalue 1e-10 -outfmt "6 qstart qend sstart send pident length evalue" | sort -k1,1n
