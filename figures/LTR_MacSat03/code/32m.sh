S=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22/seq; cd ~/ltr_tmp; cp $S/units_all.fa ua.fa
blastn -query gsrc.fa -subject ua.fa -task dc-megablast -evalue 1e-10 -outfmt "6 sseqid qstart qend sstart send pident length bitscore" > src_vs_units.tsv
cp src_vs_units.tsv /mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22/homologous/; wc -l src_vs_units.tsv
