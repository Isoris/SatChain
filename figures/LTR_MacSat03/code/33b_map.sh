set -e
W=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22; O=$W/split; cd ~/ltr_tmp
cp $O/regions.fa reg.fa; cp $W/seq/Gar_full_consensus.fa src.fa; cp $W/insert/insert_rep.fa ins.fa
makeblastdb -in reg.fa -dbtype nucl -out reg >/dev/null
F="6 sseqid sstart send qstart qend pident length evalue bitscore"
blastn -query src.fa -db reg -task dc-megablast -evalue 1e-10 -max_target_seqs 100000 -max_hsps 100000 -outfmt "$F" > src_hits.tsv
blastn -query ins.fa -db reg -task dc-megablast -evalue 1e-10 -max_target_seqs 100000 -max_hsps 100000 -outfmt "$F" > ins_hits.tsv
cp src_hits.tsv ins_hits.tsv $O/; wc -l src_hits.tsv ins_hits.tsv
