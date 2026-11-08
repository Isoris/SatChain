S=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22/seq; cd ~/ltr_tmp
cp $S/units_sample.fa units.fa; cp $S/Gar_full_consensus.fa gsrc.fa; cp $S/Mac_full_consensus.fa msrc.fa
for q in gsrc msrc; do
  echo "== $q: source segments found in each centromeric unit (qstart-qend: identity)"
  blastn -query $q.fa -subject units.fa -task dc-megablast -evalue 1e-10 -outfmt "6 sseqid qstart qend pident length" | sort -k1,1 -k2,2n | awk '{a[$1]=a[$1]" "$2"-"$3":"int($4); L[$1]+=$5; I[$1]+=$4*$5} END{for(k in a) printf "%s\tcovered %d\tid %.1f\t%s\n", k, L[k], I[k]/L[k], a[k]}' | sort
done
