W=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22; cd ~/ltr_tmp
python3 -c "
s=open('telib_noself.fa').read().split('>')[1:]
open('cand.fa','w').write(''.join('>'+x for x in s if x.split('#')[0] in ('TE_00004794_INT','TE_00003814_LTR','TE_00004794_LTR','TE_00003814_INT')))"
grep ">" cand.fa; python3 -c "
for x in open('cand.fa').read().split('>')[1:]: print(x.split()[0], len(''.join(x.split('\n')[1:])))"
echo "== vs EDTA unit consensus TE_00003409"; blastn -query cand.fa -subject $W/seq/TE_00003409_library.fa -task dc-megablast -evalue 1e-10 -outfmt "6 qseqid qstart qend sstart send pident length" | sort -k1,1 -k2,2n
echo "== vs source RTE"; blastn -query cand.fa -subject $W/seq/Gar_full_consensus.fa -task dc-megablast -evalue 1e-10 -outfmt "6 qseqid qstart qend sstart send pident length"
echo "== protein domains (tblastn refdomains vs candidates)"; tblastn -query $W/domains/refdomains.faa -subject cand.fa -evalue 1e-5 -outfmt "6 qseqid sseqid sstart send pident evalue bitscore" | sort -k7,7nr | head -8
