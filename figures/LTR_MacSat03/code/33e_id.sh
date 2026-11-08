set -e
W=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22; O=$W/split; A=/mnt/f/data/fClaHyb_func_annotation_v1.0/fClaHyb_Mac/fClaHyb_Mac_repeats_EDTA
cd ~/ltr_tmp; cp $O/flanks.fa fl.fa
python3 -c "
s=open('$A/fClaHyb_Mac.TElib.fa').read().split('>')[1:]
open('telib_noself.fa','w').write(''.join('>'+x for x in s if not x.startswith('TE_00003409')))"
cp /mnt/f/data/satellite_FISH/pipeline_v4_5_fixed/out/fClaHyb_Mac_satellites/04_satellite_family_consensuses.fa satfam.fa
F="6 qseqid sseqid qstart qend sstart send pident length evalue bitscore"
blastn -query fl.fa -subject telib_noself.fa -task dc-megablast -evalue 1e-10 -outfmt "$F" > fl_vs_te.tsv
blastn -query fl.fa -subject satfam.fa -task blastn -evalue 1e-5 -outfmt "$F" > fl_vs_sat.tsv
blastn -query fl.fa -subject $W/seq/TE_00003409_library.fa -task dc-megablast -evalue 1e-10 -outfmt "$F" > fl_vs_edta_unit.tsv
cp fl_vs_*.tsv $O/; wc -l fl_vs_*.tsv
