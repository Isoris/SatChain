set -e
D=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22/seq
mkdir -p ~/ltr_tmp && cd ~/ltr_tmp
python3 - <<'P'
seqs=open('/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22/seq/dispersed_fulllength.fa').read().split('>')[1:]
for sp in ('Gar','Mac'):
    open(f'{sp}_full.fa','w').write(''.join('>'+s for s in seqs if s.startswith(sp)))
open('all_full.fa','w').write(''.join('>'+s for s in seqs))
P
for f in Gar_full Mac_full all_full; do mafft --thread 2 --quiet --adjustdirection --localpair --maxiterate 100 $f.fa > $f.aln.fa; done
cp Gar_full.aln.fa Mac_full.aln.fa all_full.aln.fa $D/
grep -c ">" $D/*_full.aln.fa
