#!/bin/bash
# 2026-10-23 redesign: the full consensus contains ~2 lineage A satellite monomers (2,929-3,388) and internal repeats, which inflated copy chaining.
# Use the RT-bearing core 5,560-7,400 (satellite- and repeat-free) to define copies: collinear hits covering >= 70% of the core, span <= 1.5x core.
# Then: counts and centromere-lineage table per haplotype, extraction of core copies, MAFFT, IQ-TREE (GTR+F+G4, UFBoot 1000), tree statistics.
set -eo pipefail
source ~/mambaforge/etc/profile.d/conda.sh; conda activate assembly
W=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22; G=/mnt/f/data/fClaHyb_func_annotation_v1.0; T=~/ltr_tmp/core; mkdir -p $T $W/core
samtools faidx $W/seq/TE_00003409_library.fa "TE_00003409_INT#LTR/unknown:5560-7400" | sed '1s/.*/>element_core_5560_7400/' > $W/seq/element_core.fa
python3 - "$W" <<'PY'
import sys; W = sys.argv[1]; s = "".join(l.strip() for l in open(W + "/seq/TE_00003409_library.fa") if not l.startswith(">"))
s = s[:2928] + "N" * (3388 - 2928) + s[3388:]
open(W + "/seq/element_masked_satellite.fa", "w").write(">element_consensus_satellite_masked_2929_3388\n" + s + "\n")
PY
for h in Gar Mac; do
  # reuse the complete DBs from 2026-10-22 (rebuilding fails under WSL memory pressure): Gar on F:, Mac on the WSL disk
  DBP=$([ $h = Gar ] && echo $W/copies/db/Gar || echo ~/ltr_tmp/Mac)
  blastn -task dc-megablast -query $W/seq/element_core.fa -db $DBP -evalue 1e-10 -num_threads 2 -max_target_seqs 100 -max_hsps 100000 \
    -outfmt '6 sseqid sstart send qstart qend pident length evalue bitscore' > $T/$h.core_hits.tsv
  cp $T/$h.core_hits.tsv $W/core/; echo "$h core HSPs: $(wc -l < $T/$h.core_hits.tsv)"
done
python3 $W/code/06b_core_copies.py
: > $T/core_copies.fa
for h in Gar Mac; do
  while IFS=$'\t' read -r reg st name; do
    if [ "$st" = "-" ]; then samtools faidx -i $G/fClaHyb_$h/_genome/fClaHyb_${h}_LG.fa "$reg"; else samtools faidx $G/fClaHyb_$h/_genome/fClaHyb_${h}_LG.fa "$reg"; fi | sed "1s#.*#>$name#"
  done < $W/core/$h.core_regions.tsv >> $T/core_copies.fa
done
echo "core copies extracted: $(grep -c '>' $T/core_copies.fa)"
mafft --auto --thread 2 --quiet $T/core_copies.fa > $T/copies.aln.fa
python3 - $T <<'PY'
import sys
T = sys.argv[1]; names, seqs = [], []
for l in open(f"{T}/copies.aln.fa"):
    if l.startswith(">"): names.append(l[1:].strip()); seqs.append([])
    else: seqs[-1].append(l.strip())
seqs = ["".join(s) for s in seqs]; L = len(seqs[0]); n = len(seqs)
keep = [i for i in range(L) if sum(s[i] == "-" for s in seqs) <= n / 2]
with open(f"{T}/copies.trim.fa", "w") as o:
    for nm, s in zip(names, seqs): o.write(f">{nm}\n{''.join(s[i] for i in keep)}\n")
print(f"alignment {L} columns -> {len(keep)} kept (<= 50% gaps)")
PY
iqtree2 -s $T/copies.trim.fa -m GTR+F+G4 -B 1000 -nt 2 -pre $T/copies -quiet -redo
mkdir -p $W/tree; cp $T/core_copies.fa $T/copies.aln.fa $T/copies.trim.fa $T/copies.treefile $T/copies.iqtree $W/tree/
python3 $W/code/03b_tree_stats.py
