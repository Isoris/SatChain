#!/bin/bash
# Step 5 (2026-10-22): tree of the near-full copies (>= 70% of the consensus) of the lineage A-associated RTE-like element from both haplotypes.
# Intermixed Gar/Mac copies = ancestral arrangement; per-species clades = independent amplification. Work on the WSL disk, results copied to F:.
set -eo pipefail
source ~/mambaforge/etc/profile.d/conda.sh; conda activate assembly
W=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22; G=/mnt/f/data/fClaHyb_func_annotation_v1.0; T=~/ltr_tmp/tree; mkdir -p $T $W/tree
python3 - $W $T <<'PY'
import sys, pandas as pd
W, T = sys.argv[1], sys.argv[2]
for h in ("Gar", "Mac"):
    cp = pd.read_csv(f"{W}/copies/{h}.copies.tsv", sep="\t"); f = cp[cp.near_full]
    cen = pd.read_csv(f"/mnt/f/data/satellite_FISH/runs_of_record_v4.5.3_2026-10-18/{h}/centromere_check.tsv", sep="\t").set_index("chrom")
    with open(f"{T}/{h}.regions.tsv", "w") as o:
        for r in f.itertuples():
            inc = r.chrom in cen.index and r.lo < cen.loc[r.chrom, "cen_end"] + 1e6 and r.hi > cen.loc[r.chrom, "cen_start"] - 1e6
            fam = cen.loc[r.chrom, "top_family_in_cen"] if r.chrom in cen.index else "NA"
            name = f"{h}|{r.chrom.split('_')[-1]}|{r.lo}|{'cen' if inc else 'arm'}|{fam}"
            o.write(f"{r.chrom}:{r.lo}-{r.hi}\t{r.strand}\t{name}\n")
PY
: > $T/copies.fa
for h in Gar Mac; do
  while IFS=$'\t' read -r reg st name; do
    if [ "$st" = "-" ]; then samtools faidx -i $G/fClaHyb_$h/_genome/fClaHyb_${h}_LG.fa "$reg"; else samtools faidx $G/fClaHyb_$h/_genome/fClaHyb_${h}_LG.fa "$reg"; fi | sed "1s#.*#>$name#"
  done < $T/$h.regions.tsv >> $T/copies.fa
done
echo "copies: $(grep -c '>' $T/copies.fa)"
mafft --auto --thread 2 --quiet $T/copies.fa > $T/copies.aln.fa
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
cp $T/copies.fa $T/copies.aln.fa $T/copies.trim.fa $T/copies.treefile $T/copies.iqtree $T/copies.log $W/tree/
python3 $W/code/03b_tree_stats.py
