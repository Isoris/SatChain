#!/bin/bash
# 2026-10-23: IQ-TREE on all 1,560 core copies with UFBoot crashed (bus error: WSL ~1 GB free). Subsample: all Mac copies + 400 Gar copies drawn at random,
# stratified by centromere/arm and chromosome (seed 20261023); IQ-TREE -fast with SH-aLRT 1000 (low memory). Full tree to redo when memory allows.
set -eo pipefail
source ~/mambaforge/etc/profile.d/conda.sh; conda activate assembly
W=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22; T=~/ltr_tmp/core
python3 - $T <<'PY'
import sys, random, collections
T = sys.argv[1]; random.seed(20261023)
names, seqs = [], []
for l in open(f"{T}/copies.trim.fa"):
    if l.startswith(">"): names.append(l[1:].strip()); seqs.append("")
    else: seqs[-1] += l.strip()
mac = [i for i, n in enumerate(names) if n.startswith("Mac|")]
gar = [i for i, n in enumerate(names) if n.startswith("Gar|")]
strata = collections.defaultdict(list)
for i in gar: p = names[i].split("|"); strata[(p[1], p[3])].append(i)
k = 400 / len(gar); pick = []
for s, v in strata.items(): pick += random.sample(v, max(1, round(len(v) * k)) if len(v) > 1 else 1)
sel = sorted(mac + pick[:400])
with open(f"{T}/copies.sub.fa", "w") as o:
    for i in sel: o.write(f">{names[i]}\n{seqs[i]}\n")
print(f"subsample: Mac {len(mac)} + Gar {len([i for i in sel if names[i].startswith('Gar|')])}")
PY
iqtree2 -s $T/copies.sub.fa -m GTR+F+G4 -fast -alrt 1000 -nt 2 -pre $T/copies -quiet -redo
mkdir -p $W/tree; cp $T/core_copies.fa $T/copies.aln.fa $T/copies.trim.fa $T/copies.sub.fa $T/copies.treefile $T/copies.iqtree $W/tree/
python3 $W/code/03b_tree_stats.py
