#!/bin/bash
# 2026-10-23 (user: use k-mers): is the RTE-like element transcribed? Low-memory k-mer test (jellyfish --if counts only listed k-mers).
# Element k-mers: canonical 31-mers of each species' own core copies (core 5,560-7,400) present in >= 3 copies and with every genome occurrence inside a
# core copy (genome count == count in core copies). Housekeeping k-mers: 31-mers of EF1a (LOC128518327), gapdh, rpl13a transcripts with genome count 1.
# Reads: first 2 M spots of every public run. Output: expression/kmer_depth.tsv (mean k-mer count per million read ends, per target).
set -eo pipefail
source ~/mambaforge/etc/profile.d/conda.sh
W=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22; R=/mnt/f/data/rnaseq_public_v1; A=/mnt/f/data/fClaHyb_func_annotation_v1.0
T=~/ltr_tmp/kexpr; mkdir -p $T $W/expression
JF=~/mambaforge/envs/kmer/bin/jellyfish
conda activate assembly
for h in Gar Mac; do
  python3 - $h $W $A $T <<'PY'
import sys, collections
h, W, A, T = sys.argv[1:5]; K = 31
comp = str.maketrans("ACGT", "TGCA")
def canon(k): r = k.translate(comp)[::-1]; return min(k, r)
def fasta(p):
    n, s = None, []
    for l in open(p):
        if l.startswith(">"):
            if n: yield n, "".join(s).upper()
            n, s = l[1:].strip(), []
        else: s.append(l.strip())
    if n: yield n, "".join(s).upper()
occ = collections.Counter()
for n, s in fasta(f"{W}/tree/core_copies.fa"):
    if not n.startswith(h + "|"): continue
    seen = {canon(s[i:i + K]) for i in range(len(s) - K + 1) if set(s[i:i + K]) <= set("ACGT")}
    occ.update(seen)
elem = [k for k, c in occ.items() if c >= 3]
with open(f"{T}/{h}_elem_cand.fa", "w") as o:
    for i, k in enumerate(elem): o.write(f">e{i}\n{k}\n")
hk = {}
for n, s in fasta(f"{A}/fClaHyb_{h}/_genes_lifton/sequences/fClaHyb_{h}.lifton.all.transcripts.fna"):
    g = n.split("|")[0]
    if g in ("LOC128518327", "gapdh", "rpl13a") and g not in hk: hk[g] = s
with open(f"{T}/{h}_hk_cand.fa", "w") as o:
    for g, s in hk.items():
        ks = {canon(s[i:i + K]) for i in range(len(s) - K + 1) if set(s[i:i + K]) <= set("ACGT")}
        for i, k in enumerate(sorted(ks)): o.write(f">{g}_{i}\n{k}\n")
print(h, "element candidate k-mers (>= 3 copies):", len(elem), "| housekeeping genes:", ",".join(hk), file=sys.stderr)
PY
  cat $T/${h}_elem_cand.fa $T/${h}_hk_cand.fa > $T/${h}_cand.fa
  $JF count -m 31 -C -s 20M -t 2 --if $T/${h}_cand.fa -o $T/${h}_genome.jf $A/fClaHyb_$h/_genome/fClaHyb_${h}_LG.fa
  python3 -c "
import sys
s=open('$W/tree/core_copies.fa').read().split('>')[1:]
open('$T/${h}_copies.fa','w').write(''.join('>'+x for x in s if x.startswith('$h|')))"
  $JF count -m 31 -C -s 20M -t 2 --if $T/${h}_cand.fa -o $T/${h}_copies.jf $T/${h}_copies.fa
  $JF dump -c $T/${h}_genome.jf > $T/${h}_genome.counts; $JF dump -c $T/${h}_copies.jf > $T/${h}_copies.counts
done
python3 - $T <<'PY'
import sys
T = sys.argv[1]
for h in ("Gar", "Mac"):
    g = dict(l.split() for l in open(f"{T}/{h}_genome.counts")); c = dict(l.split() for l in open(f"{T}/{h}_copies.counts"))
    lab = {}
    for f in (f"{T}/{h}_elem_cand.fa", f"{T}/{h}_hk_cand.fa"):
        name = None
        for l in open(f):
            if l.startswith(">"): name = l[1:].strip()
            else: lab[l.strip()] = "element" if name.startswith("e") else name.rsplit("_", 1)[0]
    keep = {}
    for k, t in lab.items():
        gc = int(g.get(k, 0)); cc = int(c.get(k, 0))
        if t == "element" and gc > 0 and gc == cc: keep[k] = t
        elif t != "element" and gc == 1: keep[k] = t
    with open(f"{T}/{h}_selected.fa", "w") as o:
        for i, (k, t) in enumerate(keep.items()): o.write(f">{t}|{i}\n{k}\n")
    from collections import Counter
    print(h, "selected k-mers:", dict(Counter(keep.values())), file=sys.stderr)
PY
echo -e "species\trun\tread_ends\ttarget\tn_kmers\tkmers_seen\tmean_count\tmean_per_million" > $W/expression/kmer_depth.tsv
for h in Gar Mac; do
  for sra in $R/$h/*.sra; do
    run=$(basename $sra .sra)
    fastq-dump -X 2000000 --split-spot --stdout $sra 2>/dev/null > $T/$run.fq
    n=$(( $(wc -l < $T/$run.fq) / 4 ))
    $JF count -m 31 -C -s 20M -t 2 --if $T/${h}_selected.fa -o $T/$run.jf $T/$run.fq; rm -f $T/$run.fq
    $JF dump -c $T/$run.jf > $T/$run.counts
    python3 - $T $h $run $n $W <<'PY'
import sys, collections
T, h, run, n, W = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
lab = {}; name = None
for l in open(f"{T}/{h}_selected.fa"):
    if l.startswith(">"): name = l[1:].split("|")[0]
    else: lab[l.strip()] = name
cnt = dict((k, int(v)) for k, v in (l.split() for l in open(f"{T}/{run}.counts")))
tot = collections.defaultdict(list)
for k, t in lab.items(): tot[t].append(cnt.get(k, 0))
with open(f"{W}/expression/kmer_depth.tsv", "a") as o:
    for t, v in sorted(tot.items()):
        m = sum(v) / len(v); o.write(f"{h}\t{run}\t{n}\t{t}\t{len(v)}\t{sum(1 for x in v if x)}\t{m:.3f}\t{m / (n / 1e6):.3f}\n")
print(run, n, {t: round(sum(v) / len(v) / (n / 1e6), 3) for t, v in tot.items()})
PY
    rm -f $T/$run.jf
  done
done
