#!/bin/bash
# Step 7 (2026-10-22, user: check its expression): is the lineage A-associated RTE-like element transcribed?
# Per species, a small reference panel: the element consensus of that species (majority of its near-full copies in tree/copies.aln.fa),
# lineage A and B satellite consensuses, and housekeeping transcripts (EF1a LOC128518327, gapdh, rpl13a; LiftOn, same species).
# First 2 M spots of every public run (rnaseq_public_v1) -> minimap2 -ax sr -> primary hits with >= 50 aligned bp and <= 10% mismatches.
# Output: expression/rpkm.tsv (reads per kb per million read ends), expression/<run>.element_depth.tsv (coverage along the element).
# Libraries are mostly unstranded: the direction of transcription is not tested.
set -eo pipefail
source ~/mambaforge/etc/profile.d/conda.sh; conda activate assembly
W=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22; R=/mnt/f/data/rnaseq_public_v1; A=/mnt/f/data/fClaHyb_func_annotation_v1.0
S=/mnt/f/data/satellite_FISH/pipeline_v4_5_fixed/out; T=~/ltr_tmp/expr; mkdir -p $T $W/expression
python3 - $W $A $S $T <<'PY'
import sys, collections
W, A, S, T = sys.argv[1:5]
def fasta(p):
    n, s = None, []
    for l in open(p):
        if l.startswith(">"):
            if n: yield n, "".join(s)
            n, s = l[1:].strip(), []
        else: s.append(l.strip())
    if n: yield n, "".join(s)
aln = list(fasta(f"{W}/tree/copies.aln.fa"))
SAT = {"Gar": ("GarSat03", "GarSat04"), "Mac": ("MacSat03", "MacSat06")}
for h in ("Gar", "Mac"):
    seqs = [s for n, s in aln if n.startswith(h + "|")]
    cons = "".join(c for c in (collections.Counter(col).most_common(1)[0][0] for col in zip(*seqs)) if c != "-").upper()
    fam = dict((n.split()[0], s) for n, s in fasta(f"{S}/fClaHyb_{h}_satellites/04_satellite_family_consensuses.fa"))
    tx = {}
    for n, s in fasta(f"{A}/fClaHyb_{h}/_genes_lifton/sequences/fClaHyb_{h}.lifton.all.transcripts.fna"):
        g = n.split("|")[0]
        if g in ("LOC128518327", "gapdh", "rpl13a") and g not in tx: tx[g] = s
    with open(f"{T}/{h}_ref.fa", "w") as o:
        o.write(f">element_{h}\n{cons}\n")
        for f in SAT[h]:
            k = next((k for k in fam if k.split("#")[0] == f or k.startswith(f + "_") or k == f), None)
            if k: o.write(f">sat_{f}\n{fam[k]}\n")
        for g, s in tx.items(): o.write(f">hk_{'EF1a' if g == 'LOC128518327' else g}\n{s}\n")
    open(f"{W}/seq/element_{h}_consensus.fa", "w").write(f">element_{h}_majority_consensus\n{cons}\n")
    print(h, "element consensus", len(cons), "bp; housekeeping:", ",".join(tx), file=sys.stderr)
PY
echo -e "species\trun\treads\tref\tref_bp\tcount\trpkm" > $W/expression/rpkm.tsv
for h in Gar Mac; do
  for sra in $R/$h/*.sra; do
    run=$(basename $sra .sra)
    fastq-dump -X 2000000 --split-spot --stdout $sra 2>/dev/null | tee >(awk 'NR % 4 == 1' | wc -l > $T/$run.n) \
      | minimap2 -ax sr -t 2 --secondary=no $T/${h}_ref.fa - 2>/dev/null | samtools view -F 2308 - > $T/$run.sam
    python3 - $T/$run.sam $T/${h}_ref.fa $(cat $T/$run.n) $h $run $W <<'PY'
import sys, re, collections
sam, ref, n, h, run, W = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6]
L, cur = {}, None
for l in open(ref):
    if l.startswith(">"): cur = l[1:].strip(); L[cur] = 0
    else: L[cur] += len(l.strip())
cnt = collections.Counter(); depth = collections.Counter()
for l in open(sam):
    p = l.split("\t"); cig = re.findall(r"(\d+)([MIDNSHP=X])", p[5])
    aln = sum(int(k) for k, op in cig if op in "M=X"); nm = next((int(x[5:]) for x in p[11:] if x.startswith("NM:i:")), 99)
    if aln < 50 or nm > 0.1 * aln: continue
    cnt[p[2]] += 1
    if p[2].startswith("element"):
        pos = int(p[3])
        for b in range(pos // 100, (pos + aln) // 100 + 1): depth[b] += 1
with open(f"{W}/expression/rpkm.tsv", "a") as o:
    for r, bp in L.items(): o.write(f"{h}\t{run}\t{n}\t{r}\t{bp}\t{cnt[r]}\t{cnt[r] / (bp / 1e3) / (n / 1e6):.2f}\n")
with open(f"{W}/expression/{run}.element_depth.tsv", "w") as o:
    o.write("bin_start\treads\n"); [o.write(f"{b * 100}\t{depth[b]}\n") for b in range(L[f'element_{h}'] // 100 + 1)]
print(run, n, "read ends;", dict(cnt))
PY
    rm -f $T/$run.sam
  done
done
