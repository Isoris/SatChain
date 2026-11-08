#!/bin/bash
# 2026-10-23: what is the ~7.7 kb "other" sequence in type II units (Mac LG08, LG25, LG16)? Longest stretch per unit not covered by satellite runs
# or element hits; up to 6 units per chromosome; blastn vs the Mac EDTA TE library and vs each other (one insert type or several).
source ~/mambaforge/etc/profile.d/conda.sh; conda activate assembly
W=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22; G=/mnt/f/data/fClaHyb_func_annotation_v1.0/fClaHyb_Mac; T=~/ltr_tmp/other; mkdir -p $T $W/unit
python3 - $W > $T/regions.tsv <<'PY'
import sys, pandas as pd
W = sys.argv[1]
u = pd.read_csv(f"{W}/unit/units.tsv", sep="\t"); u = u[u.chrom.isin(["C_mac_LG08", "C_mac_LG25", "C_mac_LG16"]) & (u.other_bp > 5000)]
runs = pd.read_csv("/mnt/f/data/satellite_FISH/runs_of_record_v4.5.3_2026-10-18/Mac/runs_of_record.tsv", sep="\t")
runs = runs[(runs.unit_class == "satellite") & runs.tier.isin(["S1", "S2"])]
h = pd.read_csv(f"{W}/copies/Mac.hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
h = h[~((h.q0 >= 2900) & (h.q1 <= 3420))]; h["lo"] = h[["s0", "s1"]].min(1); h["hi"] = h[["s0", "s1"]].max(1)
for ch, g in u.groupby("chrom"):
    for r in g.head(6).itertuples():
        a, b = r.start, r.start + r.period
        cov = sorted([(x, y) for x, y in zip(runs[runs.chrom == ch].start0, runs[runs.chrom == ch].end) if y > a and x < b] +
                     [(x, y) for x, y in zip(h[h.chrom == ch].lo, h[h.chrom == ch].hi) if y > a and x < b])
        best, cur = (0, a, a), a
        for x, y in cov:
            if x > cur and x - cur > best[0]: best = (x - cur, cur, x)
            cur = max(cur, y)
        if b - cur > best[0]: best = (b - cur, cur, b)
        if best[0] > 1000: print(f"{ch}:{best[1] + 1}-{best[2]}\t{ch.split('_')[-1]}_{best[1]}")
PY
: > $T/other.fa
while IFS=$'\t' read -r reg name; do samtools faidx $G/_genome/fClaHyb_Mac_LG.fa "$reg" | sed "1s/.*/>$name/" >> $T/other.fa; done < $T/regions.tsv
echo "segments: $(grep -c '>' $T/other.fa); lengths: $(awk '/^>/{if(l)print l; l=0; next}{l+=length($0)} END{print l}' $T/other.fa | sort -n | tr '\n' ' ')"
blastn -task dc-megablast -query $T/other.fa -subject $G/_repeats_EDTA/fClaHyb_Mac.TElib.fa -evalue 1e-10 -outfmt '6 qseqid sseqid pident length qstart qend evalue bitscore' > $T/other_vs_TElib.tsv
echo "== best TE-library hits per segment (bp covered by each library entry)"
awk -F'\t' '{c[$1"\t"$2]+=$4} END{for(k in c) print c[k]"\t"k}' $T/other_vs_TElib.tsv | sort -k2,2 -k1,1nr | awk -F'\t' '{if(n[$2]++<3) print}'
blastn -task dc-megablast -query $T/other.fa -subject $T/other.fa -evalue 1e-20 -outfmt '6 qseqid sseqid pident length' | awk '$1!=$2' | sort -k4,4nr | awk '!s[$1]++' | head -20 > $T/other_self.tsv
echo "== segment vs segment (best partner)"; cat $T/other_self.tsv
cp $T/other.fa $T/other_vs_TElib.tsv $T/other_self.tsv $W/unit/
