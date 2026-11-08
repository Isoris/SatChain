#!/bin/bash
# 2026-10-23: satellite-derived segments inside the element consensus. Sensitive blastn (word 11) of the consensus against the lineage A/B satellite
# family consensuses of both haplotypes, plus a self-comparison for internal tandem repeats. Output: domains/element_satellite_segments.tsv
source ~/mambaforge/etc/profile.d/conda.sh; conda activate assembly
W=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22; S=/mnt/f/data/satellite_FISH/pipeline_v4_5_fixed/out
python3 - "$W" "$S" <<'PY'
import sys
W, S = sys.argv[1], sys.argv[2]
FAMS = {"Gar": ("GarSat03", "GarSat04", "GarSat09", "GarSat10", "GarSat49"), "Mac": ("MacSat03", "MacSat06", "MacSat51", "MacSat162")}
with open(W + "/seq/satellite_refs.fa", "w") as out:
    for h, fams in FAMS.items():
        keep = False
        for l in open(f"{S}/fClaHyb_{h}_satellites/04_satellite_family_consensuses.fa"):
            if l.startswith(">"): keep = l[1:].split()[0].split("#")[0] in fams
            if keep: out.write(l)
PY
echo "satellite references: $(grep -c '>' $W/seq/satellite_refs.fa)"
blastn -task blastn -word_size 11 -query $W/seq/TE_00003409_library.fa -subject $W/seq/satellite_refs.fa -evalue 1e-5 \
  -outfmt '6 sseqid pident length qstart qend evalue bitscore' > $W/domains/element_vs_lineage_satellites_sensitive.tsv
echo "== vs satellite families (sorted by element position)"; sort -k4,4n $W/domains/element_vs_lineage_satellites_sensitive.tsv | awk '{print $1, $2"%", $4"-"$5, "e="$6}' | uniq | head -30
blastn -task blastn -word_size 11 -query $W/seq/TE_00003409_library.fa -subject $W/seq/TE_00003409_library.fa -evalue 1e-10 \
  -outfmt '6 qstart qend sstart send pident length' | awk '!($1==$3 && $2==$4)' > $W/domains/element_self.tsv
echo "== internal repeats (self hits off the diagonal): $(wc -l < $W/domains/element_self.tsv)"; sort -k1,1n $W/domains/element_self.tsv | head -12
