#!/bin/bash
# 2026-10-23: the ~6-kb type II insert. (a) domains (tblastn of reference domain consensuses, as for the RT; no HMMER), (b) copies genome-wide in Mac
# and Gar (dc-megablast of a representative insert against existing BLAST DBs), (c) where the copies lie (centromere +-1 Mb, which centromeres).
source ~/mambaforge/etc/profile.d/conda.sh; conda activate assembly
W=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22; T=~/ltr_tmp/other; mkdir -p $W/insert
# representative = longest LG25 segment
python3 - $T <<'PY'
import sys
T = sys.argv[1]; seqs = {}; n = None
for l in open(f"{T}/other.fa"):
    if l.startswith(">"): n = l[1:].strip(); seqs[n] = ""
    else: seqs[n] += l.strip()
rep = max((k for k in seqs if k.startswith("LG25")), key=lambda k: len(seqs[k]))
open(f"{T}/insert_rep.fa", "w").write(f">insert_rep_{rep}\n{seqs[rep]}\n"); print("representative", rep, len(seqs[rep]), "bp")
PY
cp $T/insert_rep.fa $W/insert/
tblastn -query $W/domains/refdomains.faa -subject $T/insert_rep.fa -evalue 1e-3 -outfmt '6 qseqid pident length sstart send evalue bitscore' | sort -k7,7gr | head -12 > $W/insert/insert_domains.tsv
echo "== protein domains in the insert (tblastn, e <= 1e-3)"; cat $W/insert/insert_domains.tsv
getorf -sequence $T/insert_rep.fa -outseq $T/insert_orfs.faa -minsize 300 -find 1 >/dev/null 2>&1; echo "ORFs >= 100 aa: $(grep -c '>' $T/insert_orfs.faa)"
for h in Gar Mac; do
  DBP=$([ $h = Gar ] && echo $W/copies/db/Gar || echo ~/ltr_tmp/Mac)
  blastn -task dc-megablast -query $T/insert_rep.fa -db $DBP -evalue 1e-20 -num_threads 2 -max_target_seqs 100 -max_hsps 20000 \
    -outfmt '6 sseqid sstart send qstart qend pident length evalue bitscore' > $W/insert/$h.insert_hits.tsv
  echo "$h insert HSPs: $(wc -l < $W/insert/$h.insert_hits.tsv)"
done
python3 - $W <<'PY'
import sys, pandas as pd, numpy as np
W = sys.argv[1]; R = "/mnt/f/data/satellite_FISH/runs_of_record_v4.5.3_2026-10-18"
for h in ("Mac", "Gar"):
    t = pd.read_csv(f"{W}/insert/{h}.insert_hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
    if t.empty: print(h, "no insert hits"); continue
    t["lo"] = t[["s0", "s1"]].min(1); t["hi"] = t[["s0", "s1"]].max(1)
    # merge hits within 2 kb into loci; locus = covered query span
    loci = []
    for c, g in t.sort_values("lo").groupby("chrom"):
        cur = None
        for r in g.itertuples():
            if cur and r.lo - cur["hi"] <= 2000: cur["hi"] = max(cur["hi"], r.hi); cur["q"].append((r.q0, r.q1)); cur["pid"].append(r.pid)
            else:
                if cur: loci.append(cur)
                cur = dict(chrom=c, lo=r.lo, hi=r.hi, q=[(r.q0, r.q1)], pid=[r.pid])
        if cur: loci.append(cur)
    L = pd.DataFrame(loci)
    def qcov(q):
        q = sorted(q); tot = 0; s = e = None
        for a, b in q:
            if e is None or a > e:
                if e is not None: tot += e - s
                s, e = a, b
            else: e = max(e, b)
        return tot + (e - s)
    L["qcov_bp"] = L.q.apply(qcov); L["pid"] = L.pid.apply(np.mean)
    cen = pd.read_csv(f"{R}/{h}/centromere_check.tsv", sep="\t").set_index("chrom")
    L["in_cen"] = [c in cen.index and lo < cen.loc[c, "cen_end"] + 1e6 and hi > cen.loc[c, "cen_start"] - 1e6 for c, lo, hi in zip(L.chrom, L.lo, L.hi)]
    big = L[L.qcov_bp >= 3000]
    print(f"== {h}: insert loci {len(L)}; with >= 3 kb of the insert: {len(big)} (in centromere: {int(big.in_cen.sum())}); median identity {big.pid.median():.1f}%")
    print("   loci >= 3 kb per chromosome:", big.chrom.str.split("_").str[-1].value_counts().to_dict())
    print("   small loci (< 3 kb) in/out of centromeres:", int(L[L.qcov_bp < 3000].in_cen.sum()), "/", int((~L[L.qcov_bp < 3000].in_cen).sum()))
    L.drop(columns=["q"]).to_csv(f"{W}/insert/{h}.insert_loci.tsv", sep="\t", index=False)
PY
