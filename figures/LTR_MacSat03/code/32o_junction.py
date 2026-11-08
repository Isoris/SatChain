"""Which element end abuts the MacSat03 block in each unit, and the sequence at the 3'-tail junction."""
import pathlib, pandas as pd, numpy as np, collections
W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22")
U = pd.read_csv(W / "unit_anchored" / "units_A6500.tsv", sep="\t")
mono = pd.read_csv("F:/data/satellite_FISH/pipeline_v4_5_fixed/out/fClaHyb_Mac_satellites/36_monomers_final.tsv", sep="\t", usecols=["chrom", "genomic_start_0based", "genomic_end", "family_id"])
mono = mono[mono.family_id == "MacSat03"]
B = pd.read_csv(W / "homologous" / "src_vs_units.tsv", sep="\t", header=None, names=["unit", "q0", "q1", "s0", "s1", "pid", "len", "bits"])
seqs = {}; k = None
for l in open(W / "seq" / "units_all.fa"):
    l = l.strip()
    if l.startswith(">"): k = l[1:]; seqs[k] = []
    else: seqs[k].append(l)
seqs = {k: "".join(v) for k, v in seqs.items()}
rows = []; tails = collections.Counter()
for r in U.itertuples():
    u = f"{r.chrom.replace('C_mac_', '')}|{r.lo}"; d = B[B.unit == u]
    if d.empty: continue
    t = d.loc[d.q0.idxmin()]  # hit holding the 3' end (lowest source coordinate)
    if t.q0 > 30: continue
    # unit coordinate of the 3' tail: source runs low->high on the + strand of the unit if s1 > s0
    plus = t.s1 > t.s0; tail_pos = t.s0  # tail is at query start
    L = r.hi - r.lo
    m = mono[(mono.chrom == r.chrom) & (mono.genomic_start_0based >= r.lo) & (mono.genomic_end <= r.hi)]
    ms = (m.genomic_start_0based - r.lo).values if r.strand == "+" else (r.hi - m.genomic_end).values
    me = (m.genomic_end - r.lo).values if r.strand == "+" else (r.hi - m.genomic_start_0based).values
    # distance from the tail to the nearest monomer on the side away from the element body
    if plus: side = ms[ms < tail_pos]; dist = tail_pos - me[ms < tail_pos].max() if len(side) else np.nan
    else: side = ms[ms > tail_pos]; dist = side.min() - tail_pos if len(side) else np.nan
    s = seqs[u]; j = s[max(0, tail_pos - 60):tail_pos] if plus else s[tail_pos:tail_pos + 60]
    j = j if not plus else j.translate(str.maketrans("ACGT", "TGCA"))[::-1]  # show in LINE sense: element 3' end then flank
    rows.append(dict(unit=u, chrom=u.split("|")[0], tail_to_satellite_bp=dist, tail_seq=j))
R = pd.DataFrame(rows)
print(R.groupby("chrom").tail_to_satellite_bp.describe()[["count", "50%", "min", "max"]].round(0).to_string())
print("\nexample 3'-tail junction sequences (LINE sense, element end -> satellite side), one per chromosome:")
for c, d in R.groupby("chrom"): print(f"  {c}: {d.tail_seq.iloc[len(d)//2]}")
R.to_csv(W / "homologous" / "tail_junctions.tsv", sep="\t", index=False)
