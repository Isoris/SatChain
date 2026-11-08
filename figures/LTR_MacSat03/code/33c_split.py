"""Step 33c (2026-10-24): split arrays into clean satellite blocks and spacers; TE fragment boundaries per spacer.
Array orientation: MacSat03 majority strand = +. Source coordinates: Gar dispersed RTE consensus, 1 = LINE 3' end (tail), 3,117 = 5' end.
Rules: a monomer is TE-overlapping if >= 50% of it lies in a source-RTE or insert hit; clean monomers merge into a block when gaps < 300 bp;
a spacer is the interval between two consecutive blocks."""
import pathlib, collections, pandas as pd, numpy as np
from Bio import Align
exec(open(pathlib.Path(__file__).parent / "32e_extract.py").read().split("random.seed")[0])
O = W / "split"; fa = Fa(G["Mac"])
REG = pd.read_csv(O / "regions.tsv", sep="\t")
cols = ["reg", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"]
def hits(p):
    h = pd.read_csv(p, sep="\t", header=None, names=cols)
    h["chrom"] = h.reg.str.split(":").str[0]; off = h.reg.str.split(":").str[1].astype(int)
    h["lo"] = off + h[["s0", "s1"]].min(1) - 1; h["hi"] = off + h[["s0", "s1"]].max(1); h["plus"] = h.s1 > h.s0; return h
SRC, INS = hits(O / "src_hits.tsv"), hits(O / "ins_hits.tsv")
mono = pd.read_csv("F:/data/satellite_FISH/pipeline_v4_5_fixed/out/fClaHyb_Mac_satellites/36_monomers_final.tsv", sep="\t",
                   usecols=["chrom", "genomic_start_0based", "genomic_end", "family_id", "subfamily_id", "strand", "complete_or_truncated"])
mono = mono[mono.family_id == "MacSat03"].rename(columns={"genomic_start_0based": "lo", "genomic_end": "hi"})
cons = [l.strip() for l in open("F:/data/satellite_FISH/pipeline_v4_5_fixed/out/fClaHyb_Mac_satellites/04_satellite_families.tsv")][:0]
fam = open("F:/data/satellite_FISH/pipeline_v4_5_fixed/out/fClaHyb_Mac_satellites/04_satellite_family_consensuses.fa").read().split(">")
MS = [x.split("\n", 1)[1].replace("\n", "").upper() for x in fam if x.startswith("MacSat03 ") or x.startswith("MacSat03\n")][0]
AL = Align.PairwiseAligner(mode="local", match_score=2, mismatch_score=-3, open_gap_score=-5, extend_gap_score=-2)
def phase(seq_edge, at_end):
    """MacSat03 consensus position (1..224) of the block-edge base; seq_edge in satellite orientation"""
    a = AL.align(seq_edge, MS + MS)[0]
    if a.score < 60: return np.nan
    (s0, s1), (t0, t1) = (a.aligned[0][-1], a.aligned[1][-1]) if at_end else (a.aligned[0][0], a.aligned[1][0])
    if at_end: return int((t1 - 1 - (len(seq_edge) - s1)) % len(MS)) + 1   # extend to the last base of the edge
    return int((t0 - s0) % len(MS)) + 1
VAL = O / "monomer_validation.tsv"
if VAL.exists(): V = pd.read_csv(VAL, sep="	")
else:
    mm = mono.merge(REG[["chrom", "start", "end"]], on="chrom"); mm = mm[(mm.lo >= mm.start) & (mm.hi <= mm.end)].drop(columns=["start", "end"]).drop_duplicates(["chrom", "lo", "hi"])
    fw, rv = [], []
    for r in mm.itertuples():
        x = fa.get(r.chrom, r.lo, r.hi); n = 2 * max(1, len(x)); fw.append(AL.score(x, MS + MS) / n); rv.append(AL.score(rc(x), MS + MS) / n)
    mm["score_fw"], mm["score_rc"] = fw, rv; mm["score"] = mm[["score_fw", "score_rc"]].max(1); mm["seq_strand"] = np.where(mm.score_fw >= mm.score_rc, "+", "-")
    mm.to_csv(VAL, sep="	", index=False); V = mm
print("monomer score fraction histogram (all MacSat03 calls in the arrays):", np.histogram(V.score, bins=[0, .1, .2, .3, .4, .5, .6, .7, 1])[0].tolist())
V["valid"] = V.score >= 0.30
print("validated per centromere:", V.groupby(V.chrom.str.replace("C_mac_", "")).valid.agg(["sum", "size"]).T.to_dict())
mono = V[V.valid].copy()
def cover(iv, lo, hi):
    m = np.zeros(hi - lo, bool)
    for a, b in iv: m[max(0, a - lo):max(0, min(hi, b) - lo)] = True
    return m
rows = []; removed = collections.Counter(); kept = collections.Counter()
for g in REG.itertuples():
    c, lo, hi = g.chrom, g.start, g.end; lab = c.replace("C_mac_", "")
    m = mono[(mono.chrom == c) & (mono.lo >= lo) & (mono.hi <= hi)].sort_values("lo")
    src = SRC[(SRC.chrom == c) & (SRC.lo >= lo) & (SRC.hi <= hi)]; ins = INS[(INS.chrom == c) & (INS.lo >= lo) & (INS.hi <= hi)]
    te = cover(list(zip(src.lo, src.hi)) + list(zip(ins.lo, ins.hi)), lo, hi)
    frac = np.array([te[a - lo:b - lo].mean() for a, b in zip(m.lo, m.hi)])
    removed[lab] += int((frac >= 0.5).sum()); m = m[frac < 0.5]; kept[lab] += len(m)
    # orientation from sequence: clean monomers aligned to the MacSat03 consensus in both orientations (strand labels are not reliable)
    orient = "+" if (m.seq_strand == "+").mean() >= 0.5 else "-"
    # blocks
    B = []
    for r in m.itertuples():
        if B and r.lo - B[-1][1] < 300: B[-1][1] = max(B[-1][1], r.hi); B[-1][2] += 1
        else: B.append([r.lo, r.hi, 1])
    B = [b for b in B if b[2] >= 3]
    for (a0, a1, n0), (b0, b1, n1) in zip(B[:-1], B[1:]):
        sp0, sp1 = a1, b0
        if sp1 - sp0 < 300: continue
        s = src[(src.hi > sp0) & (src.lo < sp1)]; i_ = ins[(ins.hi > sp0) & (ins.lo < sp1)]
        rec = dict(chrom=lab, region=g.Index, orient=orient, spacer_lo=sp0, spacer_hi=sp1, spacer_bp=sp1 - sp0,
                   block_before_monomers=n0, block_after_monomers=n1,
                   src_bp=int(cover(zip(s.lo, s.hi), sp0, sp1).sum()), ins_bp=int(cover(zip(i_.lo, i_.hi), sp0, sp1).sum()))
        if len(s):
            # TE fragment ends in genomic coords; source coordinate at each end
            first, last = s.loc[s.lo.idxmin()], s.loc[s.hi.idxmax()]
            q_first = first.q0 if first.plus else first.q1; q_last = last.q1 if last.plus else last.q0
            te_lo, te_hi = int(first.lo), int(last.hi)
            # left/right in array orientation
            if orient == "+": L_q, R_q, gapL, gapR = q_first, q_last, te_lo - sp0, sp1 - te_hi
            else: L_q, R_q, gapL, gapR = q_last, q_first, sp1 - te_hi, te_lo - sp0
            plus_share = (s.plus * s.len).sum() / s.len.sum()
            rec.update(te_lo=te_lo, te_hi=te_hi, te_span=te_hi - te_lo, src_left_q=int(L_q), src_right_q=int(R_q),
                       te_strand_vs_array=("+" if (plus_share >= 0.5) == (orient == "+") else "-"), gap_sat_to_te_left=gapL, gap_te_to_sat_right=gapR)
        # satellite phase at the two block edges (array orientation)
        e1 = fa.get(c, a1 - 230, a1); e2 = fa.get(c, b0, b0 + 230)
        if orient == "-": e1, e2 = rc(fa.get(c, b0, b0 + 230)), rc(fa.get(c, a1 - 230, a1))
        rec.update(phase_block_end=phase(e1, True), phase_block_start=phase(e2, False))
        rows.append(rec)
S = pd.DataFrame(rows); S.to_csv(O / "spacers.tsv", sep="\t", index=False)
print("MacSat03 monomers removed as TE-overlapping (>= 50% in TE hits):", dict(removed), "| kept:", dict(kept))
T2 = {"LG08", "LG16", "LG25"}
def mode_share(x, tol=10):
    x = x.dropna().astype(int)
    if x.empty: return ("", 0.0)
    v = x.value_counts().idxmax(); return (v, round(((x - v).abs() <= tol).mean(), 2))
out = []
for ch, d in S.groupby("chrom"):
    t = d.dropna(subset=["src_left_q"]); t = t[t.src_bp >= 300]; br = d[~d.index.isin(t.index)]
    r = dict(chrom=ch, type="II" if ch in T2 else "I", TE_spacers=len(t), sat_breaks=len(br), break_bp_median=int(br.spacer_bp.median()) if len(br) else 0,
             TE_spacer_bp=int(t.spacer_bp.median()) if len(t) else 0,
             src_bp_median=int(t.src_bp.median()) if len(t) else 0, ins_bp_median=int(d.ins_bp.median()))
    for k in ("src_left_q", "src_right_q", "gap_sat_to_te_left", "gap_te_to_sat_right", "phase_block_end", "phase_block_start"):
        v, sh = mode_share(t[k] if k in t else pd.Series(dtype=float)); r[k] = f"{v} ({sh:.0%})"
    r["te_strand"] = t.te_strand_vs_array.value_counts().idxmax() if len(t) else ""
    out.append(r)
P = pd.DataFrame(out).sort_values(["type", "chrom"]); pd.set_option("display.width", 260)
print(P.to_string(index=False)); P.to_csv(O / "spacer_summary.tsv", sep="\t", index=False)
