"""Step 34b (2026-10-24): refine the RTE fragment in every SatChain v4.6 interruption between two MacSat03 (CmaSat001-225) runs.
Spacer = v4.6 interruption whose left and right runs are both MacSat03. Array orientation = orientation of the flanking runs
(catalogue reverse_fraction < 0.5 -> +). Fragment = local alignment of the spacer to the source RTE (Gar dispersed consensus, 3,117 bp;
1 = LINE 3' end), exact ends and internal deletions (alignment gaps >= 50 bp in the spacer relative to the source)."""
import pathlib, pandas as pd, numpy as np, collections
from Bio import Align
exec(open(pathlib.Path(__file__).parent / "32e_extract.py").read().split("random.seed")[0])
D = pathlib.Path("F:/data/fClaHyb_func_annotation_v1.0/fClaHyb_Mac/fClaHyb_Mac_satelliteDNA"); O = W / "split"; fa = Fa(G["Mac"])
C = pd.read_csv(D / "fClaHyb_Mac.satDNA_v46.catalogue.tsv", sep="\t").set_index("run_id")
I = pd.read_csv(D / "fClaHyb_Mac.satDNA_v46.interruptions_between_runs.tsv", sep="\t")
I = I[I.run_left.map(C.family).eq("MacSat03") & I.run_right.map(C.family).eq("MacSat03")].copy()
SRC = "".join(l.strip() for l in open(W / "seq" / "Gar_full_consensus.fa") if not l.startswith(">"))
AL = Align.PairwiseAligner(mode="local", match_score=2, mismatch_score=-3, open_gap_score=-6, extend_gap_score=-1)
T2 = {"LG08", "LG16", "LG25"}
rows, frags = [], collections.defaultdict(list)
for r in I.itertuples():
    rev = (C.reverse_fraction[r.run_left] + C.reverse_fraction[r.run_right]) / 2 >= 0.5
    s = fa.get(r.chrom, r.start0, r.end); s = rc(s) if rev else s
    # iterative local alignment: best segment, then search the unaligned parts of the spacer again (both strands); keep segments with
    # >= 100 aligned bp and score >= 150; the fragment is the set of kept segments on the dominant strand
    segs = []; todo = [(0, len(s))]
    for _ in range(6):
        nxt = []
        for x0, x1 in todo:
            if x1 - x0 < 100: continue
            sub = s[x0:x1]; best = None
            for strand, q in (("+", SRC), ("-", rc(SRC))):
                al = AL.align(sub, q)[0]
                if best is None or al.score > best[0].score: best = (al, strand)
            al, strand = best
            n_al = sum(e - b_ for b_, e in al.aligned[0])
            if al.score < 150 or n_al < 100: continue
            u0, u1 = al.aligned[0][0][0] + x0, al.aligned[0][-1][1] + x0
            v0, v1 = al.aligned[1][0][0], al.aligned[1][-1][1]
            if strand == "-": v0, v1 = len(SRC) - v1, len(SRC) - v0
            ids = sum(1 for (b_, e), (c_, d) in zip(al.aligned[0], al.aligned[1]) for i in range(e - b_) if sub[b_ + i] == (SRC if strand == "+" else rc(SRC))[c_ + i])
            segs.append(dict(sp0=u0, sp1=u1, q0=v0, q1=v1, strand=strand, n=n_al, pid=ids / n_al))
            nxt += [(x0, u0), (u1, x1)]
        todo = nxt
    if not segs:
        ch = r.chrom.replace("C_mac_", "")
        rows.append(dict(chrom=ch, type="II" if ch in T2 else "I", start0=r.start0, end=r.end, spacer_bp=r.bp, rev=rev, te_strand="", frag_bp=0, aligned_bp=0, deletions="", segments=""))
        continue
    dom = max(("+", "-"), key=lambda t: sum(g["n"] for g in segs if g["strand"] == t)); segs = sorted([g for g in segs if g["strand"] == dom], key=lambda g: g["sp0"])
    strand = dom; sp0, sp1 = segs[0]["sp0"], segs[-1]["sp1"]
    q_first, q_last = (segs[0]["q0"], segs[-1]["q1"]) if dom == "+" else (segs[0]["q1"], segs[-1]["q0"])
    qs = sorted((g["q0"], g["q1"]) for g in segs); dels = [(a1, b0) for (a0, a1), (b0, b1) in zip(qs[:-1], qs[1:]) if b0 - a1 >= 50]
    ident = sum(g["n"] for g in segs); pid = sum(g["n"] * g["pid"] for g in segs) / ident
    ch = r.chrom.replace("C_mac_", "")
    rows.append(dict(chrom=ch, type="II" if ch in T2 else "I", start0=r.start0, end=r.end, spacer_bp=r.bp, rev=rev, te_strand=strand,
                     frag_spacer_from=sp0, frag_spacer_to=sp1, frag_bp=sp1 - sp0, src_end_left=int(q_first), src_end_right=int(q_last),
                     aligned_bp=ident, pid=round(pid * 100, 1), segments=";".join(f"{g['q0']}-{g['q1']}" for g in segs), deletions=";".join(f"{x}-{y}" for x, y in dels), left_flank_bp=sp0, right_flank_bp=len(s) - sp1,
                     junction_left=s[max(0, sp0 - 25):sp0 + 25], junction_right=s[max(0, sp1 - 25):sp1 + 25]))
    frags[ch].append(s[sp0:sp1])
R = pd.DataFrame(rows); R.to_csv(O / "v46_fragments.tsv", sep="\t", index=False)
with open(O / "v46_fragments.fa", "w") as fh:
    for ch, L in frags.items():
        for i, x in enumerate(L): fh.write(f">{ch}_{i}\n{x}\n")
def mode_share(x, tol=10):
    x = pd.Series(x).dropna()
    if x.empty: return ""
    v = x.value_counts().idxmax()
    try: sh = ((x.astype(float) - float(v)).abs() <= tol).mean()
    except ValueError: sh = (x == v).mean()
    return f"{v} ({sh:.0%})"
out = []
R0 = R; R = R[R.aligned_bp >= 100]
print("interruptions without an RTE fragment (>= 100 aligned bp):", R0[R0.aligned_bp < 100].groupby("chrom").size().to_dict())
for ch, d in R.groupby("chrom"):
    out.append(dict(chrom=ch, type=d.type.iloc[0], spacers=len(d), spacer_bp=int(d.spacer_bp.median()), te_strand=mode_share(d.te_strand),
                    frag_bp=int(d.frag_bp.median()), pid=d.pid.median(), share_of_spacers=f"{len(d) / len(R0[R0.chrom == ch]):.0%}", src_left=mode_share(d.src_end_left), src_right=mode_share(d.src_end_right),
                    deletions=mode_share(d.deletions.replace("", "none")), left_flank=int(d.left_flank_bp.median()), right_flank=int(d.right_flank_bp.median()),
                    junction_left_shared=f"{d.junction_left.value_counts().iloc[0] / len(d):.0%}", junction_right_shared=f"{d.junction_right.value_counts().iloc[0] / len(d):.0%}"))
P = pd.DataFrame(out).sort_values(["type", "chrom"]); pd.set_option("display.width", 260); print(P.to_string(index=False))
P.to_csv(O / "v46_fragment_summary.tsv", sep="\t", index=False)
