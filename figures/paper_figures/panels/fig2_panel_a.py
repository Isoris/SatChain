"""Figure 2, panel a (house style module): the RTE LINE, 3,117-bp consensus of 37 full-length C. gariepinus copies (88 x 52 mm).
   draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None) -> dict of the numbers shown.
   Standalone:  python fig2_panel_a.py  -> panels/_standalone/fig2_panel_a.{png,pdf,svg}, fig2_panel_a_numbers.tsv
All numbers are recomputed here from the repo tables (LTR_MacSat03_2026-10-22): seq/Gar_full_consensus.fa, seq/Gar_full_orfs_minus.faa, seq/Gar_full.aln.fa (37 copies),
domains/refdomains.faa (RT: local alignment, BLOSUM62 11/1, no BLAST/HMMER needed). Coordinates: LINE sense, 5' -> 3' left to right; LINE position p = 3,118 - q, q = consensus coordinate
(the consensus file is written in satellite-array orientation, q = 1 = LINE 3' end)."""
import sys, re, pathlib, collections, functools
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "code"))
import satellite_style as ST
from matplotlib.patches import Rectangle
from matplotlib.colors import to_rgb

W = pathlib.Path(r"F:\data\satellite_FISH\LTR_MacSat03_2026-10-22")
RT_DARK = "#8a6208"            # RT region: dark shade of the RTE gold (no black, no brown block)
UTR_GREY = "#a7aeb3"


def _fa(p):
    d, k = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"): k = l[1:]; d[k] = []
        elif k: d[k].append(l.upper())
    return {k: "".join(v) for k, v in d.items()}


@functools.lru_cache(maxsize=1)
def data():
    from Bio import Align, SeqIO
    from Bio.Align import substitution_matrices
    from Bio.Seq import Seq
    cons = next(iter(_fa(W / "seq" / "Gar_full_consensus.fa").values())); L = len(cons); line = str(Seq(cons).reverse_complement())
    orf = [l.strip() for l in open(W / "seq" / "Gar_full_orfs_minus.faa")][1]
    assert L == 3117 and len(orf) == 1027 and str(Seq(line[:3 * len(orf)]).translate()) == orf
    orf_end = 3 * len(orf)                                              # LINE 1-3,081
    ape = re.search(r"[CVIL][LVI]QET", orf); first_met = orf.find("M") + 1
    # --- RT region and best reference by local alignment (BLOSUM62, gap 11/1) of the ORF against every LINE reference RT
    refs = {r.id: str(r.seq) for r in SeqIO.parse(W / "domains" / "refdomains.faa", "fasta")}; tot = sum(map(len, refs.values()))
    al = Align.PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"), open_gap_score=-11, extend_gap_score=-1)
    hits = []
    for k, v in refs.items():
        if "LINE" not in k: continue
        a = al.align(orf, v)[0]; n = sum(e - b for b, e in a.aligned[0])
        idn = sum(1 for (b, e), (c, d) in zip(*a.aligned) for i in range(e - b) if orf[b + i] == v[c + i])
        hits.append(dict(ref=k.split("/")[-1].replace(":LINE-RT", ""), score=a.score, q0=a.aligned[0][0][0] + 1, q1=a.aligned[0][-1][1], pid=100 * idn / n, aln=n,
                         E=0.041 * len(orf) * tot * np.exp(-0.267 * a.score)))
    H = pd.DataFrame(hits).sort_values("score", ascending=False).reset_index(drop=True)
    # --- 37 copies: per-position conservation and direct repeats (TSD-like) at the copy ends
    A = _fa(W / "seq" / "Gar_full.aln.fa"); M = np.array([list(s) for s in A.values()]); n, Lc = M.shape; occ = (M != "-").mean(0)
    ident, base = np.zeros(Lc), []
    for j in range(Lc):
        c = collections.Counter(x for x in M[:, j] if x != "-"); b, f = c.most_common(1)[0] if c else ("-", 0); ident[j] = f / n; base.append(b)
    w = np.convolve(ident, np.ones(50) / 50, "same"); body = np.where(w >= 0.6)[0]; c0, c1 = body[0], body[-1]
    kept = [j for j in range(c0, c1 + 1) if occ[j] >= 0.5]; assert len(kept) == L
    assert "".join(base[j] for j in kept) == cons
    carry = np.array([(M[:, j] == base[j]).mean() for j in kept])           # fraction of the 37 copies with the consensus base, per q
    cons_line = carry[::-1]                                                 # q -> LINE position
    tsd = []
    for k, s in A.items():
        left = "".join(x for x in s[:c0] if x != "-")[-25:]; right = "".join(x for x in s[c1 + 1:] if x != "-")[:25]; best = ""
        for ln in range(25, 4, -1):
            for i in range(0, 25 - ln + 1):
                sub = left[i:i + ln]
                if len(sub) == ln and sub in right[:ln + 10] and i >= 25 - ln - 10: best = sub; break
            if best: break
        tsd.append(len(best))
    tsd = np.array(tsd)
    # --- 5' end: where the 37 copies begin relative to the alignment body (in LINE bp)
    return dict(L=L, line=line, orf=orf, orf_end=orf_end, ape_aa=ape.start() + 1, ape_motif=ape.group(), first_met=first_met, H=H, n=n, cons_line=cons_line,
                tsd_n={m: int((tsd >= m).sum()) for m in (5, 6, 8, 10)}, tsd_len=sorted(int(x) for x in tsd[tsd >= 5]), tail=line[-13:], utr3=L - orf_end)


def draw(fig, left_mm, top_mm, w_mm=88.0, h_mm=52.0, letter=None):
    D = data(); L = D["L"]; H = D["H"]
    if letter: ST.letter(fig, left_mm, top_mm, letter)
    ov = ST.ax_mm(fig, left_mm, top_mm, w_mm, h_mm); ov.set_xlim(0, w_mm); ov.set_ylim(h_mm, 0); ov.axis("off"); ov.patch.set_alpha(0)
    X0, X1 = 6.0, w_mm - 4.0; sc = (X1 - X0) / L; X = lambda p: X0 + p * sc                          # p = LINE position (bp) -> mm
    ov.text(4.2, 0.2, "RTE LINE: consensus of 37 full-length copies", fontsize=7, fontweight="bold", ha="left", va="top")
    # --- element bar
    BT, BH = 13.0, 7.0
    aa = lambda a: 3 * a - 2
    rt0, rt1 = int(H.q0[0]), int(H.q1[0]); xrt0, xrt1 = X(aa(rt0) - 1), X(3 * rt1)
    ov.add_patch(Rectangle((X(0), BT), X(D["orf_end"]) - X(0), BH, fc=ST.RTE, ec="none", alpha=0.55, lw=0))                       # ORF
    ov.add_patch(Rectangle((xrt0, BT), xrt1 - xrt0, BH, fc=RT_DARK, ec="none"))                                                  # RT
    ov.add_patch(Rectangle((X(D["orf_end"]), BT + 1.6), X(L) - X(D["orf_end"]), BH - 3.2, fc=UTR_GREY, ec="none"))                # 3' end after the stop
    ov.text((X(0) + xrt0) / 2, BT + BH / 2, "N-terminal region", fontsize=7, ha="center", va="center", color=ST.NAVY)
    ov.text((xrt0 + xrt1) / 2, BT + BH / 2, "reverse\ntranscriptase", fontsize=7, ha="center", va="center", color="white", linespacing=0.95)
    ov.text(X(0) - 1.0, BT + BH / 2, "5\u2032", fontsize=7, ha="right", va="center"); ov.text(X(L) + 1.0, BT + BH / 2, "3\u2032", fontsize=7, ha="left", va="center")
    # --- APE-like motif pin (above the bar)
    xa = X(aa(D["ape_aa"]) + 1)
    ov.plot([xa, xa], [BT - 2.0, BT + BH], color=ST.NAVY, lw=0.5, solid_capstyle="butt")
    ov.text(xa - 0.4, BT - 5.6, f"APE-like motif {D['ape_motif']} (aa {D['ape_aa']})", fontsize=7, ha="left", va="bottom")
    # --- 3' end label above the right end
    ov.text(X(L) + 2.0, BT - 5.6, "(AAG)n", fontsize=7, ha="right", va="bottom")
    ov.plot([X(L) - 0.5] * 2, [BT - 2.0, BT + 1.6], color=ST.NAVY, lw=0.5, solid_capstyle="butt")
    # --- ORF bracket + label
    by = BT + BH + 1.6
    ov.plot([X(0), X(D["orf_end"])], [by, by], color=ST.NAVY, lw=0.5, solid_capstyle="butt")
    for p in (0, D["orf_end"]): ov.plot([X(p)] * 2, [by - 0.7, by + 0.7], color=ST.NAVY, lw=0.5, solid_capstyle="butt")
    ov.text((X(0) + X(D["orf_end"])) / 2, by + 0.9, f"one ORF, {len(D['orf']):,} aa (LINE 1\u2013{D['orf_end']:,}), open at the 5\u2032 end", fontsize=7, ha="center", va="top")
    # --- scale
    SY = by + 5.4
    ov.plot([X(0), X(L)], [SY] * 2, color=ST.NAVY, lw=0.5)
    for p in range(0, L + 1, 500):
        ov.plot([X(p)] * 2, [SY, SY + 1.0], color=ST.NAVY, lw=0.5, solid_capstyle="butt")
        ov.text(X(p), SY + 1.3, f"{p:,}", fontsize=7, ha="center", va="top")
    ov.text(X(L), SY + 1.3 + 3.0, "LINE position (bp)", fontsize=7, ha="right", va="top")
    # --- notes
    t2 = D["tsd_n"][5]; top = H.iloc[0]
    # --- Pfam domains (hmmscan of the six-frame ORFs against the fetched Pfam HMMs, 2026-11-01; domains/pfam_v47/pfam_domains_rte.tsv): thin bars above the ORF
    pf = W / "domains" / "pfam_v47" / "pfam_domains_rte.tsv"; PF = {}
    if pf.exists():
        for r in pd.read_csv(pf, sep="	").itertuples():
            if float(r.dom_ievalue) < 1e-4 and (r.domain not in PF or float(r.dom_ievalue) < PF[r.domain][2]): PF[r.domain] = (int(r.ali_from), int(r.ali_to), float(r.dom_ievalue), r.pfam.split(".")[0])
    PFY = BT - 0.9
    for dom, (a0, a1, ev, acc) in sorted(PF.items(), key=lambda kv: kv[1][0]):
        if dom.startswith("Exo_endo_phos_2"): continue                                     # the shorter model of the same AP-endonuclease domain
        x0, x1 = X(3 * (a0 - 1) + 1), X(3 * a1); ov.add_patch(Rectangle((x0, PFY - 0.7), x1 - x0, 0.7, fc=ST.NAVY, ec="none"))
        left = dom.startswith("Exo_endo"); ov.text(x0 if left else X(L), PFY - 0.9, f"{dom} {acc}, aa {a0}–{a1}" + (" (Pfam)" if left else ""), fontsize=5.8, ha="left" if left else "right", va="bottom", color=ST.NAVY)
    ov.text(4.2, 38.4, f"reverse transcriptase: best match RTE clade ({top.pid:.0f} % identity), next RTEX", fontsize=7, ha="left", va="top")
    ov.text(4.2, 42.2, f"no target-site duplications ({t2} of {D['n']} copies)", fontsize=7, ha="left", va="top")
    ov.text(4.2, 46.0, f"3′ end …{D['tail']}, an (AAG)n tail", fontsize=7, ha="left", va="top")
    top = H.iloc[0]
    return {"consensus_bp": L, "full_length_copies": D["n"], "ORF_aa": len(D["orf"]), "ORF_LINE_from": 1, "ORF_LINE_to": D["orf_end"], "first_Met_aa": D["first_met"], "first_ATG_LINE_pos": 3 * D["first_met"] - 2,
            "ORF_aa_from_first_Met": len(D["orf"]) - D["first_met"] + 1, "APE_motif": D["ape_motif"], "APE_motif_aa": D["ape_aa"], "APE_motif_LINE_pos": 3 * D["ape_aa"] - 2,
            "RT_aa_from": int(top.q0), "RT_aa_to": int(top.q1), "RT_LINE_from": 3 * int(top.q0) - 2, "RT_LINE_to": 3 * int(top.q1), "RT_best_reference": top.ref, "RT_identity_pct": round(float(top.pid), 1),
            "RT_raw_score": float(top.score), "RT_E_estimate": f"{top.E:.1e}", "RT_second_reference": H.ref[1], "RT_second_E_estimate": f"{H.E[1]:.1e}", "RT_third_reference": H.ref[2], "RT_third_E_estimate": f"{H.E[2]:.1e}",
            "TSD_copies_ge5bp": D["tsd_n"][5], "TSD_copies_ge6bp": D["tsd_n"][6], "TSD_copies_ge8bp": D["tsd_n"][8], "TSD_lengths_bp": ",".join(map(str, D["tsd_len"])),
            "3prime_after_stop_bp": D["utr3"], "3prime_end_LINE_sense": D["tail"]}


if __name__ == "__main__":
    ST.house_style(); Wm, Hm = 88.0, 52.0
    fig = ST.new_figure(Wm, Hm); res = draw(fig, 0, 0, Wm, Hm, "a")
    ST.save(fig, str(HERE / "_standalone" / "fig2_panel_a"))
    pd.DataFrame({"metric": list(res), "value": list(res.values())}).to_csv(HERE / "fig2_panel_a_numbers.tsv", sep="\t", index=False)
    print(res)
