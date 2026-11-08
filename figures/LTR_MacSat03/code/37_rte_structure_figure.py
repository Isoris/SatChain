"""Step 37 (2026-10-24): structural diagram of the RTE LINE interleaved with MacSat03.
Coordinates: LINE sense, 5' -> 3' left to right; position p = 3,118 - q, q = source-consensus coordinate (seq/Gar_full_consensus.fa, 3,117 bp,
written in the satellite-array orientation where q = 1 is the LINE 3' end).
A element map (ORF, RT domain, APE-like motif, UTRs, tail); B correspondence with the 8,164-bp EDTA entry TE_00003409_INT (blastn segments,
step 32); C what each centromere keeps (100-bp bins, >= 50 % of units, step 32n) + Gar dispersed full-length copies; D the element inside
the composite unit (array orientation, median lengths from step 33)."""
import pathlib, re, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/code"); import ltr_colours as LC
from matplotlib.patches import Rectangle, FancyArrow, Patch
W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22")
L = 3117
cons = "".join(l.strip() for l in open(W / "seq" / "Gar_full_consensus.fa") if not l.startswith(">"))
orf = [l.strip() for l in open(W / "seq" / "Gar_full_orfs_minus.faa")][1]
assert len(cons) == L
P = lambda q: L + 1 - q                                    # source coordinate -> LINE-sense position
aa_nt = lambda a: (3 * a - 2, 3 * a)                      # ORF starts at LINE position 1 (frame 0 of the minus strand)
orf_end = 3 * len(orf); first_M = orf.find("M") + 1
m = re.search(r"[CVIL][LVI]QET", orf); ape_aa = m.start() + 1 if m else None
tail = cons[:40].translate(str.maketrans("ACGT", "TGCA"))[::-1]   # LINE-sense 3' end
RT = (548, 954)
C_TYPE = {"I": LC.TYPE1, "II": LC.TYPE2}
fig = plt.figure(figsize=(16, 15)); fig.patch.set_facecolor("white")

# ---------------- A: element map ----------------
ax = fig.add_axes([0.06, 0.80, 0.90, 0.15]); ax.set_xlim(-60, L + 60); ax.set_ylim(-1.6, 3.5); ax.axis("off")
ax.add_patch(Rectangle((1, 0), L, 0.5, color="#cfcfcf"))
ax.add_patch(Rectangle((1, 0.6), orf_end, 0.7, color=LC.ELEM))
ax.add_patch(Rectangle((aa_nt(RT[0])[0], 0.6), aa_nt(RT[1])[1] - aa_nt(RT[0])[0], 0.7, color=LC.ELEM_CORE))
ax.add_patch(Rectangle((1, 0.6), aa_nt(RT[0])[0] - 1, 0.7, fill=False, hatch="///", edgecolor=LC.ELEM_CORE, lw=0))
if ape_aa: ax.plot([aa_nt(ape_aa)[0]] * 2, [0.55, 1.4], color="black", lw=1.5); ax.text(aa_nt(ape_aa)[0], 1.5, f"APE-like motif '{m.group()}'\n(aa {ape_aa}; sequence motif only)", ha="center", va="bottom", fontsize=10)
ax.text(aa_nt(RT[0])[0] + (aa_nt(RT[1])[1] - aa_nt(RT[0])[0]) / 2, 0.95, f"reverse transcriptase (aa {RT[0]}-{RT[1]})\nbest match RTE clade, E = 7e-88", ha="center", va="center", fontsize=10.5, color="white", fontweight="bold")
ax.text(aa_nt(RT[0])[0] / 2, 0.95, "N-terminal region, expected" + chr(10) + "APE endonuclease (not tested)", ha="center", va="center", fontsize=9, bbox=dict(facecolor=LC.ELEM, edgecolor="none", pad=0.6))
ax.text(orf_end / 2 - 300, 2.8, f"one ORF, {len(orf):,} aa (LINE positions 1-{orf_end:,}); open at the 5' end: no start codon inside the consensus (first Met at aa {first_M})", ha="center", fontsize=11)
ax.add_patch(Rectangle((orf_end + 1, 0), L - orf_end, 0.5, color="#7a7a7a"))
ax.annotate(f"3' end ...{tail[-18:]}: (AAG)n tail", xy=(L, 0.25), xytext=(L - 120, 2.0), fontsize=10, arrowprops=dict(arrowstyle="->"), ha="right", va="center")
ax.text(1, -0.55, "5'", fontsize=12, fontweight="bold", ha="center")
ax.plot([91, 91], [0.55, 1.4], color=LC.INK, lw=1.5)
ax.text(120, -1.45, "first ATG at LINE 91 (double ATG, same site in Mac); 5' UTR ~90-115 bp: copies share the 5' end within +-25 bp; "
        "no CCAAT / GC box; ORF from Met = 997 aa", fontsize=9.5, color=LC.INK, ha="left", va="center"); ax.text(L, -0.55, "3'", fontsize=12, fontweight="bold", ha="center")
for x in range(0, L + 1, 500): ax.plot([max(x, 1)] * 2, [-0.1, -0.25], color="black", lw=0.8); ax.text(max(x, 1), -0.8, f"{x:,}", ha="center", fontsize=9)
ax.text(-60, 3.4, "A  The source element: dispersed RTE LINE, 3,117 bp (consensus of 37 full-length C. gariepinus copies); no target-site duplications (2 / 37)",
        fontsize=13, fontweight="bold", ha="left")

# ---------------- B: EDTA correspondence ----------------
segs = [(6, 518, 5578, 6068), (694, 821, 6073, 6200), (887, 1574, 6203, 6844), (1570, 1699, 7265, 7388), (2681, 3115, 7387, 7795)]
ax = fig.add_axes([0.06, 0.60, 0.90, 0.14]); ax.set_xlim(-60, 8164 + 60); ax.set_ylim(-1.2, 3.6); ax.axis("off")
ax.add_patch(Rectangle((0, 0), 8164, 0.5, color="#e6e6e6")); ax.add_patch(Rectangle((2929, 0), 3388 - 2929, 0.5, color=LC.SAT))
ax.text(4100, -0.6, "EDTA TE_00003409_INT (8,164 bp, 'LTR/unknown'): built from the centromeric UNIT, not from the TE", ha="center", fontsize=10.5)
ax.text(3158, 0.75, "MacSat03\nstretch", ha="center", fontsize=8.5, color=LC.SAT)
sx = lambda p: 5550 + (p - 1) * (7820 - 5550) / L            # LINE positions drawn on a source bar above the EDTA part they match
ax.add_patch(Rectangle((sx(1), 2.0), sx(L) - sx(1), 0.5, color=LC.ELEM))
ax.text(sx(1) - 60, 2.25, "source RTE\n(5' -> 3')", ha="right", va="center", fontsize=9.5)
for q0, q1, e0, e1 in segs:
    p0, p1 = P(q1), P(q0)
    ax.fill([sx(p0), sx(p1), e0 if False else max(e0, e1), min(e0, e1)], [2.0, 2.0, 0.5, 0.5], color=LC.P.shade(LC.ELEM, 0.6), lw=0)
    ax.add_patch(Rectangle((min(e0, e1), 0), abs(e1 - e0), 0.5, color=LC.ELEM))
p0, p1 = P(2680), P(1700); ax.add_patch(Rectangle((sx(p0), 2.0), sx(p1) - sx(p0), 0.5, fill=False, hatch="xxx", edgecolor=LC.INK, lw=0))
ax.text((sx(p0) + sx(p1)) / 2, 2.65, "source 1,700-2,680 absent from EDTA", ha="center", fontsize=9, color=LC.INK)
ax.text(-60, 3.45, "B  Where the source element sits in the EDTA entry (blastn segments; reversed: the EDTA entry is written in array orientation)", fontsize=13, fontweight="bold", ha="left")

# ---------------- C: per-centromere state ----------------
B = pd.read_csv(W / "homologous" / "source_bins_per_chrom.tsv", sep="\t", index_col=0)
S = pd.read_csv(W / "split" / "v46_fragment_summary.tsv", sep="\t").set_index("chrom")
T2 = {"LG08", "LG16", "LG25"}
order = ["LG04", "LG02", "LG21", "LG23", "LG01", "LG07", "LG19", "LG15", "LG26", "LG08", "LG16", "LG25"]
order = [c for c in order if c in B.index]
ax = fig.add_axes([0.06, 0.23, 0.90, 0.33]); ax.set_xlim(-430, L + 60); ax.set_ylim(len(order) + 1.2, -2.5); ax.axis("off")
ax.add_patch(Rectangle((1, -0.95), L, 0.6, color=LC.ELEM)); ax.text(-20, -0.65, "Gar dispersed\nfull-length copies", ha="right", va="center", fontsize=9.5)
for i, c in enumerate(order):
    y = i + 0.1; typ = "II" if c in T2 else "I"
    for b, v in enumerate(B.loc[c].values[:31]):
        if v >= 0.5:
            q0, q1 = b * 100 + 1, (b + 1) * 100
            ax.add_patch(Rectangle((P(q1), y), 100, 0.6, color=C_TYPE[typ]))
    ax.plot([1, L], [y + 0.3] * 2, color="#cccccc", lw=0.6, zorder=0)
    lab = f"{c}  type {typ}"
    if c in S.index: lab += f"  ({int(S.loc[c, 'spacers'])} spacers, {S.loc[c, 'pid']:.0f} %)"
    ax.text(-20, y + 0.3, lab, ha="right", va="center", fontsize=9.5)
for q, name in ((1698, "shared deletion end: LINE 1,420\n(source 1,698)"), (2686, "deletion start: LINE 432\n(source 2,686)")):
    ax.plot([P(q)] * 2, [-0.2, len(order)], color="black", lw=0.8, ls="--"); ax.text(P(q), len(order) + 0.15, name, ha="center", va="top", fontsize=9)
for x in range(0, L + 1, 500): ax.text(max(x, 1), -1.25, f"{x:,}", ha="center", fontsize=8.5, color="#555555")
ax.text(-430, -2.4, "C  What each MacSat03 centromere keeps (LINE positions; 100-bp bins present in >= 50 % of its units; % = identity to the source)", fontsize=13, fontweight="bold", ha="left")
ax.legend(handles=[Patch(color=C_TYPE["I"], label="type I centromere"), Patch(color=C_TYPE["II"], label="type II centromere (carries the ~6-kb insert)")],
          loc="lower right", fontsize=10, frameon=False, bbox_to_anchor=(1.0, -0.12))

# ---------------- D: inside the composite unit ----------------
sp = pd.read_csv(W / "split" / "spacer_summary.tsv", sep="\t").set_index("chrom")
num = lambda s: float(str(s).split()[0])
def med(rows, col): return float(np.median([num(sp.loc[c, col]) for c in rows if c in sp.index]))
t1 = [c for c in order if c not in T2 and c in sp.index]; t2 = [c for c in order if c in T2 and c in sp.index]
ax = fig.add_axes([0.06, 0.03, 0.90, 0.16]); ax.axis("off"); ax.set_ylim(-1.6, 2.4)
def unit(y, parts, title):
    x = 0
    for name, ln, col in parts:
        if name == "RTE":
            ax.add_patch(FancyArrow(x + ln, y + 0.3, -ln, 0, width=0.6, head_width=0.6, head_length=min(500, ln * 0.25), length_includes_head=True, color=col, lw=0))
            ax.text(x + ln / 2, y + 0.3, "RTE (3' tail on the left)" if ln > 2000 else "RTE", ha="center", va="center", fontsize=9, color="black")
        else:
            ax.add_patch(Rectangle((x, y), ln, 0.6, color=col)); ax.text(x + ln / 2, y + 0.3, name, ha="center", va="center", fontsize=9, color="white" if col != "#d9d9d9" else "black")
        x += ln
    ax.text(-200, y + 0.3, title, ha="right", va="center", fontsize=10.5)
    return x
sat = 4400
xa = unit(1.2, [("MacSat03 block", sat, LC.SAT), ("flank", med(t1, "gap_sat_to_te_left"), "#d9d9d9"), ("RTE", med(t1, "src_bp_median"), LC.ELEM),
                ("flank", med(t1, "gap_te_to_sat_right"), "#d9d9d9"), ("MacSat03", 1500, LC.SAT)], "type I unit")
xb = unit(0.0, [("MacSat03 block", sat, LC.SAT), ("type II insert", med(t2, "ins_bp_median"), LC.INS), ("RTE", med(t2, "src_bp_median"), LC.ELEM),
                ("", med(t2, "gap_te_to_sat_right"), "#d9d9d9"), ("MacSat03", 1500, LC.SAT)], "type II unit")
ax.set_xlim(-3200, max(xa, xb) + 200)
ax.text(-3200, 2.2, "D  The element inside the composite centromeric unit (array orientation, median lengths; the RTE reads 3' -> 5' along the array)",
        fontsize=13, fontweight="bold", ha="left")
ax.plot([0, 5000], [-1.0, -1.0], color="black"); ax.text(2500, -1.45, "5 kb", ha="center", fontsize=9.5)
out = W / "figure" / "rte_structure"
for ext in ("png", "pdf"): fig.savefig(f"{out}.{ext}", dpi=450)
print("ok", "ORF", len(orf), "aa; first Met", first_M, "; APE-like motif", m.group() if m else None, "at aa", ape_aa)
