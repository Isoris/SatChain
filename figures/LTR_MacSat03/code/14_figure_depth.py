"""Figure 3 (2026-10-23; user: one level deeper, inside monomers and elements).
A element map: ORFs >= 100 aa (6 frames), RT domain (tblastn), built-in satellite stretch, internal repeats (self-blast pairs as arcs);
  below, share of units carrying each 100-bp node, type I vs type II.
B two consecutive units at monomer resolution, LG01 (type I) and LG25 (type II): each MacSat03 monomer (36_monomers_final, promoted) as a box coloured by
  subfamily, element hits as arrows labelled with their element coordinates, insert in red.
C MacSat03 monomers inside units: subfamily composition per centromere and identity to own consensus.
Output: figure/composite_unit_depth.png/.pdf"""
import pathlib, re, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/code"); import ltr_colours as LC
from matplotlib.patches import Rectangle, FancyArrow, Arc

W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22"); O = pathlib.Path("F:/data/satellite_FISH/pipeline_v4_5_fixed/out/fClaHyb_Mac_satellites")
T1 = ["C_mac_LG07", "C_mac_LG04", "C_mac_LG01", "C_mac_LG23", "C_mac_LG02"]; T2 = ["C_mac_LG08", "C_mac_LG25", "C_mac_LG16"]
C1, C2, CINS, CEL, CCORE = LC.TYPE1, LC.TYPE2, LC.INS, LC.ELEM, LC.ELEM_CORE
lab = lambda c: "LG" + c.split("LG")[-1]
seq = "".join(l.strip() for l in open(W / "seq" / "TE_00003409_library.fa") if not l.startswith(">")).upper(); L = len(seq)

# ---- ORFs (6 frames, >= 100 aa, stop to stop)
code = {a + b + c: aa for (a, b, c), aa in zip(((x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
rc = seq[::-1].translate(str.maketrans("ACGT", "TGCA"))
orfs = []
for strand, s in (("+", seq), ("-", rc)):
    for f in range(3):
        aa = "".join(code.get(s[i:i + 3], "X") for i in range(f, len(s) - 2, 3))
        pos = 0
        for m in re.finditer(r"[^*]+", aa):
            if m.end() - m.start() >= 100:
                a, b = f + 3 * m.start(), f + 3 * m.end()
                orfs.append((strand, a, b) if strand == "+" else (strand, L - b, L - a))
self_ = pd.read_csv(W / "domains" / "element_self.tsv", sep="\t", header=None, names=["q0", "q1", "s0", "s1", "pid", "len"])
self_ = self_[self_.q0 < self_.s0]

U = pd.read_csv(W / "unit" / "units.tsv", sep="\t")
mono = pd.read_csv(O / "36_monomers_final.tsv", sep="\t", usecols=["chrom", "genomic_start_0based", "genomic_end", "strand", "family_id", "subfamily_id",
                                                                    "best_alignment_identity", "complete_or_truncated", "promoted"])
mono = mono[(mono.family_id == "MacSat03") & (mono.chrom.isin(T1 + T2))]
eh = pd.read_csv(W / "copies" / "Mac.hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
eh["lo"] = eh[["s0", "s1"]].min(1); eh["hi"] = eh[["s0", "s1"]].max(1)
ih = pd.read_csv(W / "insert" / "Mac.insert_hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
ih["lo"] = ih[["s0", "s1"]].min(1); ih["hi"] = ih[["s0", "s1"]].max(1); ih = ih[ih.len >= 300]

# monomers inside units: subfamily share per chromosome
rows = []
for r in U[U.chrom.isin(T1 + T2)].itertuples():
    m = mono[(mono.chrom == r.chrom) & (mono.genomic_start_0based >= r.start) & (mono.genomic_end <= r.start + r.period)]
    for x in m.itertuples(): rows.append((r.chrom, x.subfamily_id, x.best_alignment_identity))
M = pd.DataFrame(rows, columns=["chrom", "subfamily", "identity"])
top_sf = M.subfamily.value_counts().index[:9].tolist()
pal = ["#1f4e9c", "#1b9e8a", "#7b4fa3", "#c0392b", "#7f9f2f", "#4a90d9", "#b5468a", "#2e8b3d", "#8c6a4f"]   # muted opaque sequence (= class_colours.CLASS_SEQ[:9]): subfamilies are categories inside one element
SFC = {s: pal[i] for i, s in enumerate(top_sf)}; sfc = lambda s: SFC.get(s, "#c8c8c8")

fig = plt.figure(figsize=(17, 17), facecolor="white")
gs = fig.add_gridspec(3, 2, height_ratios=[1.05, 1.25, 1.0], width_ratios=[1.6, 1.0], hspace=0.38, wspace=0.18)

# ---- A element map
ax = fig.add_subplot(gs[0, :]); kb = lambda x: x / 1000
ax.add_patch(Rectangle((0, 2.6), kb(L), 0.5, facecolor=CEL, edgecolor="none"))
ax.add_patch(Rectangle((kb(5560), 2.6), kb(1840), 0.5, facecolor=CCORE, edgecolor="none"))
ax.add_patch(Rectangle((kb(2929), 2.6), kb(460), 0.5, facecolor=LC.SAT, edgecolor="none"))
ax.text(kb(3158), 3.25, "built-in MacSat03\n(~2 monomers)", ha="center", fontsize=9, color=LC.SAT)
ax.text(kb(6480), 3.25, "RT core (used to find copies)", ha="center", fontsize=9, color=CCORE)
ax.add_patch(Rectangle((kb(6219), 2.05), kb(527), 0.35, facecolor=LC.INK, edgecolor="none")); ax.text(kb(6800), 2.22, "RT domain (RTE-clade, E = 3e-34)", fontsize=9, va="center")
for st, a, b in orfs:
    y = 1.5 if st == "+" else 1.0
    ax.add_patch(FancyArrow(kb(a if st == "+" else b), y, kb(b - a) * (1 if st == "+" else -1), 0, width=0.22, head_width=0.32, head_length=0.08, length_includes_head=True, color="#555"))
ax.text(-0.15, 1.5, "ORFs + strand", ha="right", va="center", fontsize=9); ax.text(-0.15, 1.0, "ORFs − strand", ha="right", va="center", fontsize=9)
for r in self_.itertuples():
    c1, c2 = kb((r.q0 + r.q1) / 2), kb((r.s0 + r.s1) / 2)
    ax.add_patch(Arc(((c1 + c2) / 2, 3.1), c2 - c1, 1.4 * min(1, (c2 - c1) / 4), theta1=0, theta2=180, color="#888", lw=1.2))
ax.text(kb(L) + 0.1, 3.55, "arcs = internal repeats\n(self-similarity ≥ 68%)", fontsize=8, color="#666", va="center")
NP = pd.read_csv(W / "unit" / "node_presence.tsv", sep="\t", index_col=0)
x = np.arange(NP.shape[1]) * 0.1 + 0.05
for typ, chs, col, y0 in (("type I", T1, C1, -0.3), ("type II", T2, C2, -1.6)):
    v = NP.loc[[c for c in chs if c in NP.index]].mean().values
    ax.fill_between(x, y0, y0 + v * 1.0, color=col, step="mid"); ax.text(-0.15, y0 + 0.5, f"kept in\n{typ} units", ha="right", va="center", fontsize=9, color=col)
    ax.add_line(plt.Line2D([0, kb(L)], [y0, y0], color="#bbb", lw=0.8))
ax.set_xlim(-1.3, kb(L) + 1.6); ax.set_ylim(-1.8, 4.6); ax.set_yticks([]); ax.set_xlabel("position along the element consensus (kb)", fontsize=11)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.set_title("A  Inside the element: reading frames, RT domain, built-in satellite, internal repeats; which parts each unit type keeps", loc="left", fontsize=13, fontweight="bold")

# ---- B monomer-resolution zoom
ax = fig.add_subplot(gs[1, :])
for i, ch in enumerate(("C_mac_LG01", "C_mac_LG25")):
    u = U[U.chrom == ch].reset_index(drop=True); j = int((u.period - u.period.median()).abs().idxmin()); j = min(j, len(u) - 2)
    a = int(u.loc[j, "start"]); b = int(u.loc[j + 1, "start"] + u.loc[j + 1, "period"]); y = 1 - i
    ax.add_patch(Rectangle((0, y - 0.32), kb(b - a), 0.64, facecolor="#eceae4", edgecolor="none"))
    m = mono[(mono.chrom == ch) & (mono.genomic_end > a) & (mono.genomic_start_0based < b)]
    for r in m.itertuples():
        ax.add_patch(Rectangle((kb(r.genomic_start_0based - a), y - 0.3), kb(r.genomic_end - r.genomic_start_0based), 0.25, facecolor=sfc(r.subfamily_id), edgecolor="white", lw=0.3))
    for r in ih[(ih.chrom == ch) & (ih.hi > a) & (ih.lo < b)].itertuples():
        ax.add_patch(Rectangle((kb(max(r.lo, a) - a), y + 0.0), kb(min(r.hi, b) - max(r.lo, a)), 0.28, facecolor=CINS, edgecolor="none"))
    e = eh[(eh.chrom == ch) & (eh.hi > a) & (eh.lo < b) & (eh.len >= 150)].sort_values("lo")
    for r in e.itertuples():
        core = r.q0 < 7400 and r.q1 > 5560; sat = r.q0 >= 2900 and r.q1 <= 3420
        if sat: continue
        fwd = r.s1 >= r.s0; x0, x1 = kb(max(r.lo, a) - a), kb(min(r.hi, b) - a)
        ax.add_patch(FancyArrow(x0 if fwd else x1, y + 0.14, (x1 - x0) * (1 if fwd else -1), 0, width=0.24, head_width=0.3, head_length=min(0.12, (x1 - x0) * 0.4),
                                length_includes_head=True, color=CCORE if core else CEL, ec="none"))
        if r.len >= 600: ax.text((x0 + x1) / 2, y + 0.36, f"{r.q0 / 1000:.1f}–{r.q1 / 1000:.1f}", ha="center", fontsize=7, color="#333")
    for xb in (a, int(u.loc[j + 1, "start"]), b):
        ax.add_line(plt.Line2D([kb(xb - a)] * 2, [y - 0.45, y + 0.5], color="black", lw=0.9, zorder=6))
    ax.text(-0.3, y, f"{lab(ch)}\n({'type I' if ch in T1 else 'type II'})", ha="right", va="center", fontsize=11, fontweight="bold", color=C1 if ch in T1 else C2)
    ax.text(kb(b - a) + 0.3, y, f"2 units, {kb(b - a):.1f} kb\n{len(m)} MacSat03 monomers", va="center", fontsize=9, color="#444")
ax.set_xlim(-3.5, 34); ax.set_ylim(-0.7, 1.75); ax.set_yticks([]); ax.set_xlabel("kb (two consecutive units)", fontsize=11)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.set_title("B  Two units at monomer resolution: each box = one MacSat03 monomer (colour = subfamily); arrows = element pieces (label = element kb)", loc="left", fontsize=13, fontweight="bold")
for k, s in enumerate(top_sf):
    ax.add_patch(Rectangle((k * 3.4, 1.55), 0.5, 0.12, facecolor=SFC[s], clip_on=False)); ax.text(k * 3.4 + 0.6, 1.61, s, fontsize=8, va="center")
ax.add_patch(Rectangle((31, 1.55), 0.5, 0.12, facecolor="#c8c8c8", clip_on=False)); ax.text(31.6, 1.61, "other", fontsize=8, va="center")

# ---- C subfamily composition per centromere
ax = fig.add_subplot(gs[2, 0])
order = T1 + T2; bottom = np.zeros(len(order))
tab = M.assign(sf=M.subfamily.where(M.subfamily.isin(top_sf), "other")).groupby(["chrom", "sf"]).size().unstack(fill_value=0).reindex(order).fillna(0)
tab = tab.div(tab.sum(1), axis=0)
for s in top_sf + ["other"]:
    if s in tab: v = tab[s].values; ax.bar(range(len(order)), v, bottom=bottom, color=sfc(s) if s != "other" else "#c8c8c8", width=0.75); bottom += v
ax.set_xticks(range(len(order))); ax.set_xticklabels([f"{lab(c)}\n{'I' if c in T1 else 'II'}" for c in order], fontsize=10)
for t, c in zip(ax.get_xticklabels(), order): t.set_color(C1 if c in T1 else C2)
ax.set_ylabel("share of MacSat03 monomers in the units", fontsize=10); ax.set_ylim(0, 1)
ax.set_title("C  Monomer subfamilies inside the units, per centromere", loc="left", fontsize=13, fontweight="bold")

# ---- D monomer identity
ax = fig.add_subplot(gs[2, 1])
data = [M[M.chrom == c].identity.values for c in order]
bp = ax.boxplot(data, widths=0.6, showfliers=False, patch_artist=True)
for p, c in zip(bp["boxes"], order): p.set_facecolor(C1 if c in T1 else C2)
ax.set_xticks(range(1, len(order) + 1)); ax.set_xticklabels([lab(c) for c in order], fontsize=9, rotation=45)
ax.set_ylabel("monomer identity to its class consensus (%)", fontsize=10)
ax.set_title("D  Monomer identity", loc="left", fontsize=13, fontweight="bold")
fig.suptitle("Inside the composite unit: element structure and satellite monomers", fontsize=16, fontweight="bold", x=0.02, ha="left", y=0.995)
for ext in ("png", "pdf"): fig.savefig(W / "figure" / f"composite_unit_depth.{ext}", dpi=450, bbox_inches="tight")
print("ORFs >= 100 aa:", orfs); print("monomers in units:", len(M), "| subfamilies:", M.subfamily.nunique())
print(tab.round(2).to_string())
