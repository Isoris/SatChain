"""Figure (2026-10-23): the composite centromeric unit of C. macrocephalus lineage A centromeres.
A representative unit per chromosome to scale (median-period unit; oriented so the element core reads forward): MacSat03 satellite, element segments
coloured by element coordinate, the type II insert. B element node presence (100-bp bins) per chromosome. C unit period vs satellite content.
D MacSat03 chromosome-consensus identity by unit type. E within-species copy distance vs genomic separation (homogenization).
Output: figure/composite_unit.png/.pdf"""
import pathlib, itertools, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/code"); import ltr_colours as LC
from matplotlib.patches import Rectangle
from matplotlib import cm, colors

W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22"); R = pathlib.Path("F:/data/satellite_FISH/runs_of_record_v4.5.3_2026-10-18")
(W / "figure").mkdir(exist_ok=True)
T1 = ["C_mac_LG07", "C_mac_LG04", "C_mac_LG01", "C_mac_LG23", "C_mac_LG02"]; T2 = ["C_mac_LG08", "C_mac_LG25", "C_mac_LG16"]; ORDER = T1 + T2
C1, C2, CSAT, CINS, COTH = LC.TYPE1, LC.TYPE2, LC.SAT, LC.INS, LC.OTHER
lab = lambda c: "LG" + c.split("LG")[-1]
U = pd.read_csv(W / "unit" / "units.tsv", sep="\t")
cp = pd.read_csv(W / "core" / "Mac.core_copies.tsv", sep="\t").set_index(["chrom", "lo"])
runs = pd.read_csv(R / "Mac" / "runs_of_record.tsv", sep="\t"); runs = runs[(runs.unit_class == "satellite") & runs.tier.isin(["S1", "S2"])]
eh = pd.read_csv(W / "copies" / "Mac.hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
eh = eh[~((eh.q0 >= 2900) & (eh.q1 <= 3420))]; eh["lo"] = eh[["s0", "s1"]].min(1); eh["hi"] = eh[["s0", "s1"]].max(1)
ih = pd.read_csv(W / "insert" / "Mac.insert_hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
ih["lo"] = ih[["s0", "s1"]].min(1); ih["hi"] = ih[["s0", "s1"]].max(1)
cmap = colors.LinearSegmentedColormap.from_list("elem_amber", ["#f6e2a0", LC.ELEM, LC.ELEM_CORE]); norm = colors.Normalize(0, 8164)   # element coordinate: light amber -> amber -> dark amber (the RT core), one hue = the element

fig = plt.figure(figsize=(16, 13), facecolor="white")
gs = fig.add_gridspec(3, 3, height_ratios=[1.35, 1.0, 1.0], hspace=0.42, wspace=0.3)

# ---- A: representative units
ax = fig.add_subplot(gs[0, :])
for i, ch in enumerate(ORDER):
    u = U[U.chrom == ch]; r = u.iloc[(u.period - u.period.median()).abs().argsort().iloc[0]]
    a, b = int(r.start), int(r.start + r.period)
    st = cp.loc[(ch, a), "strand"] if (ch, a) in cp.index else "+"
    X = (lambda p: p - a) if st == "+" else (lambda p: b - p)
    y = len(ORDER) - 1 - i
    ax.add_patch(Rectangle((0, y - 0.3), b - a, 0.6, facecolor=COTH, edgecolor="none"))
    for s, e in zip(runs[(runs.chrom == ch) & (runs.end > a) & (runs.start0 < b)].start0, runs[(runs.chrom == ch) & (runs.end > a) & (runs.start0 < b)].end):
        x0, x1 = sorted((X(max(s, a)), X(min(e, b)))); ax.add_patch(Rectangle((x0, y - 0.3), x1 - x0, 0.6, facecolor=CSAT, edgecolor="none"))
    g = ih[(ih.chrom == ch) & (ih.hi > a) & (ih.lo < b) & (ih.len >= 300)]
    for s, e in zip(g.lo, g.hi):
        x0, x1 = sorted((X(max(s, a)), X(min(e, b)))); ax.add_patch(Rectangle((x0, y - 0.22), x1 - x0, 0.44, facecolor=CINS, edgecolor="none"))
    g = eh[(eh.chrom == ch) & (eh.hi > a) & (eh.lo < b)]
    for s, e, q0, q1 in zip(g.lo, g.hi, g.q0, g.q1):
        x0, x1 = sorted((X(max(s, a)), X(min(e, b)))); ax.add_patch(Rectangle((x0, y - 0.3), x1 - x0, 0.6, facecolor=cmap(norm((q0 + q1) / 2)), edgecolor="none"))
    ax.text(-250, y, f"{lab(ch)}  ({'I' if ch in T1 else 'II'})", ha="right", va="center", fontsize=11, color=C1 if ch in T1 else C2, fontweight="bold")
    ax.text(b - a + 200, y, f"{(b - a) / 1000:.1f} kb", va="center", fontsize=10, color="#444")
ax.set_xlim(-2200, 19500); ax.set_ylim(-0.8, len(ORDER) - 0.2); ax.set_yticks([]); ax.set_xlabel("position in the unit (bp), from the element core onward", fontsize=11)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.set_title("A  One representative unit per centromere, to scale", loc="left", fontsize=13, fontweight="bold")
hx = 15700
for k, (c, t) in enumerate(((CSAT, "MacSat03 satellite"), (CINS, "type II insert (~6 kb)"), (COTH, "other"))):
    ax.add_patch(Rectangle((hx, len(ORDER) - 1.3 - k * 0.7), 300, 0.4, facecolor=c, clip_on=False)); ax.text(hx + 380, len(ORDER) - 1.1 - k * 0.7, t, fontsize=9, va="center")
cax = ax.inset_axes([0.87, 0.5, 0.11, 0.035]); fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
cax.set_title("position in the element (bp)", fontsize=9); cax.tick_params(labelsize=8)

# ---- B: node presence
ax = fig.add_subplot(gs[1, :2])
NP = pd.read_csv(W / "unit" / "node_presence.tsv", sep="\t", index_col=0).loc[ORDER]
im = ax.imshow(NP.values, aspect="auto", cmap=colors.LinearSegmentedColormap.from_list("node_amber", ["#ffffff", LC.ELEM, LC.ELEM_CORE]), vmin=0, vmax=1, interpolation="nearest")
ax.set_yticks(range(len(ORDER))); ax.set_yticklabels([f"{lab(c)} ({'I' if c in T1 else 'II'})" for c in ORDER], fontsize=10)
for t, c in zip(ax.get_yticklabels(), ORDER): t.set_color(C1 if c in T1 else C2)
ax.set_xticks(range(0, NP.shape[1], 10)); ax.set_xticklabels([f"{i * 100}" for i in range(0, NP.shape[1], 10)], fontsize=9)
[ax.plot([x_, x_], [-0.5, len(ORDER) - 0.5], color=c_, lw=1.6, ls="--", zorder=5) for x_, c_ in ((29, LC.SAT), (33.9, LC.SAT), (62, LC.INK), (67.5, LC.INK))]
ax.text(31.5, -0.6, "satellite stretch", ha="center", va="bottom", clip_on=False, fontsize=8, color=LC.SAT, bbox=dict(facecolor="white", edgecolor="none", pad=1)); ax.text(64.7, -0.6, "RT domain", ha="center", va="bottom", clip_on=False, fontsize=8, color=LC.INK, bbox=dict(facecolor="white", edgecolor="none", pad=1))
ax.set_xlabel("element node (100-bp bin; start position in the element consensus, bp)", fontsize=10)
ax.set_title("B  Element nodes carried by each centromere's units (share of units)", loc="left", fontsize=13, fontweight="bold", pad=18)
fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)

# ---- C: period vs satellite
ax = fig.add_subplot(gs[1, 2])
for ch in ORDER:
    u = U[U.chrom == ch]; ax.scatter(u.MacSat03_monomers, u.period / 1000, s=18, color=C1 if ch in T1 else C2, edgecolor="white", linewidths=0.3)
    ax.text(u.MacSat03_monomers.median(), u.period.median() / 1000 + 0.35, lab(ch), fontsize=8, ha="center", color=C1 if ch in T1 else C2)
ax.set_xlabel("MacSat03 monomers in the unit", fontsize=10); ax.set_ylabel("unit length (kb)", fontsize=10); ax.set_ylim(10, 18)
ax.set_title("C  Unit length vs satellite", loc="left", fontsize=13, fontweight="bold")
ax.scatter([], [], color=C1, label="type I: complete element"); ax.scatter([], [], color=C2, label="type II: truncated + insert"); ax.legend(fontsize=9, frameon=False, loc="upper left")

# ---- D: satellite variant similarity
ax = fig.add_subplot(gs[2, 0])
ci = pd.read_csv(R / "chromosome_specific_variants" / "chrom_consensus_identity.tsv", sep="\t"); ci = ci[ci.family == "MacSat03"]
grp = {"II–II": [], "I–I": [], "I–II": [], "with other\ncentromeres": []}
for r in ci.itertuples():
    ta = "II" if r.chrom_a in T2 else "I" if r.chrom_a in T1 else "o"; tb = "II" if r.chrom_b in T2 else "I" if r.chrom_b in T1 else "o"
    k = "-".join(sorted([ta, tb])); k = {"II-II": "II–II", "I-I": "I–I", "I-II": "I–II"}.get(k, "with other\ncentromeres")
    grp[k].append(r.consensus_identity * 100)
ks = list(grp); ax.boxplot([grp[k] for k in ks], widths=0.55, showfliers=False)
for i, k in enumerate(ks): ax.scatter(np.random.default_rng(1).normal(i + 1, 0.06, len(grp[k])), grp[k], s=12, color=[C2, C1, LC.REL[2], LC.REL[3]][i], zorder=3)
ax.set_xticks(range(1, len(ks) + 1)); ax.set_xticklabels(ks, fontsize=9); ax.set_ylabel("MacSat03 consensus identity (%)", fontsize=10)
ax.set_title("D  Satellite variant shared by type", loc="left", fontsize=13, fontweight="bold")

# ---- E: homogenization
ax = fig.add_subplot(gs[2, 1])
from Bio import Phylo
t = Phylo.read(str(W / "tree" / "copies.treefile"), "newick"); t.root_at_midpoint()
tips = [x for x in t.get_terminals() if x.name.startswith("Mac|") and x.name.split("|")[3] == "cen"]
info = [(x.name.split("|")[1], int(x.name.split("|")[2])) for x in tips]
B = {"≤25 kb": [], "25–100 kb": [], "0.1–1 Mb": [], ">1 Mb": [], "other\nchromosome": []}
for i, j in itertools.combinations(range(len(tips)), 2):
    (ci_, pi), (cj, pj) = info[i], info[j]; d = t.distance(tips[i], tips[j])
    if ci_ != cj: B["other\nchromosome"].append(d); continue
    gap = abs(pi - pj); B["≤25 kb" if gap <= 25e3 else "25–100 kb" if gap <= 1e5 else "0.1–1 Mb" if gap <= 1e6 else ">1 Mb"].append(d)
ks = list(B); ax.bar(range(len(ks)), [np.median(B[k]) for k in ks], color=[LC.ELEM] * 3 + [LC.REL[2], LC.REL[2]])
for i, k in enumerate(ks): ax.text(i, np.median(B[k]) + 0.004, f"{np.median(B[k]):.3f}\nn={len(B[k])}", ha="center", fontsize=8)
ax.set_ylim(0, max(np.median(B[k]) for k in ks) * 1.25); ax.set_xticks(range(len(ks))); ax.set_xticklabels(ks, fontsize=9); ax.set_ylabel("median tree distance\nbetween element copies", fontsize=10)
ax.set_xlabel("separation of two centromeric copies", fontsize=10)
ax.set_title("E  Homogenized within each centromere", loc="left", fontsize=13, fontweight="bold")

# ---- F: summary text
ax = fig.add_subplot(gs[2, 2]); ax.axis("off")
txt = ("F  Summary\n\n"
       "• Element: RTE-type LINE (RT E = 3e-34), present in\n  both species; copies amplified separately per species.\n"
       "• Gar: 1,244 dispersed copies (9% centromeric).\n"
       "• Mac: 316 copies, 67% in lineage A centromeres,\n  as a 12–15 kb satellite + element unit.\n"
       "• Type I units (LG07, 04, 01, 23, 02): complete element.\n"
       "• Type II units (LG08, 25, 16): element cut to its RT\n  core + a ~6 kb insert, found ONLY in these 3\n  centromeres (69 copies), not in Gar.\n"
       "• Type II centromeres also share their MacSat03\n  variant (95% vs 83% to type I): whole units spread\n  between non-homologous centromeres, then were\n  homogenized within each.")
ax.text(0, 1, txt, va="top", fontsize=10.5, family="DejaVu Sans")
fig.suptitle("The composite centromeric unit of C. macrocephalus lineage A centromeres", fontsize=16, fontweight="bold", x=0.02, ha="left", y=0.995)
for ext in ("png", "pdf"): fig.savefig(W / "figure" / f"composite_unit.{ext}", dpi=450, bbox_inches="tight")
print("saved")
