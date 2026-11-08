"""Figure 2 (2026-10-23; user: show what is around, prove it repeats, show whether type II exists on other chromosomes).
A: 120-kb window of each of the 8 Mac centromeres (start at a centromeric element core copy): satellite, element, element RT core, type II insert -> the unit repeats.
B: genome-wide insert copies on every chromosome of both haplotypes (>= 3 kb of the 6-kb insert = slate; 1-3 kb = grey; centromere = black bar).
Output: figure/composite_unit_context.png/.pdf"""
import pathlib, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/code"); import ltr_colours as LC
from matplotlib.patches import Rectangle

W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22"); R = pathlib.Path("F:/data/satellite_FISH/runs_of_record_v4.5.3_2026-10-18")
A = pathlib.Path("F:/data/fClaHyb_func_annotation_v1.0")
T1 = ["C_mac_LG07", "C_mac_LG04", "C_mac_LG01", "C_mac_LG23", "C_mac_LG02"]; T2 = ["C_mac_LG08", "C_mac_LG25", "C_mac_LG16"]
CSAT, CEL, CCORE, CINS = LC.SAT, LC.ELEM, LC.ELEM_CORE, LC.INS
lab = lambda c: "LG" + c.split("LG")[-1]
WIN = 120_000
cp = pd.read_csv(W / "core" / "Mac.core_copies.tsv", sep="\t"); cp = cp[cp.in_cen_1Mb].sort_values(["chrom", "lo"])
runs = pd.read_csv(R / "Mac" / "runs_of_record.tsv", sep="\t"); runs = runs[(runs.unit_class == "satellite") & runs.tier.isin(["S1", "S2"])]
eh = pd.read_csv(W / "copies" / "Mac.hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
eh = eh[~((eh.q0 >= 2900) & (eh.q1 <= 3420))]; eh["lo"] = eh[["s0", "s1"]].min(1); eh["hi"] = eh[["s0", "s1"]].max(1)
ch_ = pd.read_csv(W / "core" / "Mac.core_hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
ch_["lo"] = ch_[["s0", "s1"]].min(1); ch_["hi"] = ch_[["s0", "s1"]].max(1)
ih = pd.read_csv(W / "insert" / "Mac.insert_hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
ih["lo"] = ih[["s0", "s1"]].min(1); ih["hi"] = ih[["s0", "s1"]].max(1); ih = ih[ih.len >= 300]

fig = plt.figure(figsize=(16, 19), facecolor="white")
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.55], hspace=0.16, wspace=0.12)

# ---- A: repeating units in context
ax = fig.add_subplot(gs[0, :])
ORDER = T1 + T2
for i, c in enumerate(ORDER):
    g = cp[cp.chrom == c]; lo = g.lo.values
    # start at the copy that begins the densest 120-kb stretch
    k = int(np.argmax([((lo >= x) & (lo < x + WIN)).sum() for x in lo])); a = lo[k] - 2000; b = a + WIN
    y = len(ORDER) - 1 - i; h = 0.62
    ax.add_patch(Rectangle((0, y - h / 2), WIN / 1000, h, facecolor="#eceae4", edgecolor="none"))
    def draw(df, s, e, col, hh, z):
        d = df[(df.chrom == c) & (df[e] > a) & (df[s] < b)]
        for x0, x1 in zip(d[s], d[e]): ax.add_patch(Rectangle(((max(x0, a) - a) / 1000, y - hh / 2), max((min(x1, b) - max(x0, a)) / 1000, 0.08), hh, facecolor=col, edgecolor="none", zorder=z))
    draw(runs, "start0", "end", CSAT, h, 1); draw(eh, "lo", "hi", CEL, h, 2); draw(ch_, "lo", "hi", CCORE, h, 3); draw(ih, "lo", "hi", CINS, h * 0.7, 4)
    n = ((lo >= a) & (lo < b)).sum()
    for x in lo[(lo >= a) & (lo < b)]:                      # unit boundaries = start of each centromeric element core copy
        ax.add_line(plt.Line2D([(x - a) / 1000] * 2, [y - 0.42, y + 0.42], color="black", lw=0.8, zorder=6))
    ax.text(-1.5, y, f"{lab(c)} ({'I' if c in T1 else 'II'})", ha="right", va="center", fontsize=12, fontweight="bold", color=LC.TYPE1 if c in T1 else LC.TYPE2)
    ax.text(WIN / 1000 + 1.5, y, f"{n} units\n{c.split('_')[-1]} {a / 1e6:.2f} Mb", va="center", fontsize=9, color="#444")
ax.set_xlim(-12, WIN / 1000 + 12); ax.set_ylim(-1.7, len(ORDER) - 0.2); ax.set_yticks([])
ax.set_xticks(range(0, WIN // 1000 + 1, 10)); ax.set_xlabel("kb along the centromere (120-kb window)", fontsize=12)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.set_title("A  The unit repeats along each centromere (C. macrocephalus)", loc="left", fontsize=14, fontweight="bold")
for j, (col, t) in enumerate(((CSAT, "MacSat03 satellite"), (CEL, "RTE-like element"), (CCORE, "element RT core"), (CINS, "type II insert"), ("black", "| unit boundary"))):
    ax.add_patch(Rectangle((j * 26, -1.45), 2.5 if col != "black" else 0.3, 0.3, facecolor=col, clip_on=False)); ax.text(j * 26 + 3.3, -1.3, t.lstrip("| "), fontsize=11, va="center")

# ---- B: insert copies genome-wide, both haplotypes
for k, h in enumerate(("Mac", "Gar")):
    ax = fig.add_subplot(gs[1, k])
    sz = pd.read_csv(A / f"fClaHyb_{h}/fClaHyb_{h}_genome/fClaHyb_{h}_LG.chrom.sizes", sep="\t", header=None, names=["chrom", "len"])
    sz = sz[sz.chrom.str.contains("_LG")].sort_values("chrom").reset_index(drop=True)
    L = pd.read_csv(W / "insert" / f"{h}.insert_loci.tsv", sep="\t")
    cen = pd.read_csv(R / h / "centromere_check.tsv", sep="\t").set_index("chrom")
    n = len(sz)
    for i, r in sz.iterrows():
        y = n - 1 - i
        ax.add_line(plt.Line2D([0, r.len / 1e6], [y, y], color=LC.P.CHROM_LINE, lw=3, solid_capstyle="round", zorder=1))
        if r.chrom in cen.index: ax.add_patch(Rectangle((cen.loc[r.chrom, "cen_start"] / 1e6, y - 0.32), max((cen.loc[r.chrom, "cen_end"] - cen.loc[r.chrom, "cen_start"]) / 1e6, 0.25), 0.64, facecolor=LC.P.CEN, edgecolor=LC.P.CEN_EDGE, lw=0.6, zorder=2))
        l = L[L.chrom == r.chrom]; small = l[(l.qcov_bp >= 1000) & (l.qcov_bp < 3000)]; big = l[l.qcov_bp >= 3000]
        ax.scatter((small.lo + small.hi) / 2e6, [y] * len(small), marker="|", s=60, color="#7a7a7a", zorder=3, linewidths=1)
        [ax.add_patch(Rectangle((m_ - 0.4, y - 0.46), 0.8, 0.92, facecolor=CINS, edgecolor="white", lw=1.0, zorder=4)) for m_ in (big.lo + big.hi) / 2e6]
        ax.text(-0.8, y, lab(r.chrom), ha="right", va="center", fontsize=9)
        ax.text(r.len / 1e6 + 0.8, y, f"{len(big)}" if len(big) else "", va="center", fontsize=10, color=CINS, fontweight="bold")
        if len(small): ax.text(r.len / 1e6 + (3.2 if len(big) else 0.8), y, f"({len(small)})", va="center", fontsize=8, color="#777")
    ax.set_xlim(-4, 58); ax.set_ylim(-1, n); ax.set_yticks([]); ax.set_xlabel("Mb", fontsize=11)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    nb = int((L.qcov_bp >= 3000).sum()); ns = int(((L.qcov_bp >= 1000) & (L.qcov_bp < 3000)).sum())
    sp = "C. macrocephalus" if h == "Mac" else "C. gariepinus"
    ax.set_title(f"{'B' if h == 'Mac' else 'C'}  Type II insert copies, {sp}: {nb} large, {ns} partial", loc="left", fontsize=13, fontweight="bold")
ax = fig.axes[2]
ax.add_patch(Rectangle((0, -5), 0, 0, facecolor=CINS, label="≥ 3 kb of the 6-kb insert (count in slate)"))
ax.scatter([], [], marker="|", s=60, color="#7a7a7a", label="1–3 kb fragment (count in brackets)")
ax.add_patch(Rectangle((0, -5), 0, 0, facecolor=LC.P.CEN, edgecolor=LC.P.CEN_EDGE, label="centromere"))
ax.legend(loc="upper right", bbox_to_anchor=(1.0, -0.06), ncol=3, fontsize=10, frameon=False)
fig.suptitle("Composite centromeric unit: repeats in context, and where the type II insert occurs", fontsize=16, fontweight="bold", x=0.02, ha="left", y=0.995)
for ext in ("png", "pdf"): fig.savefig(W / "figure" / f"composite_unit_context.{ext}", dpi=450, bbox_inches="tight")
print("saved")
