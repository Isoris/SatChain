"""Step 35c (2026-10-24): enlarged centromere composition, both haplotypes (assembly only, no reads / ChIP).
Per chromosome, top to bottom = region (v4.6 centromeric cluster +-150 kb). Columns per chromosome:
 strand (copy-level strand blocks, step-03 copies inside v4.6 runs; switches = orientation changes, unvalidated)
 satellite (v4.6 runs by family; full colour inside the homogenized core, pale outside)
 TE insertions (EDTA, < 50% inside satellite runs)
 homogeneity h(x) (W = 20 kb, k = 21) and core bracket (tau = 0.5; thin bracket = tau 0.3)."""
import pathlib, pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
O = pathlib.Path("F:/data/satellite_FISH/centromere_composition_2026-10-24")
R = pd.read_csv(O / "regions.tsv", sep="\t"); S = pd.read_csv(O / "satellite_runs.tsv", sep="\t"); T = pd.read_csv(O / "te_insertions.tsv", sep="\t")
H = pd.read_csv(O / "homogeneity.tsv.gz", sep="\t"); K = pd.read_csv(O / "cores.tsv", sep="\t"); GP = pd.read_csv(O / "gaps.tsv", sep="\t")
H = H[(H.W == 20000) & (H.k == 21)]
K5 = K[(K.W == 20000) & (K.k == 21) & (K.tau == 0.5)]; K3 = K[(K.W == 20000) & (K.k == 21) & (K.tau == 0.3)]
FAM = {"CmaSat001-225": ("#9e1b32", "CmaSat001-225 (MacSat03, lineage A)"), "CmaSat002-214": ("#2f5d9e", "CmaSat002-214 (lineage B, 5S-derived)"),
       "CmaSat006-212": ("#7a4fa0", "CmaSat006-212 (dimer 216+243)"), "CmaSat011-215": ("#d9822b", "CmaSat011-215 (lineage A)"),
       "CgaSat001-216": ("#3a7bbf", "CgaSat001-216 (lineage B)"), "CgaSat002-226": ("#c0392b", "CgaSat002-226 (lineage A)"),
       "CgaSat004-225": ("#e67e22", "CgaSat004-225 (lineage A)"), "CgaSat005-352": ("#8e5a2b", "CgaSat005-352 (lineage A, LG24, 352 bp)")}
LIN = {"A": "#e8a0a0", "B": "#9cbbe0", "cen_minisatellite": "#e5c84a", "unresolved_dimer_216+243": "#7a4fa0", "none": "#8a9a5b"}


def scol(r):
    if r.family in FAM:
        return FAM[r.family][0]
    if isinstance(r.gene, str) and r.gene != "none":
        return "#111111"
    return LIN.get(r.lineage, "#8a9a5b")


TEC = {"LINE": "#e6a817", "LTR": "#2e86c1", "Gypsy": "#2e86c1", "TIR": "#7d3c98", "helitron": "#1e8449"}


def tcol(t):
    for k, v in TEC.items():
        if k in t:
            return v
    return "#9e9e9e"


cops = {}
for sp in ("Mac", "Gar"):
    c = pd.read_csv(f"F:/data/satellite_FISH/pipeline_v4_5_fixed/out/fClaHyb_{sp}_satellites/03_monomer_copies.tsv", sep="\t",
                    usecols=["chrom", "genomic_start_0based", "genomic_end", "strand", "best_alignment_identity"])
    cops[sp] = c[c.best_alignment_identity >= 60]


def strand_blocks(sp, ch, lo, hi):
    s = S[(S.sp == sp) & (S.chrom == ch)]
    c = cops[sp]; c = c[(c.chrom == ch) & (c.genomic_end > lo) & (c.genomic_start_0based < hi)].sort_values("genomic_start_0based")
    ins = np.zeros(len(c), bool)
    for r in s.itertuples():
        ins |= ((c.genomic_start_0based >= r.s) & (c.genomic_end <= r.e)).values
    c = c[ins]; out = []
    for r in c.itertuples():
        if out and out[-1][2] == r.strand and r.genomic_start_0based - out[-1][1] < 2000:
            out[-1][1] = r.genomic_end
        else:
            out.append([r.genomic_start_0based, r.genomic_end, r.strand])
    plus = sum(b - a for a, b, t in out if t == "+"); minus = sum(b - a for a, b, t in out if t == "-")
    maj = "+" if plus >= minus else "-"
    return [(a, b, t == maj) for a, b, t in out]


MB = 1e6
CM = plt.get_cmap("viridis")
def Y(p, lo): return (p - lo) / MB


fig = plt.figure(figsize=(26, 22)); fig.patch.set_facecolor("white")
maxh = (R.hi - R.lo).max() / MB
for row, sp, title in ((0, "Mac", "C. macrocephalus"), (1, "Gar", "C. gariepinus")):
    ax = fig.add_axes([0.04, 0.55 - row * 0.47, 0.94, 0.40])
    rr = R[R.sp == sp].copy(); rr["n"] = rr.chrom.str.extract(r"LG(\d+)")[0].astype(int); rr = rr.sort_values("n")
    for i, r in enumerate(rr.itertuples()):
        x0 = float(i); lo, hi = r.lo, r.hi
        ax.add_patch(Rectangle((x0 + 0.28, 0), 0.22, Y(hi, lo), color="#d9d9d9", lw=0))
        cores = K5[(K5.sp == sp) & (K5.chrom == r.chrom)]
        def inside(a, b):
            return any(min(b, z) - max(a, y) > 0.5 * (b - a) for y, z in zip(cores.s, cores.e))
        for s in S[(S.sp == sp) & (S.chrom == r.chrom)].itertuples():
            a, b = max(s.s, lo), min(s.e, hi)
            if b > a:
                ax.add_patch(Rectangle((x0 + 0.28, Y(a, lo)), 0.22, Y(b, lo) - Y(a, lo), color=scol(s), lw=0))
                if isinstance(s.gene, str) and s.gene != "none" and s.family in FAM:
                    ax.add_patch(Rectangle((x0 + 0.28, Y(a, lo)), 0.22, Y(b, lo) - Y(a, lo), fill=False, hatch="///", edgecolor="white", lw=0))
        for a, b, maj in strand_blocks(sp, r.chrom, lo, hi):
            ax.add_patch(Rectangle((x0 + (0.17 if maj else 0.22), Y(a, lo)), 0.045, max(Y(b, lo) - Y(a, lo), 0.002), color="#444444" if maj else "#d62728", lw=0))
        for t in T[(T.sp == sp) & (T.chrom == r.chrom)].itertuples():
            ax.add_patch(Rectangle((x0 + 0.52, Y(t.s, lo)), 0.07, max(Y(t.e, lo) - Y(t.s, lo), 0.0015), color=tcol(t.type), lw=0))
        for g in GP[(GP.sp == sp) & (GP.chrom == r.chrom)].itertuples():
            ax.plot([x0 + 0.15, x0 + 0.62], [Y(g.s, lo)] * 2, color="black", lw=1.5, ls=":")
        h = H[(H.sp == sp) & (H.chrom == r.chrom)].sort_values("pos")
        for p_, hv in zip(h.pos, h.h):
            if hv == hv:
                ax.add_patch(Rectangle((x0 + 0.62, Y(p_ - 2500, lo)), 0.13, 0.005, color=CM(hv), lw=0))
        for kk, lw, off in ((K3[(K3.sp == sp) & (K3.chrom == r.chrom)], 1.0, 0.79), (cores, 2.5, 0.83)):
            for c in kk.itertuples():
                ax.plot([x0 + off] * 2, [Y(c.s, lo), Y(c.e, lo)], color="black", lw=lw, solid_capstyle="butt")
        ax.text(x0 + 0.4, -0.03, str(r.n), ha="center", va="bottom", fontsize=13, fontweight="bold")
    ax.set_xlim(-0.1, 28.2); ax.set_ylim(maxh + 0.05, -0.1); ax.set_xticks([])
    ax.set_ylabel("Mb from region start", fontsize=12); ax.tick_params(axis="y", labelsize=11)
    for s_ in ("top", "right", "bottom"):
        ax.spines[s_].set_visible(False)
    ax.text(-0.03, 1.04, "AB"[row], transform=ax.transAxes, fontsize=20, fontweight="bold")
    ax.text(0.0, 1.04, title, transform=ax.transAxes, fontsize=16, style="italic", fontweight="bold")
    fams = [f for f in FAM if f.startswith("Cma" if sp == "Mac" else "Cga")]
    h1 = [Patch(color=FAM[f][0], label=FAM[f][1]) for f in fams] + [Patch(color=LIN["A"], label="other lineage A"), Patch(color=LIN["B"], label="other lineage B"),
          Patch(color=LIN["cen_minisatellite"], label="centromeric minisatellite"), Patch(color="#111111", label="tandem gene array (5S / tRNA)"),
          Patch(color="#d9d9d9", label="non-satellite"), Patch(facecolor="#777777", hatch="///", edgecolor="white", label="family also flagged as 5S gene array")]
    ax.legend(handles=h1, loc="lower right", fontsize=11, ncol=2, frameon=False)
h2 = [Patch(color="#444444", label="copies, majority strand"), Patch(color="#d62728", label="copies, minority strand (orientation switch, unvalidated)"),
      Patch(color="#e6a817", label="TE: LINE"), Patch(color="#2e86c1", label="TE: LTR"), Patch(color="#7d3c98", label="TE: DNA/TIR"), Patch(color="#1e8449", label="TE: Helitron"),
            Line2D([], [], color="black", lw=2.5, label="homogenized core, h >= 0.5 over >= 50 kb"),
      Line2D([], [], color="black", lw=1, label="core at h >= 0.3"), Line2D([], [], color="black", lw=1.5, ls=":", label="assembly gap")]
cax = fig.add_axes([0.80, 0.025, 0.15, 0.012]); cb = fig.colorbar(plt.cm.ScalarMappable(cmap=CM), cax=cax, orientation="horizontal")
cb.set_label("homogeneity h: k-mer (k = 21) Jaccard of adjacent 20-kb windows", fontsize=11)
fig.legend(handles=h2, loc="lower center", ncol=4, fontsize=12, frameon=False, bbox_to_anchor=(0.42, 0.0))
fig.suptitle("Centromere composition from the assembly alone: satellite families, strand, TE insertions and the homogenized core",
             fontsize=20, fontweight="bold", x=0.03, ha="left", y=0.995)
for ext in ("png", "pdf"):
    fig.savefig(O / f"centromere_composition.{ext}", dpi=130)
print("ok")
