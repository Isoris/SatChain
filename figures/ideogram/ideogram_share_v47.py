#!/usr/bin/env python3
"""Chromosome ideograms with satellite runs, candidate centromeres and telomeres, v4.7 (2026-11-01; from ideogram_share_v46.py) - replaces share_LITE_v4.5.3/ideograms/{Gar,Mac}_ideogram_with_locations.png.
User: "update these ones as well". Changes against v4.5.3: satellite runs of record (tier S1, unit >= 100 bp; gene arrays black, micro / minisatellites not drawn),
family colours = the shared palette (palette_v47: cluster C2.1 reds / oranges, cluster C2.2 blues, other families of C2.1 / C2.2 light red / light blue), v4.6 family names (CgaSat / CmaSat),
candidate centromeres from centromere_check.tsv, 45S rDNA sites, and NO text on top of the chromosomes: candidate-centromere and large-array labels sit in a fixed right-hand column, so nothing overlaps.
Numbered-story layout (2026-11-01, user: "it is the best for understanding, it's for noobs"): two banners above the ideogram, "1  Satellite runs of record ..." (the chromosomes) and
"2  Candidate centromeres ..." (the right-hand column), in the banner contract of figures/paper_figures/code/assemble_v2.py (BAN 7 mm, bold number + bold sentence, 0.5-pt navy hairline,
house style: Arial, navy #17252D, no top / right spines).  The former title is kept as the technical line below the banners.
Out: <run>/ideograms/{Gar,Mac}_ideogram_with_locations_v4.7.{png,pdf}
usage: python ideogram_share_v47.py"""
import os, re, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import palette_v47 as PAL
S = "F:/data/satellite_FISH"; RUN = os.environ.get("SATCHAIN_RUN", f"{S}/satchain_v4.7_2026-11-01"); V = f"{RUN}/runs_of_record"; R0 = "F:/data/fClaHyb_func_annotation_v1.0"; OUT = f"{RUN}/ideograms"; os.makedirs(OUT, exist_ok=True)
NOMT = pd.read_csv(f"{V}/satellite_family_nomenclature.tsv", sep="\t")                      # names are per haplotype (MacSat03 runs exist in both)
TEL, RDNA45, CEN_BOX = "#1e8449", "#6c3483", "#111111"
NAVY = "#17252D"; MM = 1 / 25.4; BAN = 7.0                                                  # banner contract of assemble_v2.py (height in mm)
plt.rcParams.update({"font.family": "Arial", "text.color": NAVY, "axes.labelcolor": NAVY, "xtick.color": NAVY, "ytick.color": NAVY, "axes.edgecolor": NAVY,
                     "axes.linewidth": 0.5, "xtick.major.width": 0.5, "pdf.fonttype": 42})
BANNERS = [("1", "Satellite runs of record (unit >= 100 bp) along the chromosomes; gene arrays in black"),
           ("2", "Candidate centromeres: the synteny-inferred interval and the satellite share inside it")]


def banner(fig, W, H, top, num, text):
    """assemble_v2.banner: hairline at `top` mm from the top of the figure, bold number at 0.5 mm, bold sentence at 5 mm, 8 pt navy"""
    fig.add_artist(Line2D([0, 1], [1 - top / H, 1 - top / H], color=NAVY, lw=0.5, transform=fig.transFigure))
    fig.text(0.5 / W, 1 - (top + 1.0) / H, num, fontsize=8, fontweight="bold", ha="left", va="top", color=NAVY)
    fig.text(5.0 / W, 1 - (top + 1.0) / H, text, fontsize=8, fontweight="bold", ha="left", va="top", color=NAVY)


for sp in ("Gar", "Mac"):
    NOM = NOMT[NOMT.haplotype == sp].set_index("family").name.to_dict()
    size = {l.split()[0]: int(l.split()[1]) for l in open(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.fa.fai")}
    runs = pd.read_csv(f"{V}/{sp}/runs_of_record.tsv", sep="\t"); runs = runs[runs.tier != "EXCLUDED"]
    sat = runs[runs.unit_class == "satellite"].copy(); sat["bp"] = sat.end - sat.start0; gene = runs[runs.unit_class == "gene_array"]
    cen = pd.read_csv(f"{V}/{sp}/centromere_check.tsv", sep="\t").set_index("chrom")
    tel = pd.read_csv(f"{S}/pipeline_v4_5_fixed/out/fClaHyb_{sp}_satellites/11_telomere_blocks.tsv", sep="\t", usecols=["chrom", "start_0based", "end", "block_length_bp"]); tel = tel[tel.block_length_bp >= 500]
    hits = pd.read_csv(f"{S}/rDNA_45S/hits_{sp}.tsv", sep="\t", header=None, names="q chr pid len qs qe ss se ev bits ql".split()); hits = hits[hits.q.str.endswith("18S") & (hits.len >= 1500) & (hits.pid >= 95)]
    chroms = sorted(size, key=lambda c: int(re.search(r"(\d+)$", c).group(1))); n = len(chroms); H = 0.55; xmax = max(size.values()) / 1e6
    top = sat.groupby("family").bp.sum().sort_values(ascending=False).head(12).index.tolist()
    fig, ax = plt.subplots(figsize=(15, 0.5 * n + 3.6)); fig.patch.set_facecolor("white")       # +0.6 in = the two banner rows (2 x 7 mm) and the technical line
    W_mm, H_mm = fig.get_size_inches() / MM
    for i, (num, text) in enumerate(BANNERS): banner(fig, W_mm, H_mm, i * BAN, num, text)
    fig.add_artist(Line2D([0, 1], [1 - 2 * BAN / H_mm] * 2, color=NAVY, lw=0.5, transform=fig.transFigure))
    fig.subplots_adjust(left=0.03, right=0.995, top=1 - (2 * BAN + 9.0) / H_mm, bottom=0.085)  # banners + technical line reserved at the top, legend at the bottom
    for i, c in enumerate(chroms):
        y = n - i
        ax.add_patch(Rectangle((0, y - H / 2), size[c] / 1e6, H, fc="#f2f2f2", ec="#555555", lw=0.6, zorder=1))
        for r in sat[sat.chrom == c].itertuples(): ax.add_patch(Rectangle((r.start0 / 1e6, y - H / 2), max(r.bp / 1e6, 0.04), H, fc=PAL.family(r.family), lw=0, zorder=2))
        for r in gene[gene.chrom == c].itertuples(): ax.add_patch(Rectangle((r.start0 / 1e6, y - H / 2), max((r.end - r.start0) / 1e6, 0.04), H, fc=CEN_BOX, lw=0, zorder=2))
        for t in tel[tel.chrom == c].itertuples(): ax.add_patch(Rectangle((t.start_0based / 1e6, y - H / 2), max((t.end - t.start_0based) / 1e6, 0.06), H, fc=TEL, lw=0, zorder=3))
        for h in hits[hits.chr == c].itertuples(): ax.add_patch(Rectangle((min(h.ss, h.se) / 1e6, y - H / 2), 0.08, H, fc=RDNA45, lw=0, zorder=3))
        line2 = ""
        if c in cen.index:
            k = cen.loc[c]; a, b = k.cen_start0 / 1e6, k.cen_end / 1e6
            ax.add_patch(Rectangle((a, y - H / 2 - 0.08), b - a, H + 0.16, fill=False, ec="black", lw=1.4, zorder=5))
            fam = k.dominant_family; line1 = f"candidate cen {a:.2f}-{b:.2f} Mb  " + (f"{NOM.get(fam, fam)}  {k.support:.0%} satellite" if isinstance(fam, str) and fam != "none" and k.support > 0 else "no satellite run in the call")
        else: line1 = "no candidate centromere call"
        big = sat[(sat.chrom == c) & (sat.bp >= 100_000)].groupby("family").bp.sum().sort_values(ascending=False)
        if c in cen.index and isinstance(cen.loc[c, "dominant_family"], str): big = big[[f for f in big.index if f != cen.loc[c, "dominant_family"]]]
        if len(big): line2 = "also >= 100 kb: " + ", ".join(f"{NOM.get(f, f)} ({v / 1e3:.0f} kb)" for f, v in big.head(2).items())
        ax.text(xmax + 1.2, y + 0.08, line1, fontsize=7.5, va="center", color="#111"); ax.text(xmax + 1.2, y - 0.2, line2, fontsize=6.5, va="center", color="#555")
        ax.text(-0.6, y, c.replace("C_gar_", "").replace("C_mac_", ""), ha="right", va="center", fontsize=9)
    handles = [Rectangle((0, 0), 1, 1, fc=PAL.family(f)) for f in top] + [Rectangle((0, 0), 1, 1, fc=PAL.LINEAGE_COL["A"]), Rectangle((0, 0), 1, 1, fc=PAL.LINEAGE_COL["B"]), Rectangle((0, 0), 1, 1, fc=PAL.OTHER_SAT),
                                                                           Rectangle((0, 0), 1, 1, fc=CEN_BOX), Rectangle((0, 0), 1, 1, fc=TEL), Rectangle((0, 0), 1, 1, fc=RDNA45), Rectangle((0, 0), 1, 1, fill=False, ec="black", lw=1.4)]
    labels = [f"{NOM.get(f, f)}  ({f}, {PAL.cluster_of(f)})" for f in top] + ["other family of " + PAL.GROUP_LABEL["A"], "other family of " + PAL.GROUP_LABEL["B"], "other satellite family (no cluster)", "tandem gene array (5S / tRNA / snRNA / miR-430)", "telomere (TTAGGG)n >= 500 bp", "45S rDNA (18S >= 1.5 kb, >= 95 %)", "candidate centromere (synteny-inferred interval)"]
    ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.46, -0.035), ncol=3, fontsize=8, frameon=False, title="satellite runs of record v4.7 (unit >= 100 bp; gene arrays apart); the 12 largest families are listed; clusters = k-mer composition clusters of step 55 (a grouping, not a history)", title_fontsize=8)
    ax.set_xlim(-1.5, xmax + 26); ax.set_ylim(0.3, n + 0.9); ax.set_yticks([])
    for side in ("left", "right", "top"): ax.spines[side].set_visible(False)
    ax.set_xlabel("position (Mb)")
    ax.set_title(f"{'Clarias gariepinus (Gar)' if sp == 'Gar' else 'Clarias macrocephalus (Mac)'} haplotype of fClaHyb: satellite runs, candidate centromeres, telomeres and 45S rDNA (satellite runs of record v4.7)", fontsize=9, loc="left", pad=6)
    fig.savefig(f"{OUT}/{sp}_ideogram_with_locations_v4.7.png", dpi=450, bbox_inches="tight"); fig.savefig(f"{OUT}/{sp}_ideogram_with_locations_v4.7.pdf", bbox_inches="tight"); plt.close(fig); print("wrote", sp)
