#!/usr/bin/env python3
"""Chromosome ideograms with satellite runs, centromeres and telomeres, v4.6 (2026-10-27) - replaces share_LITE_v4.5.3/ideograms/{Gar,Mac}_ideogram_with_locations.png.
User: "update these ones as well". Changes against v4.5.3: satellite runs of record v4.6 (tier S1, unit >= 100 bp; gene arrays black, micro / minisatellites not drawn),
family colours = the shared palette (palette_v46: lineage A reds / oranges, lineage B blues, other lineage A / B light red / light blue), v4.6 family names (CgaSat / CmaSat),
centromeres from centromere_check.tsv, 45S rDNA sites, and NO text on top of the chromosomes: centromere and large-array labels sit in a fixed right-hand column, so nothing overlaps.
Out: FISH_package/ideograms/{Gar,Mac}_ideogram_with_locations_v4.6.{png,pdf}
usage: python ideogram_share_v46.py"""
import os, re, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import palette_v46 as PAL
S = "F:/data/satellite_FISH"; V = f"{S}/runs_of_record_v4.6_2026-10-27"; R0 = "F:/data/fClaHyb_func_annotation_v1.0"; OUT = f"{V}/FISH_package/ideograms"; os.makedirs(OUT, exist_ok=True)
NOM = pd.read_csv(f"{V}/satellite_family_nomenclature_v46.tsv", sep="\t").set_index("pipeline_family_id").name_v46.to_dict()
TEL, RDNA45, CEN_BOX = "#1e8449", "#6c3483", "#111111"
plt.rcParams["font.family"] = "DejaVu Sans"
for sp in ("Gar", "Mac"):
    size = {l.split()[0]: int(l.split()[1]) for l in open(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.fa.fai")}
    runs = pd.read_csv(f"{V}/{sp}/runs_of_record.tsv", sep="\t"); runs = runs[runs.tier != "EXCLUDED"]
    sat = runs[runs.unit_class == "satellite"].copy(); sat["bp"] = sat.end - sat.start0; gene = runs[runs.unit_class == "gene_array"]
    cen = pd.read_csv(f"{V}/{sp}/centromere_check.tsv", sep="\t").set_index("chrom")
    tel = pd.read_csv(f"{S}/pipeline_v4_5_fixed/out/fClaHyb_{sp}_satellites/11_telomere_blocks.tsv", sep="\t", usecols=["chrom", "start_0based", "end", "block_length_bp"]); tel = tel[tel.block_length_bp >= 500]
    hits = pd.read_csv(f"{S}/rDNA_45S/hits_{sp}.tsv", sep="\t", header=None, names="q chr pid len qs qe ss se ev bits ql".split()); hits = hits[hits.q.str.endswith("18S") & (hits.len >= 1500) & (hits.pid >= 95)]
    chroms = sorted(size, key=lambda c: int(re.search(r"(\d+)$", c).group(1))); n = len(chroms); H = 0.55; xmax = max(size.values()) / 1e6
    top = sat.groupby("family").bp.sum().sort_values(ascending=False).head(12).index.tolist()
    fig, ax = plt.subplots(figsize=(15, 0.5 * n + 3.0)); fig.patch.set_facecolor("white")
    for i, c in enumerate(chroms):
        y = n - i
        ax.add_patch(Rectangle((0, y - H / 2), size[c] / 1e6, H, fc="#f2f2f2", ec="#555555", lw=0.6, zorder=1))
        for r in sat[sat.chrom == c].itertuples(): ax.add_patch(Rectangle((r.start0 / 1e6, y - H / 2), max(r.bp / 1e6, 0.04), H, fc=PAL.family(r.family), lw=0, zorder=2))
        for r in gene[gene.chrom == c].itertuples(): ax.add_patch(Rectangle((r.start0 / 1e6, y - H / 2), max((r.end - r.start0) / 1e6, 0.04), H, fc=CEN_BOX, lw=0, zorder=2))
        for t in tel[tel.chrom == c].itertuples(): ax.add_patch(Rectangle((t.start_0based / 1e6, y - H / 2), max((t.end - t.start_0based) / 1e6, 0.06), H, fc=TEL, lw=0, zorder=3))
        for h in hits[hits.chr == c].itertuples(): ax.add_patch(Rectangle((min(h.ss, h.se) / 1e6, y - H / 2), 0.08, H, fc=RDNA45, lw=0, zorder=3))
        line2 = ""
        if c in cen.index:
            k = cen.loc[c]; a, b = k.cen_start / 1e6, k.cen_end / 1e6
            ax.add_patch(Rectangle((a, y - H / 2 - 0.08), b - a, H + 0.16, fill=False, ec="black", lw=1.4, zorder=5))
            fam = k.top_family_in_cen; line1 = f"cen {a:.2f}-{b:.2f} Mb  " + (f"{NOM.get(fam, fam)}  {k.sat_frac_in_cen:.0%} satellite" if isinstance(fam, str) and k.sat_frac_in_cen > 0 else "no satellite run in the call")
        else: line1 = "no centromere call"
        big = sat[(sat.chrom == c) & (sat.bp >= 100_000)].groupby("family").bp.sum().sort_values(ascending=False)
        if c in cen.index and isinstance(cen.loc[c, "top_family_in_cen"], str): big = big[[f for f in big.index if f != cen.loc[c, "top_family_in_cen"]]]
        if len(big): line2 = "also >= 100 kb: " + ", ".join(f"{NOM.get(f, f)} ({v / 1e3:.0f} kb)" for f, v in big.head(2).items())
        ax.text(xmax + 1.2, y + 0.08, line1, fontsize=7.5, va="center", color="#111"); ax.text(xmax + 1.2, y - 0.2, line2, fontsize=6.5, va="center", color="#555")
        ax.text(-0.6, y, c.replace("C_gar_", "").replace("C_mac_", ""), ha="right", va="center", fontsize=9)
    lin = lambda f: sat[sat.family == f].shape[0] and PAL.lineage_of(f)
    handles = [Rectangle((0, 0), 1, 1, fc=PAL.family(f)) for f in top] + [Rectangle((0, 0), 1, 1, fc=PAL.LINEAGE_COL["A"]), Rectangle((0, 0), 1, 1, fc=PAL.LINEAGE_COL["B"]), Rectangle((0, 0), 1, 1, fc=PAL.OTHER_SAT),
                                                                           Rectangle((0, 0), 1, 1, fc=CEN_BOX), Rectangle((0, 0), 1, 1, fc=TEL), Rectangle((0, 0), 1, 1, fc=RDNA45), Rectangle((0, 0), 1, 1, fill=False, ec="black", lw=1.4)]
    labels = [f"{NOM.get(f, f)}  ({f}, lineage {PAL.lineage_of(f)})" for f in top] + ["other lineage A", "other lineage B", "other satellite family", "tandem gene array (5S / tRNA / snRNA / miR-430)", "telomere (TTAGGG)n >= 500 bp", "45S rDNA (18S >= 1.5 kb, >= 95 %)", "putative centromere (synteny-inferred)"]
    ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.46, -0.035), ncol=3, fontsize=8, frameon=False, title="satellite runs of record v4.6 (unit >= 100 bp); the 12 largest families are listed", title_fontsize=8)
    ax.set_xlim(-1.5, xmax + 26); ax.set_ylim(0.3, n + 0.9); ax.set_yticks([])
    for side in ("left", "right", "top"): ax.spines[side].set_visible(False)
    ax.set_xlabel("position (Mb)")
    ax.set_title(f"{'Clarias gariepinus (Gar)' if sp == 'Gar' else 'Clarias macrocephalus (Mac)'} haplotype of fClaHyb: satellite runs, centromeres, telomeres and 45S rDNA (satellite runs of record v4.6)", fontsize=11, loc="left")
    fig.tight_layout(); fig.savefig(f"{OUT}/{sp}_ideogram_with_locations_v4.6.png", dpi=450, bbox_inches="tight"); fig.savefig(f"{OUT}/{sp}_ideogram_with_locations_v4.6.pdf", bbox_inches="tight"); plt.close(fig); print("wrote", sp)
