"""Genome-wide ideogram of the step-02 intervals and their problems (2026-10-27).
Per species, every chromosome to scale. Ticks below the bar: multi-run step-02 intervals (SAT/NEW ids) coloured by what the run-based definition
(compare_array_methods.py, G = 20 kb) does to them: identical, split because it mixes families, split at an assembly gap, split only by an
interruption > 20 kb, merged with a neighbour. Black lines on the bar: assembly gaps. Above the bar: the major centromere (largest
centromeric interval), labelled with its family and spacer type (RTE / no-RTE / few interruptions).
usage: python -P ideogram_intervals.py"""
import os, pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "F:/data/satellite_FISH/runs_of_record_v4.6_2026-10-27/code"); import palette_v46 as PAL
from matplotlib.patches import Rectangle, Patch
B = os.path.dirname(os.path.abspath(__file__)); R = "F:/data/fClaHyb_func_annotation_v1.0"
D0 = "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/figure/diagnostics/centromeres"
CEN = pd.read_csv(f"{D0}/index.tsv", sep="\t").merge(pd.read_csv(f"{D0}/rte_in_interruptions.tsv", sep="\t")[["array_id", "interruptions", "share"]], on="array_id")
STR = pd.read_csv(f"{D0}/array_strand.tsv", sep="\t")[["array_id", "strand_vs_family"]]; CEN = CEN.merge(STR, on="array_id", how="left")
KCOL = {"1:1": "#bdbdbd", "family": "#3d4f66", "gap": "#7a7467", "interruption": "#b5468a", "merged": "#2a9d8f", "complex": "#2a9d8f"}   # call types: no family / class hue (reds, oranges, blues, purples) is reused
TRIE = {"RTE spacer": "#7a5200", "no-RTE spacer": PAL.CEN_EDGE, "few interruptions": "#5b5b5b"}
TRI = {"RTE spacer": PAL.TE["LINE"], "no-RTE spacer": "white", "few interruptions": "#9a9a9a"}   # major-centromere marker: amber = RTE (LINE) spacer, open = no RTE, grey = few interruptions
C = pd.read_csv(f"{R}/fClaHyb_Mac/fClaHyb_Mac_satelliteDNA/fClaHyb_Mac.satDNA_v46.catalogue.tsv", sep="\t")
for sp in ("Mac", "Gar"):
    A = pd.read_csv(f"{B}/compare_{sp}_20k.tsv", sep="\t"); A = A[A.runs >= 2]
    cat = pd.read_csv(f"{R}/fClaHyb_{sp}/fClaHyb_{sp}_satelliteDNA/fClaHyb_{sp}.satDNA_v46.catalogue.tsv", sep="\t")
    def why(r):
        if r.kind == "1:1": return "1:1"
        if r.kind in ("A merged by B", "complex"): return "merged"
        if r.n_fam > 1: return "family"
        if r.gap_inside: return "gap"
        return "interruption"
    A["why"] = A.apply(why, axis=1)
    L = {l.split()[0]: int(l.split()[1]) for l in open(f"{R}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.fa.fai")}
    chroms = sorted(L, key=lambda c: int(c.split("LG")[-1]))
    G = pd.read_csv(f"{R}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.gaps.detgaps.bed", sep="\t", header=None, usecols=[0, 1, 2], names=["c", "a", "b"])
    mx = max(L.values()) / 1e6; n = len(chroms)
    fig, ax = plt.subplots(figsize=(18, 0.62 * n + 2)); fig.patch.set_facecolor("white")
    for i, c in enumerate(chroms):
        y = n - 1 - i
        ax.add_patch(Rectangle((0, y - 0.12), L[c] / 1e6, 0.24, color="#eceae4", ec="#9a9890", lw=0.6))
        for g in G[(G.c == c) & (G.a > 0)].itertuples(): ax.plot([g.a / 1e6] * 2, [y - 0.14, y + 0.14], color="black", lw=0.8)
        for r in A[A.chrom == c].itertuples():
            w = max((r.hi - r.lo) / 1e6, 0.06); ax.add_patch(Rectangle((r.lo / 1e6, y - 0.36), w, 0.2, color=KCOL[r.why], lw=0))
        k = CEN[CEN.chrom == c]
        if len(k):
            k = k.iloc[0]; mid = (cat[cat.array_id == k.array_id].start0.min() + cat[cat.array_id == k.array_id].end.max()) / 2e6
            spc = "few interruptions" if k.interruptions < 3 else ("RTE spacer" if k.share >= 0.5 else "no-RTE spacer")
            ax.plot([mid], [y + 0.22], marker="v", ls="", mfc=TRI[spc], mec=TRIE[spc], mew=1.0, ms=9)
            ax.text(mid + 0.4, y + 0.22, f"{k.family} ({spc}; {k.strand_vs_family})", fontsize=7.5, va="center", color=PAL.INK2)
        ax.text(-0.6, y, c.split("_")[-1], ha="right", va="center", fontsize=9)
    ax.set_xlim(-2.5, mx + 6); ax.set_ylim(-0.7, n - 0.3); ax.set_yticks([]); ax.set_xlabel("Mb")
    for s_ in ("left", "right", "top"): ax.spines[s_].set_visible(False)
    cnt = A.why.value_counts().to_dict()
    lg = [Patch(color=KCOL[k], label=f"{lab} ({cnt.get(k, 0)})") for k, lab in (("1:1", "interval = one run-based array"), ("family", "mixes satellite families"),
          ("gap", "crosses an assembly gap"), ("interruption", "split only by an interruption > 20 kb"), ("merged", "merged / complex"))]
    lg += [plt.Line2D([], [], marker="v", ls="", mfc=TRI[k_], mec=TRIE[k_], mew=1.0, ms=9, label=l_) for k_, l_ in (("RTE spacer", "major centromere, RTE spacer"), ("no-RTE spacer", "major centromere, no-RTE spacer"),
           ("few interruptions", "major centromere, few interruptions"))] + [plt.Line2D([], [], color="black", lw=1, label="assembly gap")]
    ax.legend(handles=lg, loc="lower right", fontsize=9, frameon=False)
    ax.set_title(f"fClaHyb {sp}: multi-run step-02 intervals (SAT/NEW ids, n = {len(A)}) and what the run-based array definition does to them;" + chr(10) + "major centromeres with family, "
                 f"spacer type and strand vs family consensus", loc="left", fontsize=12, fontweight="bold")
    fig.tight_layout(); fig.savefig(f"{B}/ideogram_intervals_{sp}.png", dpi=450); fig.savefig(f"{B}/ideogram_intervals_{sp}.png".replace(".png", ".pdf")); print("wrote", sp, cnt)