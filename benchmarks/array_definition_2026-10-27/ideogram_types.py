"""Array types from sequence-level metrics, coloured on the ideogram (2026-10-27). No copies needed (works for NEW arrays too).
Per multi-run step-02 interval: relative span = run bp / interval span; interruptions per 100 kb = (runs - 1) / span; unit_bp (catalogue
mode); k-mer homogeneity = top-unit coverage (canonical 21-mers of the runs, the unit_bp most frequent k-mers' share); gene-array flag
(name_v46 'gene_array:'). Type rules (relative-span cut R swept 0.7 / 0.8 / 0.9; R = 0.8 drawn):
  gene array                     gene_array flag
  micro / mini-satellite         unit <= 60 bp
  colinear, homogeneous          relative span >= R and top-unit coverage >= 0.5
  colinear, multi-variant        relative span >= R and top-unit coverage < 0.5
  composite periodic             0.3 <= relative span < R and >= 5 interruptions per 100 kb
  interrupted, irregular         0.3 <= relative span < R, fewer interruptions
  scattered                      relative span < 0.3
Outputs: array_types_<sp>.tsv, ideogram_types_<sp>.png, counts per R printed.
usage: python -P ideogram_types.py"""
import os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
B = os.path.dirname(os.path.abspath(__file__)); R0 = "F:/data/fClaHyb_func_annotation_v1.0"
sys.path.insert(0, B); from kmer_homogeneity import Fasta, kcodes                      # reuse (module guarded below)
TCOL = {"gene array": "#7f7f7f", "micro / mini-satellite": "#bcbd22", "colinear, homogeneous": "#1f4e9c", "colinear, multi-variant": "#6aa1e0",
        "composite periodic": "#c0392b", "interrupted, irregular": "#e67e22", "scattered": "#b39ddb"}
def typ(r, R):
    if r.gene: return "gene array"
    if r.unit_bp <= 60: return "micro / mini-satellite"
    if r.rel_span >= R: return "colinear, homogeneous" if r.top_unit_cov >= 0.5 else "colinear, multi-variant"
    if r.rel_span >= 0.3: return "composite periodic" if r.int_per_100kb >= 5 else "interrupted, irregular"
    return "scattered"
for sp in ("Mac", "Gar"):
    C = pd.read_csv(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_satelliteDNA/fClaHyb_{sp}.satDNA_v46.catalogue.tsv", sep="\t")
    fa = Fasta(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.fa"); rows = []
    for aid, g in C.groupby("array_id"):
        if len(g) < 2: continue
        lo, hi = int(g.start0.min()), int(g.end.max()); L = max(1, int(g.unit_bp.mode().iloc[0]))
        seq = b"N".join(fa.get(r.chrom, int(r.start0), int(r.end)) for r in g.itertuples())
        _, cnt = np.unique(kcodes(seq, 21), return_counts=True); cnt = np.sort(cnt)[::-1]
        cov = float(cnt[:L].sum() / cnt.sum()) if len(cnt) else np.nan
        rows.append(dict(array_id=aid, chrom=g.chrom.iloc[0], lo=lo, hi=hi, runs=len(g), sat_kb=round(g.bp.sum() / 1e3, 1), span_kb=round((hi - lo) / 1e3, 1),
                         rel_span=round(g.bp.sum() / max(1, hi - lo), 3), int_per_100kb=round((len(g) - 1) / max(1, hi - lo) * 1e5, 2), unit_bp=L,
                         top_unit_cov=round(cov, 3), gene=str(g.name_v46.iloc[0]).startswith("gene_array"), family=g.family_v46.iloc[0],
                         centromeric=(g.putative_centromeric_array_L1 == "yes").any()))
    T = pd.DataFrame(rows)
    for R in (0.7, 0.8, 0.9): T[f"type_R{int(R * 10)}"] = [typ(r, R) for r in T.itertuples()]
    T.to_csv(f"{B}/array_types_{sp}.tsv", sep="\t", index=False)
    print(f"== {sp}: {len(T)} multi-run intervals")
    print(pd.DataFrame({f"R = {R}": T[f"type_R{int(R * 10)}"].value_counts() for R in (0.7, 0.8, 0.9)}).fillna(0).astype(int).to_string())
    print("   centromeric intervals by type (R = 0.8):", T[T.centromeric]["type_R8"].value_counts().to_dict())
    Lc = {l.split()[0]: int(l.split()[1]) for l in open(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.fa.fai")}
    chroms = sorted(Lc, key=lambda c: int(c.split("LG")[-1])); n = len(chroms)
    G = pd.read_csv(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.gaps.detgaps.bed", sep="\t", header=None, usecols=[0, 1], names=["c", "a"])
    fig, ax = plt.subplots(figsize=(18, 0.55 * n + 2)); fig.patch.set_facecolor("white")
    for i, c in enumerate(chroms):
        y = n - 1 - i; ax.add_patch(Rectangle((0, y - 0.1), Lc[c] / 1e6, 0.2, color="#eceae4", ec="#9a9890", lw=0.6))
        for g in G[(G.c == c) & (G.a > 0)].itertuples(): ax.plot([g.a / 1e6] * 2, [y - 0.12, y + 0.12], color="black", lw=0.7)
        for r in T[T.chrom == c].itertuples():
            ax.add_patch(Rectangle((r.lo / 1e6, y - 0.36), max((r.hi - r.lo) / 1e6, 0.07), 0.24, color=TCOL[r.type_R8], lw=0))
            if r.centromeric: ax.plot([(r.lo + r.hi) / 2e6], [y + 0.22], marker="v", ms=5, color="black")
        ax.text(-0.6, y, c.split("_")[-1], ha="right", va="center", fontsize=9)
    ax.set_xlim(-2.5, max(Lc.values()) / 1e6 + 1); ax.set_ylim(-0.7, n - 0.3); ax.set_yticks([]); ax.set_xlabel("Mb")
    for s_ in ("left", "right", "top"): ax.spines[s_].set_visible(False)
    vc = T["type_R8"].value_counts()
    ax.legend(handles=[Patch(color=v, label=f"{k} ({vc.get(k, 0)})") for k, v in TCOL.items()] + [plt.Line2D([], [], marker="v", ls="", color="black", label="centromeric (L1)")],
              loc="lower right", fontsize=9, frameon=False)
    ax.set_title(f"fClaHyb {sp}: multi-run satellite intervals (n = {len(T)}) coloured by array type (relative span = satellite / span, cut 0.8; k-mer homogeneity; unit; gene arrays)",
                 loc="left", fontsize=12, fontweight="bold")
    fig.tight_layout(); fig.savefig(f"{B}/ideogram_types_{sp}.png", dpi=450); plt.close(fig)