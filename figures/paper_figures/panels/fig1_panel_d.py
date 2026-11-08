"""Figure 1d: organisation of the major centromeres to scale in both species (92 x 84 mm).  Redraw of diagnostics/centromeres/centromeres_rlvl1_overview_{Mac,Gar}.png (made by diagnostics/centromere_levels.py --rlvl1-only)
in the house style, from the same tables, one shared kb scale for both species.
Data:
  F:/data/satellite_FISH/satchain_v4.7_2026-11-01/cat59/59_catalogue_{Mac,Gar}.tsv                                        SatChain v4.7 catalogue = runs of record (runs of each array; matched on array_id_v46)
  F:/data/satellite_FISH/satchain_v4.7_2026-11-01/runs_of_record/satellite_family_nomenclature.tsv                        family -> name, element_class (gene arrays), seq_cluster
  F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/figure/diagnostics/centromeres/index.tsv                              the 27 + 27 major centromere arrays (v4.6 array_id per chromosome)
  F:/data/SatChain/benchmarks/array_definition_2026-10-27/array_types_rlvl2_{Mac,Gar}.tsv                              relative span, interruptions / 100 kb, top-unit coverage (array kind)
  F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/split/rte_genome_hits_both.tsv                                          RTE LINE genome hits (gold)
  F:/data/fClaHyb_func_annotation_v1.0/fClaHyb_{Mac,Gar}/..._genome/*_LG.gaps.detgaps.bed                                assembly gaps (hatched)
Array kind (rules of ideogram_types.py): gene array if the array is a gene array (index name or v4.7 unit_class gene_array); single run if 1 run; colinear (C homogeneous / V multi-variant) if relative span >= R;
composite periodic (P) if 0.3 <= span < R and >= 5 interruptions / 100 kb; interrupted irregular (I) if 0.3 <= span < R with fewer; scattered (S) < 0.3.  R = 0.80 drawn; swept 0.80 / 0.85 / 0.90
(fig1_panel_d_numbers.tsv, columns kind_R080 / kind_R085 / kind_R090).  Arrays whose chromosome has no row in the v4.7 tables (Mac LG14 / 17 / 18, Gar LG23 / 27 / 28) are shown as empty rows ("no v4.7 runs").
usage: python fig1_panel_d.py"""
import sys, pathlib
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "code")); import satellite_style as ST
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

D = "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/figure/diagnostics/centromeres"; BM = "F:/data/SatChain/benchmarks/array_definition_2026-10-27"; R0 = "F:/data/fClaHyb_func_annotation_v1.0"
RTE_HITS = "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/split/rte_genome_hits_both.tsv"
GENE = "#e6a100"; GENE_PALE = "#f4e2a8"; BAND = "#e3e0d7"; GAPFILL = "#dcd8cc"; GAPHATCH = "#8A949C"
CODE = {"composite periodic": "P", "interrupted, irregular": "I", "colinear, homogeneous": "C", "colinear, multi-variant": "V", "scattered": "S", "single run": "R", "gene array": "G", "absent": ""}
RS = (0.80, 0.85, 0.90)


def kind(r, R):
    if r["runs"] == 0: return "absent"
    if r["gene"]: return "gene array"
    if r["runs"] == 1: return "single run"
    if r["rel_span"] >= R: return "colinear, homogeneous" if r["top_unit_cov"] >= 0.5 else "colinear, multi-variant"
    if r["rel_span"] >= 0.3: return "composite periodic" if r["int_per_100kb"] >= 5 else "interrupted, irregular"
    return "scattered"


def load():
    I = pd.read_csv(f"{D}/index.tsv", sep="\t")
    H = pd.read_csv(RTE_HITS, sep="\t"); H["a"] = H[["sstart", "send"]].min(axis=1); H["b"] = H[["sstart", "send"]].max(axis=1)
    rows = []
    for sp in ("Mac", "Gar"):
        C = ST.catalogue(sp)
        TY = pd.read_csv(f"{BM}/array_types_rlvl2_{sp}.tsv", sep="\t").set_index("array_id")
        G = pd.read_csv(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.gaps.detgaps.bed", sep="\t", header=None, usecols=[0, 1, 2], names=["c", "a", "b"])
        for r in I[I.species == sp].itertuples():
            x = C[C.array_id_v46 == r.array_id].sort_values("start0")
            absent = len(x) == 0
            lo, hi = (int(x.start0.min()), int(x.end.max())) if not absent else (0, 0)
            gene = str(r.name).startswith("gene_array") or (not absent and x.element_class.astype(str).str.startswith("gene_array").any())   # cat59: element_class carries the gene-array call
            t = TY.loc[r.array_id] if r.array_id in TY.index else None
            fam = str(x.groupby("family").bp.sum().idxmax()) if not absent else ST.pipeline_name(r.name)      # dominant family of the array by bp in the v4.7 rows
            d = dict(species=sp, chrom=r.chrom, lg=int(r.chrom.split("LG")[-1]), array_id_v46=r.array_id, array_id_v47=(str(x.array_id.iloc[0]) if not absent else ""), family=fam,
                     manuscript_name=("gene array" if gene else fam), gene=gene, runs=len(x), start=lo, end=hi, span_kb=round((hi - lo) / 1e3, 1),
                     sat_kb=round(x.bp.sum() / 1e3, 1), sat_pct=(round(100 * x.bp.sum() / (hi - lo)) if hi > lo else 0), rel_span=(t.rel_span if t is not None else 1.0), int_per_100kb=(t.int_per_100kb if t is not None else 0.0),
                     top_unit_cov=(t.top_unit_cov if t is not None else 1.0), table_span_kb=(t.span_kb if t is not None else np.nan), table_type_R8=(t.type_R8 if t is not None else ""), index_span_kb=r.span_kb, index_runs=r.runs,
                     seq_cluster_v47=(str(x.seq_cluster.iloc[0]) if not absent else ST.family_cluster(fam)), element_class_v47=(str(x.element_class.iloc[0]) if not absent else ""),
                     in_v47_tables=int(not absent), rte_hits=int(((H.chrom == r.chrom) & (H.b > lo) & (H.a < hi)).sum()), assembly_gaps=int(((G.c == r.chrom) & (G.b > lo) & (G.a < hi)).sum()))
            for R in RS: d[f"kind_R{int(R * 100):03d}"] = kind(d, R)
            rows.append(d); d["_runs"] = list(zip(x.start0, x.end))
            h = H[(H.chrom == r.chrom) & (H.b > lo) & (H.a < hi)]; d["_rte"] = list(zip(h.a, h.b))
            g = G[(G.c == r.chrom) & (G.b > lo) & (G.a < hi)]; d["_gaps"] = list(zip(g.a, g.b))
    return pd.DataFrame(rows)


def draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None):
    T = load()
    if letter: ST.letter(fig, left_mm, top_mm, letter)
    _W, _H = fig.get_size_inches() / ST.MM                                                    # small bold header, as in the other panels: what is shown
    fig.text((left_mm + 4.5) / _W, 1 - (top_mm + 0.4) / _H, "Major satellite arrays to scale, both species", fontsize=7, fontweight="bold", ha="left", va="top", color=ST.NAVY)
    W_, H_ = fig.get_size_inches() / ST.MM
    NAMEH = 4.0; PITCH = 2.24; LAB = 4.2; PW = (w_mm - 2.0) / 2 - LAB; XM = {"Mac": 700.0, "Gar": 1330.0}          # kb; one scale PER species (runs of the short Mac arrays stay separable); Mac LG20 (1.2 Mb) is cut at the axis end with a break mark
    col0 = {"Mac": 0.0, "Gar": (w_mm - 2.0) / 2 + 2.0}; y0 = top_mm + 4 + NAMEH
    axes = {}; kx = {}
    _m = T[(T.species == "Mac") & (~T.gene.astype(bool)) & (T.runs > 0)]; MAC_DOM = _m.manuscript_name.value_counts().index[0]
    _dm = _m[_m.manuscript_name == MAC_DOM]; MAC_KEYROW = _dm.loc[_dm.end.sub(_dm.start).idxmin(), "lg"]
    for sp, title in (("Mac", "C. macrocephalus"), ("Gar", "C. gariepinus")):
        XMAX = XM[sp]; d = T[T.species == sp].sort_values("lg").reset_index(drop=True); n = len(d); hh = n * PITCH
        ax = ST.ax_mm(fig, left_mm + col0[sp] + LAB, y0, PW, hh); axes[sp] = ax
        ax.set_xlim(0, XMAX); ax.set_ylim(n - 0.5, -0.5); ax.set_facecolor("none")
        for s_ in ("top", "right", "left"): ax.spines[s_].set_visible(False)
        ax.spines["bottom"].set_position(("outward", 1.5)); ax.spines["bottom"].set_bounds(0, 600 if sp == "Mac" else 1200)
        ax.set_yticks([]); ax.set_xticks([0, 200, 400, 600] if sp == "Mac" else [0, 400, 800, 1200]); ax.tick_params(axis="x", labelsize=7, length=2.0, pad=1.5)
        ax.set_xlabel("kb from the first run", fontsize=7, labelpad=1.0); ax.xaxis.set_label_coords(0.45, -(4.9 / hh))
        ax.annotate(title, xy=(-LAB / PW * XMAX, -0.5), xytext=(0, 3), textcoords="offset points", fontsize=7, fontweight="bold", fontstyle="italic", ha="left", va="bottom", color=ST.NAVY, annotation_clip=False)
        for i, r in d.iterrows():
            ax.text(-0.8 / PW * XMAX, i, str(r.lg), fontsize=7, ha="right", va="center", color=ST.NAVY if r.runs else ST.MUTED, clip_on=False)
            if r.runs == 0:                                                                # chromosome absent from the v4.7 tables
                ax.text(8.0 * XMAX / 1330.0, i, "no v4.7 rows", fontsize=6, ha="left", va="center", color=ST.MUTED, zorder=5); continue
            lo = r.start; L = (r.end - lo) / 1e3
            ax.add_patch(Rectangle((0, i - 0.17), L, 0.34, facecolor=BAND, edgecolor="none", lw=0, zorder=1))
            for a, b in r._runs:
                x0, w = (a - lo) / 1e3, (b - a) / 1e3
                if r.gene: ax.add_patch(Rectangle((x0, i - 0.36), w, 0.72, facecolor=GENE_PALE, edgecolor=GENE, hatch="////", lw=0, zorder=2))
                else: ST.block(ax, x0, i - 0.36, w, 0.72, ST.family_colour(r.family), zorder=2)
            for a, b in r._rte: ax.add_patch(Rectangle(((a - lo) / 1e3 - 1.0, i - 0.50), (b - a) / 1e3 + 2.0 * XMAX / 1330.0, 0.14, facecolor=ST.RTE, edgecolor="none", lw=0, zorder=3))
            for a, b in r._gaps: ax.add_patch(Rectangle(((a - lo) / 1e3 - 3.0, i - 0.46), (b - a) / 1e3 + 6.0, 0.92, facecolor=GAPFILL, edgecolor=GAPHATCH, hatch="/////", lw=0, zorder=4))
            if L > XMAX:                                                                   # long array cut at the axis end: break mark; code + true length written inside the pale interruption
                for dx in (-0.045, 0.0): ax.plot([XMAX * (1 + dx) - 6, XMAX * (1 + dx) + 6], [i + 0.55, i - 0.55], color="white", lw=1.2, zorder=6, solid_capstyle="butt")
                ax.text(345.0, i, f"{CODE[r.kind_R080]}, {L / 1e3:.1f} Mb", fontsize=6.5, fontweight="bold", ha="center", va="center", color=ST.NAVY, zorder=7); continue
            ax.text(L + 14 * XMAX / 1330.0, i, CODE[r.kind_R080], fontsize=7, fontweight="bold", ha="left", va="center", color=ST.NAVY, zorder=5)
            if sp == "Mac" and not r.gene and r.manuscript_name != MAC_DOM:                  # families other than the dominant one are named on their own row (no key needed)
                ax.text(L + (14 + 72) * XMAX / 1330.0, i, r.manuscript_name, fontsize=6.5, ha="left", va="center", color=ST.family_colour(r.family), zorder=5)
            if sp == "Mac" and r.gene and r.family in ("MacSat06", "GarSat49"):             # the 5S rDNA array is named (gene array, not a satellite family)
                ax.text(L + (14 + 60) * XMAX / 1330.0, i, r.family, fontsize=6.5, ha="left", va="center", color=GENE, zorder=5)   # G = gene array; "5S rDNA array" in the key / caption
            if sp == "Mac" and not r.gene and r.manuscript_name == MAC_DOM and r.lg == MAC_KEYROW:   # the dominant family is named once, on its shortest array
                cwk = 2.2 / PW * XMAX; ST.block(ax, L + 14 * XMAX / 1330.0 + 72, i - 0.36, cwk, 0.72, ST.family_colour(r.family), zorder=5)
                ax.text(L + 14 * XMAX / 1330.0 + 72 + cwk + 0.9 / PW * XMAX, i, r.manuscript_name, fontsize=6.5, ha="left", va="center", color=ST.NAVY, zorder=5)
    # overlay in mm for the keys
    ov = ST.ax_mm(fig, left_mm, top_mm, w_mm, h_mm); ov.set_xlim(left_mm, left_mm + w_mm); ov.set_ylim(top_mm + h_mm, top_mm); ov.axis("off"); ov.patch.set_alpha(0)
    lh = PITCH; cw, chh = 2.2, 1.4
    def chip(x, y, col, label, kind="rect"):
        if kind == "rect": ST.block(ov, x, y - chh / 2, cw, chh, col)
        elif kind == "gene": ov.add_patch(Rectangle((x, y - chh / 2), cw, chh, facecolor=GENE_PALE, edgecolor=GENE, hatch="////", lw=0))
        elif kind == "gap": ov.add_patch(Rectangle((x, y - chh / 2), cw, chh, facecolor=GAPFILL, edgecolor=GAPHATCH, hatch="/////", lw=0))
        elif kind == "band": ov.add_patch(Rectangle((x, y - 0.35), cw, 0.7, facecolor=BAND, edgecolor="none", lw=0))
        elif kind == "tick": ov.add_patch(Rectangle((x, y - chh / 2), cw, 0.5, facecolor=col, edgecolor="none", lw=0))
        ov.text(x + cw + 0.9, y, label, fontsize=7, ha="left", va="center", color=ST.NAVY)
    # family keys inside the free lower right of the Gar plot (families present in the v4.7 rows, ordered by bp drawn)
    fams = {sp: [] for sp in ("Mac", "Gar")}
    for r in T.sort_values("sat_kb", ascending=False).itertuples():
        if not r.gene and r.runs and r.manuscript_name not in [f for f, _ in fams[r.species]]: fams[r.species].append((r.manuscript_name, ST.family_colour(r.family)))
    for sp in ("Gar",):
        fd = dict(fams[sp]); order = [f for f, _ in fams[sp]]; n = len(order)
        xk = col0[sp] + LAB + PW * (700 / XM[sp]) + 2.0                      # right of the short arrays of LG23-28 (rows 23 / 27 / 28 carry the 'no v4.7 rows' note)
        nrows = 27
        for k, f in enumerate(order):
            yk = top_mm + (y0 - top_mm) + (nrows - n + k + 0.5) * PITCH
            chip(left_mm + xk, yk, fd[f], f)
    # shared key (chips on one line, kind codes on two lines), positions from the measured text widths
    rend = fig.canvas.get_renderer(); mm_px = fig.dpi / 25.4
    def tw(t): return t.get_window_extent(rend).width / mm_px
    ylow = y0 + 27 * PITCH + 9.2
    x = left_mm + 0.4
    for kd, c, lab in (("band", BAND, "interruption"), ("tick", ST.RTE, "RTE LINE"), ("gap", None, "assembly gap"), ("gene", None, "gene array")):
        chip(x, ylow, c, lab, kd); t = ov.texts[-1]; x += cw + 0.9 + tw(t) + 3.4
    codes = [["P composite periodic", "I interrupted irregular", "C colinear homogeneous"], ["V colinear multi-variant", "S scattered", "R single run", "G gene array"]]
    for li, line in enumerate(codes):
        xx = left_mm + 0.4; yy = ylow + (li + 1) * 2.65
        for e in line:
            a = ov.text(xx, yy, e[0], fontsize=7, fontweight="bold", ha="left", va="center", color=ST.NAVY); xx += tw(a) + 0.5
            b = ov.text(xx, yy, e[2:], fontsize=7, ha="left", va="center", color=ST.NAVY); xx += tw(b) + 3.4
    return dict(table=T.drop(columns=[c for c in T.columns if c.startswith("_")]))


if __name__ == "__main__":
    ST.house_style(); plt.rcParams["hatch.linewidth"] = 0.35
    W, H = 92.0, 84.0
    fig = ST.new_figure(W, H); res = draw(fig, 0, 0, W, H, "d")
    ST.save(fig, str(HERE / "_standalone" / "fig1_panel_d")); res["table"].to_csv(HERE / "fig1_panel_d_numbers.tsv", sep="\t", index=False)
