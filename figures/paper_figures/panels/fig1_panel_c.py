"""Figure 1c: reference centromeres, composition of every centromere of the assembly in both haplotypes (88 x 84 mm).
Redraw (house style, directly from tables) of reference_figures/centromere_composition.png.  Data (all read here, nothing pasted):
  F:/data/satellite_FISH/centromere_composition_2026-10-24/{regions,satellite_runs,te_insertions,gaps,cores}.tsv   (built by reference_figures/35b_build.py from the v4.6 catalogue)
  panels/_data/fig1c_strand_blocks.tsv                                                              (code/fig1_panel_c_prepare.py: copy-level strand blocks, same rule as 35c_figure.py)
  F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/split/rte_genome_hits_both.tsv                      (RTE LINE genome hits, gold)
  SatChain v4.7 (satchain_v4.7_2026-11-01): cat59/59_catalogue_{Mac,Gar}.tsv (name_v46 -> family) and runs_of_record/satellite_family_nomenclature.tsv (family, element_class, seq_cluster):
  family ids for the colours, gene-array status (MacSat06 and GarSat49 = 5S rDNA arrays, drawn as gene arrays) and the sequence-cluster cue above each column
  (C2.1 = the GarSat03 / MacSat03 group, C2.2 = GarSat04; the v4.6 'lineage' column of the composition tables is no longer used).
y axis = Mb from the first run of the centromeric cluster (0), flanks +-150 kb of the original region are kept (negative values = the 150 kb before the first run).
Per chromosome, left to right: strand strip (copy-level strand blocks, majority / minority) | satellite runs by family on a pale-grey centromere | transposable element strip | homogenised-core bars.
usage: python fig1_panel_c.py   (renders the panel alone: panels/_standalone/fig1_panel_c.png / .pdf at 500 dpi + fig1_panel_c_numbers.tsv)"""
import sys, os, pathlib
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "code")); import satellite_style as ST
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

CC = pathlib.Path("F:/data/satellite_FISH/centromere_composition_2026-10-24")
RTE_HITS = "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/split/rte_genome_hits_both.tsv"
MB = 1e6
GAPCOL = "#7a7467"; BACK = "#e7e4db"; GENE = "#e6a100"
TECOL = {"LTR": "#6F93B8", "DNA": "#A58BB8", "Helitron": "#78A68B", "other": ST.GREY}
STRAND_MAJ = "#7C8A93"; STRAND_MIN = ST.TERRACOTTA


def te_class(t):
    if "LTR" in t or "DIRS" in t or "Penelope" in t or "pararetro" in t: return "LTR"
    if "TIR" in t or "Sola" in t or "polinton" in t: return "DNA"
    if "helitron" in t.lower(): return "Helitron"
    return "other"                              # repeat fragments and LINE calls that are not the RTE LINE


def load():
    R = pd.read_csv(CC / "regions.tsv", sep="\t"); S = pd.read_csv(CC / "satellite_runs.tsv", sep="\t"); T = pd.read_csv(CC / "te_insertions.tsv", sep="\t")
    G = pd.read_csv(CC / "gaps.tsv", sep="\t"); K = pd.read_csv(CC / "cores.tsv", sep="\t"); SB = pd.read_csv(HERE / "_data" / "fig1c_strand_blocks.tsv", sep="\t")
    S["fid"] = [ST.pipeline_name(f) for f in S.family]                                   # v4.6 manuscript code -> pipeline family (v4.7 catalogue name_v46 + nomenclature)
    S["isgene"] = S.family.str.startswith("gene_array") | S.fid.map(ST.is_gene_array)   # v4.7: MacSat06 / GarSat49 are 5S rDNA gene arrays
    S["colour"] = [ST.family_colour(f) if not str(f).startswith("gene_array") else ST.GREY for f in S.fid]
    S["mname"] = S.fid                                                                   # one name everywhere: pipeline name
    S["cluster"] = S.fid.map(ST.family_cluster)                                         # v4.7 seq_cluster (C2.1 / C2.2 / none / na)
    H = pd.read_csv(RTE_HITS, sep="\t"); H["a"] = H[["sstart", "send"]].min(axis=1); H["b"] = H[["sstart", "send"]].max(axis=1)
    return R, S, T, G, K, SB, H


def merge(iv, gap=0):
    iv = sorted(iv); out = []
    for a, b in iv:
        if out and a <= out[-1][1] + gap: out[-1][1] = max(out[-1][1], b)
        else: out.append([a, b])
    return out


def numbers_per_chrom(R, S, T, G, K, SB, H):
    rows = []
    K5 = K[(K.W == 20000) & (K.k == 21) & (K.tau == 0.5)]; K3 = K[(K.W == 20000) & (K.k == 21) & (K.tau == 0.3)]
    for r in R.itertuples():
        s = S[(S.sp == r.sp) & (S.chrom == r.chrom)]; sb = SB[(SB.sp == r.sp) & (SB.chrom == r.chrom)]; core = K5[(K5.sp == r.sp) & (K5.chrom == r.chrom)]
        fam = s.assign(l=s.e - s.s).groupby("mname").l.sum().sort_values(ascending=False)
        h = H[(H.chrom == r.chrom) & (H.b > r.lo) & (H.a < r.hi)]
        t = T[(T.sp == r.sp) & (T.chrom == r.chrom)]
        clu = s[s.mname == fam.index[0]].cluster.iloc[0] if len(fam) else ""
        rows.append(dict(species=r.sp, chrom=r.chrom, region_start=r.lo, region_end=r.hi, first_run=r.cen_lo, last_run=r.cen_hi, region_kb=round((r.hi - r.lo) / 1e3, 1), n_runs=len(s), satellite_kb=round((s.e - s.s).sum() / 1e3, 1),
                         main_family=fam.index[0] if len(fam) else "", main_family_kb=round(fam.iloc[0] / 1e3, 1) if len(fam) else 0, main_family_seq_cluster_v47=clu,
                         main_family_class_v47=("gene_array" if len(fam) and ST.is_gene_array(fam.index[0]) else "satellite"), n_families=len(fam), gene_array_kb=round((s[s.isgene].e - s[s.isgene].s).sum() / 1e3, 1),
                         minority_strand_kb=round((sb[~sb.major].e - sb[~sb.major].s).sum() / 1e3, 1), majority_strand_kb=round((sb[sb.major].e - sb[sb.major].s).sum() / 1e3, 1),
                         core_h05_kb=round((core.e - core.s).sum() / 1e3, 1), core_h03_kb=round((K3[(K3.sp == r.sp) & (K3.chrom == r.chrom)].e - K3[(K3.sp == r.sp) & (K3.chrom == r.chrom)].s).sum() / 1e3, 1),
                         n_assembly_gaps=len(G[(G.sp == r.sp) & (G.chrom == r.chrom)]), n_RTE_hits=len(h), n_TE_features=len(t)))
    return pd.DataFrame(rows)


def sweep_core(K):
    """stability of the homogenised-core call over window W, k and threshold tau: chromosomes with a core (>= 50 kb) and total core Mb"""
    rows = []
    for (W, k, tau, sp), d in K.groupby(["W", "k", "tau", "sp"]): rows.append(dict(W=W, k=k, tau=tau, species=sp, chromosomes_with_core=d.chrom.nunique(), core_Mb=round((d.e - d.s).sum() / 1e6, 2)))
    return pd.DataFrame(rows)


def draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None):
    R, S, T, G, K, SB, H = load()
    ST.letter(fig, left_mm, top_mm, letter) if letter else None
    num = numbers_per_chrom(R, S, T, G, K, SB, H)
    K5 = K[(K.W == 20000) & (K.k == 21) & (K.tau == 0.5)]; K3 = K[(K.W == 20000) & (K.k == 21) & (K.tau == 0.3)]
    W_, H_ = fig.get_size_inches() / ST.MM
    LM, RM = 9.0, 0.4; NCOL = 28.4; plot_w = w_mm - LM - RM; mm_per_col = plot_w / NCOL
    KEY_H = 14.5; HDR = 7.0; GAP = 1.5
    # scale: mm per Mb common to both species
    rng = {}
    for sp in ("Mac", "Gar"):
        rr = R[R.sp == sp]; rng[sp] = (-(rr.cen_lo - rr.lo).max() / MB, ((rr.hi - rr.cen_lo)).max() / MB)
    tot = sum(b - a for a, b in rng.values()); avail = h_mm - 4 - 2 * HDR - KEY_H - GAP - 1.5
    s_mm = avail / (tot + 0.0)
    extra = {"Mac": 2.0, "Gar": 3.5}                    # room for the key inside the plots (mm)
    s_mm = (avail - sum(extra.values())) / tot
    yy = top_mm + 4
    axes = {}
    out = {}
    for sp, title, ncol in (("Mac", "C. macrocephalus", 27), ("Gar", "C. gariepinus", 28)):
        a, b = rng[sp]; hh = (b - a) * s_mm + extra[sp]
        ax = ST.ax_mm(fig, left_mm + LM, yy + HDR, plot_w, hh); axes[sp] = ax
        ax.set_xlim(-0.12, NCOL - 0.12); ax.set_ylim(b + extra[sp] / s_mm, a); ax.set_facecolor("none")
        for s_ in ("top", "right", "bottom"): ax.spines[s_].set_visible(False)
        ax.spines["left"].set_bounds(0, np.floor(b * 2) / 2); ax.spines["left"].set_position(("data", -0.12))
        ax.set_xticks([]); ax.set_yticks(np.arange(0, b + 0.01, 0.5)); ax.set_yticklabels(["0" if v == 0 else f"{v:.1f}" for v in np.arange(0, b + 0.01, 0.5)], fontsize=7)
        ax.set_ylabel("Mb from first run", fontsize=7); ax.yaxis.set_label_coords(-(LM - 3.6) / plot_w, 0.4)
        rr = R[R.sp == sp].copy(); rr["n"] = rr.chrom.str.extract(r"LG(\d+)")[0].astype(int)
        slot = {n: i for i, n in enumerate(range(1, ncol + 1))}
        patches = {k: [] for k in ("back", "sat", "gene", "strandM", "strandm", "te_other", "te_LTR", "te_DNA", "te_Helitron", "te_RTE")}; hatch = []; cols = {}
        for r in rr.itertuples():
            x0 = slot[r.n] + 0.08; c0 = r.cen_lo; yv = lambda p: (p - c0) / MB
            patches["back"].append(Rectangle((x0 + 0.14, yv(r.lo)), 0.44, yv(r.hi) - yv(r.lo)))
            for q in S[(S.sp == sp) & (S.chrom == r.chrom)].itertuples():
                a_, b_ = max(q.s, r.lo), min(q.e, r.hi)
                if b_ <= a_: continue
                if q.isgene: ax.add_patch(Rectangle((x0 + 0.14, yv(a_)), 0.44, yv(b_) - yv(a_), facecolor="#f4e2a8", edgecolor=GENE, hatch="////", lw=0, zorder=2))
                else:
                    ax.add_patch(Rectangle((x0 + 0.14, yv(a_)), 0.44, yv(b_) - yv(a_), facecolor=q.colour, edgecolor="none", lw=0, zorder=2))
                    if isinstance(q.gene, str) and q.gene != "none": ax.add_patch(Rectangle((x0 + 0.14, yv(a_)), 0.44, yv(b_) - yv(a_), fill=False, hatch="////", edgecolor="white", lw=0, zorder=3))
            for q in SB[(SB.sp == sp) & (SB.chrom == r.chrom)].itertuples():
                patches["strandM" if q.major else "strandm"].append(Rectangle((x0 + (0.0 if q.major else 0.06), yv(q.s)), 0.07, max(yv(q.e) - yv(q.s), 0.004)))
            t = T[(T.sp == sp) & (T.chrom == r.chrom)]
            by = {}
            for q in t.itertuples(): by.setdefault(te_class(q.type), []).append((q.s, q.e))
            for k_, iv in by.items():
                for a_, b_ in merge(iv, 300): patches["te_" + k_].append(Rectangle((x0 + 0.62, yv(a_)), 0.16, max(yv(b_) - yv(a_), 0.004)))
            h = H[(H.chrom == r.chrom) & (H.b > r.lo) & (H.a < r.hi)]
            for a_, b_ in merge(list(zip(h.a, h.b)), 300): patches["te_RTE"].append(Rectangle((x0 + 0.62, yv(a_)), 0.16, max(yv(b_) - yv(a_), 0.004)))
            for g in G[(G.sp == sp) & (G.chrom == r.chrom)].itertuples(): ax.plot([x0 + 0.02, x0 + 0.80], [yv(g.s)] * 2, color=GAPCOL, lw=0.5, ls=(0, (1.2, 1.2)), zorder=5)
            for kk, lw, off, col in ((K3, 0.5, 0.88, "#A9C9C4"), (K5, 0.9, 0.95, ST.ROH[3])):
                for c in kk[(kk.sp == sp) & (kk.chrom == r.chrom)].itertuples(): ax.plot([x0 + off] * 2, [yv(c.s), yv(c.e)], color=col, lw=lw, solid_capstyle="butt", zorder=4)
            # cue (thin bar above the column) = v4.7 sequence cluster of the main family (C2.1 = GarSat03 / MacSat03 group, C2.2 = GarSat04), gold when gene arrays dominate
            nr = num[(num.species == sp) & (num.chrom == r.chrom)].iloc[0]
            lc = GENE if nr.gene_array_kb > 0.5 * nr.satellite_kb else {"C2.1": ST.CEN, "C2.2": ST.FIVE_S}.get(nr.main_family_seq_cluster_v47, ST.GREY)
            ax.add_patch(Rectangle((x0 + 0.14, a - 0.05), 0.44, 0.035, facecolor=lc, edgecolor="none", lw=0, clip_on=False))
        for k_, col in (("back", BACK), ("strandM", STRAND_MAJ), ("strandm", STRAND_MIN), ("te_other", TECOL["other"]), ("te_LTR", TECOL["LTR"]), ("te_DNA", TECOL["DNA"]), ("te_Helitron", TECOL["Helitron"]), ("te_RTE", ST.RTE)):
            if patches[k_]: ax.add_collection(PatchCollection(patches[k_], facecolor=col, edgecolor="none", lw=0, zorder=1 if k_ == "back" else 3))
        # chromosome numbers, staggered in two rows (odd up, even down)
        for n in range(1, ncol + 1):
            have = ((rr.n == n).any())
            ax.annotate(str(n), xy=(slot[n] + 0.08 + 0.36, a - 0.05), xytext=(0, 3.4 if n % 2 == 0 else 11.4), textcoords="offset points", ha="center", va="bottom",
                        fontsize=7, fontweight="bold", color=ST.NAVY if have else ST.MUTED, annotation_clip=False)
        out[sp] = dict(n_centromeres=len(rr), mm_per_Mb=round(s_mm, 2))
        yy += HDR + hh + GAP
    # ---- overlay axes in mm for keys, species names
    ov = ST.ax_mm(fig, left_mm, top_mm, w_mm, h_mm); ov.set_xlim(0, w_mm); ov.set_ylim(h_mm, 0); ov.axis("off"); ov.patch.set_alpha(0)
    ch_w, ch_h, lh = 2.2, 1.5, 2.65
    def chip(x, y, col, label, kind="rect", fs=7, hatchc=None):
        yc = y
        if kind == "rect": ov.add_patch(Rectangle((x, yc - ch_h / 2), ch_w, ch_h, facecolor=col, edgecolor="none", lw=0))
        elif kind == "hatch": ov.add_patch(Rectangle((x, yc - ch_h / 2), ch_w, ch_h, facecolor=col, edgecolor=hatchc, hatch="////", lw=0))
        elif kind == "line": ov.plot([x + ch_w / 2] * 2, [yc - ch_h / 2 - 0.2, yc + ch_h / 2 + 0.2], color=col, lw=hatchc or 1.0, solid_capstyle="butt")
        elif kind == "dots": ov.plot([x - 0.3, x + ch_w + 0.3], [yc, yc], color=col, lw=0.5, ls=(0, (1.2, 1.2)))
        ov.text(x + ch_w + 0.9, yc, label, fontsize=fs, va="center", ha="left", color=ST.NAVY)
    # species names + family keys inside the empty lower-right of each plot
    # family keys: the satellite families that cover the most bp in the regions (gene arrays are keyed once in the shared key: MacSat06 / GarSat49 = 5S rDNA arrays)
    fam_keys = {}
    for sp in ("Mac", "Gar"):
        s = S[(S.sp == sp) & (~S.isgene)]; fam_keys[sp] = list(s.assign(l=s.e - s.s).groupby("mname").l.sum().sort_values(ascending=False).index[:4])
    fam_col = S.drop_duplicates("mname").set_index("mname").colour
    yb = {}
    for sp, title in (("Mac", "C. macrocephalus"), ("Gar", "C. gariepinus")):
        pos = axes[sp].get_position(); y_top = (1 - pos.y1) * H_ - top_mm; y_bot = (1 - pos.y0) * H_ - top_mm; yb[sp] = (y_top, y_bot)
        xr = w_mm - RM; x0k = xr - 43.8
        ov.text(xr, y_bot - 2 * lh - 1.7, title, fontsize=7, fontstyle="italic", fontweight="bold", ha="right", va="center", color=ST.NAVY)
        for i, f in enumerate(fam_keys[sp]):
            col = i // 2; row = i % 2
            chip(x0k + col * 22.2, y_bot - (2 - row) * lh + 1.0, fam_col[f], f)
    # shared key block; the two cue bars = v4.7 sequence clusters of the main family (C2.1 = GarSat03 / MacSat03 group, C2.2 = GarSat04)
    ky = h_mm - KEY_H + lh / 2 + 0.3; cx = [0.4, 26.0, 55.0]
    items = [
        [("rect", ST.RTE, "RTE LINE", None), ("rect", TECOL["LTR"], "LTR element", None), ("rect", TECOL["DNA"], "DNA (TIR) element", None), ("rect", TECOL["Helitron"], "Helitron", None), ("rect", TECOL["other"], "other repeat", None)],
        [("rect", STRAND_MAJ, "copies, majority strand", None), ("rect", STRAND_MIN, "copies, minority strand", None), ("line", ST.ROH[3], "core, h ≥ 0.5 (≥ 50 kb)", 0.9), ("line", "#A9C9C4", "core, h ≥ 0.3", 0.5), ("dots", GAPCOL, "assembly gap", None)],
        [("hatch", "#f4e2a8", "gene array (5S rDNA, tRNA)", GENE), ("hatch", ST.CENTROMERE, "5S gene call in a run", "white"), ("rect", BACK, "non-satellite", None), ("bar", ST.CEN, "GarSat03 / MacSat03 group", None), ("bar", ST.FIVE_S, "GarSat04 (5S-derived)", None)]]
    for c, col_items in enumerate(items):
        for r_, (kind, colr, label, extra) in enumerate(col_items):
            y = ky + r_ * lh; x = cx[c]; fs = 6.5
            if kind == "bar": ov.add_patch(Rectangle((x, y - 0.35), ch_w, 0.7, facecolor=colr, edgecolor="none")); ov.text(x + ch_w + 0.9, y, label, fontsize=fs, va="center", ha="left", color=ST.NAVY)
            elif kind == "hatch" and extra == "white":
                ov.add_patch(Rectangle((x, y - ch_h / 2), ch_w, ch_h, facecolor=ST.family_colour("GarSat04"), edgecolor="none")); ov.add_patch(Rectangle((x, y - ch_h / 2), ch_w, ch_h, fill=False, hatch="////", edgecolor="white", lw=0))
                ov.text(x + ch_w + 0.9, y, label, fontsize=fs, va="center", ha="left", color=ST.NAVY)
            else: chip(x, y, colr, label, kind, fs=fs, hatchc=extra)
    out["numbers"] = num
    return out


if __name__ == "__main__":
    plt.rcParams["hatch.linewidth"] = 0.35
    ST.house_style(); plt.rcParams["hatch.linewidth"] = 0.35
    W, H = 88.0, 84.0
    fig = ST.new_figure(W, H); res = draw(fig, 0, 0, W, H, "c")
    base = str(HERE / "_standalone" / "fig1_panel_c"); ST.save(fig, base)
    res["numbers"].to_csv(HERE / "fig1_panel_c_numbers.tsv", sep="\t", index=False)
    sw = sweep_core(pd.read_csv(CC / "cores.tsv", sep="\t")); sw.to_csv(HERE / "fig1_panel_c_core_sweep.tsv", sep="\t", index=False)
    print({k: v for k, v in res.items() if k != "numbers"})
