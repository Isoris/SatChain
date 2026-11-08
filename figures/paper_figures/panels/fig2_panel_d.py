"""Figure 2, panel d (house style module): circular tree of 716 RTE core copies (316 C. macrocephalus + 400 C. gariepinus), 66 x 66 mm.
   draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None) -> dict of the numbers shown.   Standalone:  python fig2_panel_d.py  -> panels/_standalone/fig2_panel_d.{png,pdf} + fig2_panel_d_numbers.tsv
Data code = panels_rte/panel_d.py (imported, not edited): tree reading, midpoint re-rooting (display only), tree classes, pairwise distances, parsimony (Sankoff), clade statistics.
The label-shuffle null (1,000 permutations) is NOT re-run here: its results are read from panels_rte/panel_d_stats.tsv; the observed values are recomputed from the tree and checked against that table.
Tree = LTR_MacSat03_2026-10-22/tree/copies.treefile (IQ-TREE -fast, SH-aLRT 1000; read, not re-run)."""
import sys, pathlib, importlib.util, functools
import numpy as np, pandas as pd
from matplotlib.collections import LineCollection, PolyCollection
import matplotlib.patches as mp
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "code"))
import satellite_style as ST

OLDDIR = pathlib.Path(r"F:\data\satellite_FISH\LTR_MacSat03_2026-10-22\figure\diagnostics\docs\mockups_figure1\panels_rte")
OLIVE = "#8a9a5b"; OUTSIDE = "#D9DCDF"
POS = {"A": ST.CEN, "B": ST.FIVE_S, "other": OLIVE, "arm": OUTSIDE}      # tip categories of panel_d.py (v4.5.3 intervals): A = GarSat03 / MacSat03-group centromere, B = GarSat04 / 5S rDNA-array centromere
SUP, RMAX = 80, 0.225


def _old():
    spec = importlib.util.spec_from_file_location("_old_panel_d", OLDDIR / "panel_d.py"); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


@functools.lru_cache(maxsize=1)
def prepare():
    o = _old()
    names, adj = o.read_unrooted(o.TREE); root0_deg = len(adj[0]); names, adj, m, D0 = o.midpoint_root(names, adj)
    par, kids, blen, sup, order = o.rooted(names, adj, m); t = o.T(names, par, kids, blen, sup, order, m)
    nm = [names[i] for i in t.tips]; sp = np.array([x.split("|")[0] for x in nm]); pos = np.array([x.split("|")[3] for x in nm]); fam = np.array([x.split("|")[4] for x in nm])
    cat = np.array([("arm" if p == "arm" else o.LIN.get(f, "other")) for p, f in zip(pos, fam)])
    # observed statistics recomputed (the null comes from the table)
    D = o.pdist(t); Ds = D.copy(); np.fill_diagonal(Ds, np.inf); nn = Ds.argmin(1)
    lab = (sp == "Mac").astype(int)[:, None]
    obs = dict(NN_other_species_share=float(np.mean(sp[nn] != sp)), parsimony_species_changes=int(o.sankoff(t, lab, np.ones(len(sp), bool))[0]))
    # ladderize, angles
    size = t.M.sum(1)
    for u in range(t.n): t.kids[u].sort(key=lambda c: size[c])
    leaves = []; sys.setrecursionlimit(10000)
    def walk(u):
        if t.names[u] is not None: leaves.append(u); return
        for c in t.kids[u]: walk(c)
    walk(t.root); n = len(leaves); ang = {}
    for i, u in enumerate(leaves): ang[u] = 2 * np.pi * (i + 0.5) / n
    for u in t.post:
        if t.names[u] is None: ks = [ang[c] for c in t.kids[u]]; ang[u] = 0.5 * (min(ks) + max(ks))
    pure = {}
    for u in range(t.n):
        ss = set(sp[t.M[u]]); pure[u] = ss.pop() if len(ss) == 1 else None
    stats = pd.read_csv(OLDDIR / "panel_d_stats.tsv", sep="\t")
    return dict(o=o, t=t, sp=sp, cat=cat, leaves=leaves, ang=ang, pure=pure, n=n, obs=obs, stats=stats, D0=D0, root0_deg=root0_deg)


def stat(S, scope, metric, col="observed"):
    r = S[(S.scope == scope) & (S.metric == metric)]; return float(r[col].iloc[0])


def draw(fig, left_mm, top_mm, w_mm=66.0, h_mm=66.0, letter=None):
    P = prepare(); t, sp, cat, leaves, ang, pure, n, S = P["t"], P["sp"], P["cat"], P["leaves"], P["ang"], P["pure"], P["n"], P["stats"]
    if letter: ST.letter(fig, left_mm, top_mm, letter)
    ov = ST.ax_mm(fig, left_mm, top_mm, w_mm, h_mm); ov.set_xlim(0, w_mm); ov.set_ylim(h_mm, 0); ov.axis("off"); ov.patch.set_alpha(0)
    ov.text(4.2, 0.2, "Tree of RTE core copies, both species", fontsize=7, fontweight="bold", ha="left", va="top")
    # ---- circle geometry: diameter DIA mm incl. rings, centred horizontally, below the title
    DIA, CY = 49.0, 4.6 + 49.0 / 2
    ax = ST.ax_mm(fig, left_mm + (w_mm - DIA) / 2, top_mm + 4.6, DIA, DIA); ax.set_aspect("equal"); ax.axis("off"); ax.patch.set_alpha(0)
    R = RMAX; r_sp0, r_sp1 = R * 1.03, R * 1.085; r_ps0, r_ps1 = R * 1.095, R * 1.155
    lim = r_ps1 * 1.005; ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    xy = lambda r, a: (r * np.cos(a), r * np.sin(a)); dep = lambda u: min(t.depth[u], RMAX)
    segs, cols = [], []
    for u in t.order:
        if t.names[u] is None:
            ks = t.kids[u]; a0, a1 = min(ang[c] for c in ks), max(ang[c] for c in ks); aa = np.linspace(a0, a1, max(3, int((a1 - a0) / 0.02) + 2))
            if len(ks) > 1: segs.append(np.c_[t.depth[u] * np.cos(aa), t.depth[u] * np.sin(aa)]); cols.append(ST.SPECIES[pure[u]] if pure[u] else ST.MUTED)
        if t.par[u] >= 0:
            p = t.par[u]; x0, y0 = xy(t.depth[p], ang[u]); x1, y1 = xy(dep(u), ang[u])
            segs.append(np.array([[x0, y0], [x1, y1]])); cols.append(ST.SPECIES[pure[u]] if pure[u] else ST.MUTED)
    guide = [np.array([xy(dep(u), ang[u]), xy(r_sp0, ang[u])]) for u in leaves if r_sp0 - dep(u) > 0.01 * R]
    ax.add_collection(LineCollection(guide, colors="#ECEDEE", linewidths=0.12, zorder=0, rasterized=True))
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=0.22, capstyle="butt", zorder=1, rasterized=True))
    da = 2 * np.pi / n
    def ring(r0, r1, colf):
        polys, cc = [], []
        for k, u in enumerate(leaves):
            a0, a1 = ang[u] - da / 2 * 1.04, ang[u] + da / 2 * 1.04
            polys.append([xy(r0, a0), xy(r1, a0), xy(r1, a1), xy(r0, a1)]); cc.append(colf(u))
        ax.add_collection(PolyCollection(polys, facecolors=cc, edgecolors="none", linewidths=0, zorder=2, rasterized=True))
    ring(r_sp0, r_sp1, lambda u: ST.SPECIES[sp[t.tix[u]]]); ring(r_ps0, r_ps1, lambda u: POS[cat[t.tix[u]]])
    clip = [u for u in leaves if t.depth[u] > RMAX]
    for u in clip:
        a = ang[u]; r = R * 0.93; hw = 0.004
        ax.plot([xy(r, a - hw)[0], xy(r, a + hw)[0]], [xy(r, a - hw)[1], xy(r, a + hw)[1]], color="white", lw=0.7, solid_capstyle="butt", zorder=4, rasterized=True)
    dots = np.array([xy(t.depth[u], ang[u]) for u in t.order if t.names[u] is None and u != t.root and t.sup[u] >= SUP])
    ax.scatter(dots[:, 0], dots[:, 1], s=0.45, c=ST.NAVY, linewidths=0, zorder=3, rasterized=True)
    # ---- key block (bottom, 3 rows) ; mm per data unit gives the true length of the scale bar
    mm_per = DIA / (2 * lim); sb = 0.05 * mm_per
    ky0, dy = 56.4, 2.85; C2 = 33.0
    def sw(x, y, col, lab, italic=False, fs=7):
        ov.add_patch(mp.Rectangle((x, y - 1.05), 2.1, 2.1, fc=col, ec="none")); ov.text(x + 3.0, y, lab, fontsize=fs, va="center", ha="left", style="italic" if italic else "normal")
    sw(0, ky0, ST.SPECIES["Mac"], "C. macrocephalus", True); sw(C2, ky0, ST.SPECIES["Gar"], "C. gariepinus", True)
    sw(0, ky0 + dy, POS["A"], "in GarSat03/MacSat03 cen.", fs=6.5); sw(C2, ky0 + dy, POS["B"], "in GarSat04/5S rDNA cen.", fs=6.5)
    sw(0, ky0 + 2 * dy, POS["other"], "in other-family centromere"); sw(C2, ky0 + 2 * dy, POS["arm"], "outside centromere")
    y4 = ky0 + 3 * dy
    ov.plot([1.05], [y4], marker="o", ms=1.2, mfc=ST.NAVY, mec="none", ls=""); ov.text(3.0, y4, "branch SH-aLRT ≥ 80", fontsize=7, va="center", ha="left")
    ov.plot([C2, C2 + sb], [y4] * 2, color=ST.NAVY, lw=0.6, solid_capstyle="butt"); ov.text(C2 + sb + 1.0, y4, "0.05 subst./site", fontsize=7, va="center", ha="left")
    # ---- numbers shown
    allS = "all copies"; arm = "arm copies only (outside centromere +-1 Mb)"
    out = dict(copies=len(sp), copies_Mac=int((sp == "Mac").sum()), copies_Gar=int((sp == "Gar").sum()),
               Mac_in_cen_group_centromere=int(((sp == "Mac") & (cat == "A")).sum()), Gar_in_cen_group_centromere=int(((sp == "Gar") & (cat == "A")).sum()),
               Mac_in_5S_group_centromere=int(((sp == "Mac") & (cat == "B")).sum()), Gar_in_5S_group_centromere=int(((sp == "Gar") & (cat == "B")).sum()),
               Mac_other_family_centromere=int(((sp == "Mac") & (cat == "other")).sum()), Gar_other_family_centromere=int(((sp == "Gar") & (cat == "other")).sum()),
               Mac_outside=int(((sp == "Mac") & (cat == "arm")).sum()), Gar_outside=int(((sp == "Gar") & (cat == "arm")).sum()),
               supported_nodes_SH_ge_80_drawn=len(dots), branches_cut_at_0p225=len(clip), max_depth_from_midpoint=round(float(max(t.depth[u] for u in t.tips)), 4),
               tree_diameter_midpoint_path=round(float(P["D0"]), 4), scale_bar_substitutions_per_site=0.05,
               rooting="midpoint of the longest tip-to-tip path, display only (file root = trifurcation, no outgroup)")
    # verdict numbers (caption + numbers table, not in the picture): observed recomputed, null from panel_d_stats.tsv
    out.update(NN_other_species_observed_recomputed=round(P["obs"]["NN_other_species_share"], 5), parsimony_changes_observed_recomputed=P["obs"]["parsimony_species_changes"],
               NN_other_species_all=stat(S, allS, "NN_other_species_share"), NN_other_species_all_null=stat(S, allS, "NN_other_species_share", "null_mean"),
               parsimony_all=stat(S, allS, "parsimony_species_changes"), parsimony_all_null=stat(S, allS, "parsimony_species_changes", "null_mean"),
               NN_other_species_arm=stat(S, arm, "NN_other_species_share"), NN_other_species_arm_null=stat(S, arm, "NN_other_species_share", "null_mean"),
               parsimony_arm=stat(S, arm, "parsimony_species_changes"), parsimony_arm_null=stat(S, arm, "parsimony_species_changes", "null_mean"))
    out["verdict_matches_table"] = bool(abs(out["NN_other_species_observed_recomputed"] - out["NN_other_species_all"]) < 1e-4 and out["parsimony_changes_observed_recomputed"] == out["parsimony_all"])
    return out


if __name__ == "__main__":
    ST.house_style(); W = H = 66.0
    fig = ST.new_figure(W, H); res = draw(fig, 0, 0, W, H, "d")
    ST.save(fig, str(HERE / "_standalone" / "fig2_panel_d"))
    P = prepare(); S = P["stats"]
    pd.DataFrame({"metric": list(res), "value": list(res.values())}).to_csv(HERE / "fig2_panel_d_numbers.tsv", sep="\t", index=False)
    S.to_csv(HERE / "fig2_panel_d_stats_full.tsv", sep="\t", index=False)
    print({k: v for k, v in res.items() if "recomputed" in k or "verdict" in k or "NN" in k or "pars" in k})
