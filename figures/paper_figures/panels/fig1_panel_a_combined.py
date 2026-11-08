#!/usr/bin/env python3
"""Figure 1 panel a, COMBINED variant (92 x 56 mm): the two trees (5S-related families; centromeric GarSat03 / MacSat03 group) as ONE figure on ONE common scale (identical mm per
substitution/site) with one scale bar and one key.  The two groups are NOT joined into one tree on purpose: no joint alignment exists (family consensus max local identity between the
two groups is 53-58 %, i.e. background, versus 78-82 % inside a group; 04_family_cross_identity.tsv), so a single tree would have an arbitrary, unsupported join.  Same data / layout cache as fig1_panel_a.py.
usage: draw(fig, left_mm, top_mm, w_mm, h_mm, letter) -> dict"""
import sys, pathlib
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np
import fig1_panel_a as A
ST = A.ST


def draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None):
    ax = ST.ax_mm(fig, left_mm, top_mm, w_mm, h_mm); ax.set_xlim(0, w_mm); ax.set_ylim(h_mm, 0); ax.axis("off")
    if letter: ST.letter(fig, left_mm, top_mm, letter)
    trees = {k: A.load_tree(k) for k in ("B", "A")}; bh = {k: A.CELL_Y1 - A.CELL_Y0 - 2 * A.TREE_MY[k] for k in 'BA'}
    fit = []
    for k, d in trees.items():
        lo, hi = d["segs"].reshape(-1, 2).min(0), d["segs"].reshape(-1, 2).max(0)
        fit.append(min((A.COLS[k][1] - 2 * A.TREE_MX[k]) / (hi[0] - lo[0]), bh[k] / (hi[1] - lo[1])) * 0.98)
    sc = min(fit); rows = []; labs = {}; models = {}
    for k, key in enumerate(("B", "A")):
        x0, colw = A.COLS[key]
        ax.text(x0 + 0.3, 5.6, A.LIN[key][1], fontsize=7, fontweight="bold", color=A.GROUPCOL[key], ha="left", va="center")
        r = A.draw_lineage(ax, key, (x0, A.CELL_Y0, colw, A.CELL_Y1 - A.CELL_Y0), (x0 + A.TREE_MX[key], A.CELL_Y0 + A.TREE_MY[key], colw - 2 * A.TREE_MX[key], bh[key]), scale=sc)
        labs[key] = r["labels"]; models[key] = r["model"]; A.family_key(ax, x0 + 0.3, A.FAM_Y, r["fams"], w=colw - 0.6)
        rows += A.numbers_rows("a_combined", key, trees[key])
    L = A.scalebar(ax, 0.8, A.SCALE_Y, sc, label=False); ax.text(0.8 + L + 0.8, A.SCALE_Y, "0.05 substitutions/site (both trees, same scale)", fontsize=6, ha="left", va="center")
    xe = A.key_row(ax, 0.8, A.K1_Y); A.support_key(ax, xe + 1.0, A.K1_Y); A.ring_key(ax, 0.8, A.K2_Y)
    return dict(rows=rows, mm_per_unit=sc, label_positions_mm=labs, models=models)


if __name__ == "__main__":
    ST.house_style(); W, H = 92, 56
    fig = ST.new_figure(W, H); num = draw(fig, 0, 0, W, H, "a")
    ST.save(fig, str(HERE / "_standalone" / "fig1_panel_a_combined")); A.write_tables(num, "fig1_panel_a_combined"); print(num["mm_per_unit"], num["label_positions_mm"])
