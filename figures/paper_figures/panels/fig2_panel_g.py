"""Figure 2, panel g (house style): divergence of the RTE copies, a RELATIVE age proxy (2026-10-27). Size 38 x 46 mm.
Reads the table written by panel_h.py (no work at import, no recomputation):
  F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/figure/diagnostics/docs/mockups_figure1/panels_rte/panel_h_table.tsv  (+ panel_h_sweep.tsv for the numbers dict)
Contract: draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None) -> dict of the numbers shown.  Evidence L1, parameters and sweeps are in fig2_panel_g_caption.md.
Standalone:  python fig2_panel_g.py  -> panels/_standalone/fig2_panel_g.png/.pdf (500 dpi)"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "code"))
import satellite_style as ST

SRC = pathlib.Path(r"F:\data\satellite_FISH\LTR_MacSat03_2026-10-22\figure\diagnostics\docs\mockups_figure1\panels_rte")
# table groups keep their 2026-10 names ("lineage A / B" = the GarSat03 / MacSat03 group / the GarSat04-and-5S-rDNA group of the v4.5.3 centromere table); the figure labels do not
ROWS = [("Mac, outside centromeres", "Mac\noutside", None), ("Mac, in centromeres (pooled)", "Mac cen.\npooled", "cen"), ("Mac, in centromeres (within each)", "Mac cen.\neach alone", "cen"),
        ("Gar, all copies", "Gar\nall copies", None), ("Gar, in lineage A centromeres", "Gar cen.\nGarSat03 grp", "cen"), ("Gar, in lineage B centromeres", "Gar cen.\nGarSat04", "5S")]


def draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None):
    import pandas as pd, numpy as np
    T = pd.read_csv(SRC / "panel_h_table.tsv", sep="\t").set_index("group"); SW = pd.read_csv(SRC / "panel_h_sweep.tsv", sep="\t").drop_duplicates("setting")
    if letter: ST.letter(fig, left_mm, top_mm, letter)
    lab_w = 16.5; x0 = left_mm + lab_w; ax = ST.ax_mm(fig, x0, top_mm + 5.0, w_mm - lab_w - 1.0, h_mm - 5.0 - 7.5)
    fig.text((left_mm + 4.5) / (fig.get_size_inches()[0] / ST.MM), 1 - (top_mm + 0.4) / (fig.get_size_inches()[1] / ST.MM), "RTE copy divergence", fontsize=7, fontweight="bold", ha="left", va="top", color=ST.NAVY)
    ypos = np.array([0, 1, 2, 3.3, 4.3, 5.3]); ax.set_ylim(5.8, -0.6)
    mg = float(T.mac_vs_gar_median.iloc[0]); xmax = 0.2
    for y, (g, lab, lin) in zip(ypos, ROWS):
        r = T.loc[g]; ax.barh(y, r.diversity, height=0.62, color=ST.RTE, alpha=0.55 if "each" in lab else 0.95, lw=0)
        ax.plot([r.boot_lo, r.boot_hi], [y, y], color=ST.NAVY, lw=0.5, solid_capstyle="butt")
        ax.text(-0.04, y, lab, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7, linespacing=0.95)
        if lin: ax.plot([-0.015], [y], marker="s", ms=2.6, mfc=ST.GROUP[lin], mec="none", transform=ax.get_yaxis_transform(), clip_on=False)
    ax.axvline(mg, color=ST.MUTED, lw=0.5, ls=(0, (3, 2)))
    ax.set_xlim(0, xmax); ax.set_yticks([]); ax.spines["left"].set_visible(False); ax.set_xticks([0, 0.1, 0.2]); ax.set_xticklabels(["0", "0.1", "0.2"])
    ax.set_xlabel("pairwise distance", fontsize=7, labelpad=1.5)
    ax.tick_params(length=2, pad=1.5)
    return {"groups": {g: dict(copies=int(T.loc[g].copies), median=float(T.loc[g].diversity), lo=float(T.loc[g].boot_lo), hi=float(T.loc[g].boot_hi)) for g, _, _ in ROWS},
            "mac_vs_gar_median": mg, "settings_swept": int(len(SW)), "order_unchanged_3_groups": int(SW.same_order_as_base_3.sum()), "order_unchanged_all_groups": int(SW.same_order_as_base_all.sum())}


if __name__ == "__main__":
    ST.house_style(); out = pathlib.Path(__file__).resolve().parent / "_standalone"; out.mkdir(exist_ok=True)
    fig = ST.new_figure(38, 46); nums = draw(fig, 0, 0, 38, 46, "g"); fig.savefig(out / "fig2_panel_g.png", dpi=500); fig.savefig(out / "fig2_panel_g.pdf"); print(nums)
