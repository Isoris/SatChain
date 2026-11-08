"""Figure 2, panel h (house style): homogenisation of the RTE composite unit within each centromere (2026-10-27). Size 39 x 46 mm.
H = 1 - within / between (units of one centromere vs units of the other centromeres; 15 units sampled per centromere, 200 rarefactions), with the label-shuffle null band.
Reads panel_i_table.tsv / panel_i_sweep.tsv written by panel_i.py (no work at import, no recomputation):
  F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/figure/diagnostics/docs/mockups_figure1/panels_rte/
Contract: draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None) -> dict.  Evidence L1; parameters and sweeps in fig2_panel_h_caption.md.
Standalone:  python fig2_panel_h.py -> panels/_standalone/fig2_panel_h.png/.pdf (500 dpi)"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "code"))
import satellite_style as ST

SRC = pathlib.Path(r"F:\data\satellite_FISH\LTR_MacSat03_2026-10-22\figure\diagnostics\docs\mockups_figure1\panels_rte")


def draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None):
    import pandas as pd, numpy as np
    T = pd.read_csv(SRC / "panel_i_table.tsv", sep="\t"); SW = pd.read_csv(SRC / "panel_i_sweep.tsv", sep="\t")
    T = T.sort_values(["type", "H"], ascending=[True, False]).reset_index(drop=True)
    if letter: ST.letter(fig, left_mm, top_mm, letter)
    W, Hh = fig.get_size_inches() / ST.MM
    fig.text((left_mm + 4.5) / W, 1 - (top_mm + 0.4) / Hh, "Unit homogenisation", fontsize=7, fontweight="bold", ha="left", va="top", color=ST.NAVY)
    lab_w = 8.0; ax = ST.ax_mm(fig, left_mm + lab_w, top_mm + 4.6, w_mm - lab_w - 1.5, h_mm - 4.6 - 7.6)
    n = len(T); y = np.arange(n) + (T.type.values == "II") * 0.8                      # small gap between type I and type II
    for yy, r in zip(y, T.itertuples()):
        ax.add_patch(__import__("matplotlib").patches.Rectangle((r.null_lo, yy - 0.42), r.null_hi - r.null_lo, 0.84, color=ST.GRID, lw=0, zorder=0))
        c = ST.CEN if r.type == "I" else "#3d4f66"
        ax.plot([r.H_lo, r.H_hi], [yy, yy], color=c, lw=0.5, zorder=2); ax.plot([r.H], [yy], marker="o", ms=2.2, mfc=c, mec="none", zorder=3)
        ax.text(-0.05, yy, r.chrom, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7)
    ax.set_ylim(y.max() + 0.7, -0.7); ax.set_xlim(-0.2, 1.0); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xticks([0, 0.5, 1.0]); ax.set_xticklabels(["0", "0.5", "1"]); ax.tick_params(length=2, pad=1.5); ax.axvline(0, color=ST.MUTED, lw=0.5, ls=(0, (3, 2)))
    ax.set_xlabel("H (1 = homogenised)", fontsize=7, labelpad=1.5)
    b = SW.iloc[0]
    return {"centromeres": int(b.centromeres), "units": int(b.units), "units_sampled_per_centromere": int(b.units_per_centromere_min), "H_mean": float(b.H_mean), "null_mean": float(b.null_mean),
            "p_mean_H": float(b.p_mean_H), "H_per_centromere": {r.chrom: float(r.H) for r in T.itertuples()}, "settings_swept": int(len(SW))}


if __name__ == "__main__":
    ST.house_style(); out = pathlib.Path(__file__).resolve().parent / "_standalone"; out.mkdir(exist_ok=True)
    fig = ST.new_figure(39, 46); nums = draw(fig, 0, 0, 39, 46, "h"); fig.savefig(out / "fig2_panel_h.png", dpi=500); fig.savefig(out / "fig2_panel_h.pdf"); print(nums)
