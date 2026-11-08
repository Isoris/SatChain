"""Figure 1 / Figure 2 in the numbered-story layout of the hatchery population figure (2026-10-30): each row carries a numbered banner with its one-line conclusion,
the panels are re-lettered in reading order, panel modules are unchanged (they draw into the box they are given).
usage: python assemble_v2.py [1|2|both]  -> figures/Figure{1,2}_v2.png (500 dpi) + pdf + svg.   The old assemble.py and Figure*_draft.* are untouched."""
import sys, pathlib, importlib, traceback
HERE = pathlib.Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "panels"))
import satellite_style as ST
from matplotlib.lines import Line2D

BAN = 7.0   # banner height (mm)

# (new letter, module, left, top, w, h)
# banners (SatChain v4.7 wording, 2026-11-01: no "lineage"; the two groups are GarSat04 (5S-derived satellite, sequence cluster C2.2) and the GarSat03 / MacSat03 group (C2.1))
FIG = {
 1: dict(w=183, banners=[(0, "1", "Two unrelated satellite groups, GarSat04 and the GarSat03 / MacSat03 group, each expanded separately in the two species"),
                         (68, "2", "What the units carry: GarSat04 a 5S rRNA gene segment, MacSat06 a whole 5S gene (a 5S rDNA array), MacSat03 an RTE LINE"),
                         (132, "3", "GarSat03 and MacSat03 form the candidate centromeric arrays, organised differently in the two species")],
         panels=[("a", "fig1_panel_a", 0, 8, 92, 56), ("b", "fig1_panel_b", 95, 8, 88, 56),
                 ("c", "fig1_panel_g", 0, 76, 55, 52), ("d", "fig1_panel_f", 59, 76, 56, 52), ("e", "fig1_panel_e", 118, 76, 65, 52),   # 2026-11-01: the unit panel joins row 2 (was alone at the bottom)
                 ("f", "fig1_panel_c", 0, 140, 88, 84), ("g", "fig1_panel_d", 91, 140, 92, 84)], h=228),
}
# Figure 2: story = 1 the element and the unit it builds (a, b, e); 2 where it sits and how it spread (c, d); 3 one array, and age and homogenisation (f, g, h)
FIG[2] = dict(w=183, banners=[(0, "1", "An RTE-clade LINE is built into the MacSat03 centromeric unit"),
                              (62, "2", "The element is concentrated in GarSat03 / MacSat03 centromeres and amplified separately in each species"),
                              (142, "3", "Within a centromere the units and the RTE copies are near-identical")],
              panels=[("a", "fig2_panel_a", 0, 8, 88, 52), ("b", "fig2_panel_b", 91, 8, 92, 52),
                      ("c", "fig2_panel_c", 0, 70, 56, 62), ("d", "fig2_panel_d", 59, 70, 66, 62), ("e", "fig2_panel_e", 128, 70, 55, 62),
                      ("f", "fig2_panel_f", 0, 150, 100, 46), ("g", "fig2_panel_g", 103, 150, 38, 46), ("h", "fig2_panel_h", 144, 150, 39, 46)], h=198)


def banner(fig, W, H, top, num, text):
    fig.add_artist(Line2D([0, 1], [1 - top / H, 1 - top / H], color=ST.NAVY, lw=0.5, transform=fig.transFigure))
    fig.text(0.5 / W, 1 - (top + 1.0) / H, num, fontsize=8, fontweight="bold", ha="left", va="top", color=ST.NAVY)
    fig.text(5.0 / W, 1 - (top + 1.0) / H, text, fontsize=8, fontweight="bold", ha="left", va="top", color=ST.NAVY)


def build(n, tag="v2"):
    ST.house_style(); cfg = FIG[n]; W, H = cfg["w"], cfg["h"]; fig = ST.new_figure(W, H); nums = {}
    for top, num, text in cfg["banners"]:
        banner(fig, W, H, top, num, text)
    for L, mod, l, t, w, h in cfg["panels"]:
        try:
            m = importlib.import_module(mod); nums[L] = m.draw(fig, l, t, w, h, L)
        except Exception:
            print("PANEL", L, mod, "FAILED"); traceback.print_exc(); ST.letter(fig, l, t, L)
    out = ROOT / "figures"; out.mkdir(exist_ok=True); ST.save(fig, str(out / f"Figure{n}_{tag}")); return nums


if __name__ == "__main__":
    a = sys.argv[1] if len(sys.argv) > 1 else "both"
    for n in ([1, 2] if a == "both" else [int(a)]):
        build(n); print("figure", n, "done")
