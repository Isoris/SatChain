"""Assemble Figure 1 / Figure 2 from the panel modules (LAYOUT.md). Missing panel modules are skipped (grey placeholder box on the draft).
usage: python -I assemble.py [1|2|both]   -> figures/Figure{1,2}_draft.png (500 dpi) + pdf"""
import sys, pathlib, importlib, traceback, json
HERE = pathlib.Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "panels"))
import satellite_style as ST
from matplotlib.patches import Rectangle
FIG = {
 1: dict(size=(183, 205), panels=[("a","fig1_panel_a",0,0,92,56),("b","fig1_panel_b",95,0,88,56),("c","fig1_panel_c",0,60,88,84),("d","fig1_panel_d",91,60,92,84),
                                  ("e","fig1_panel_e",0,148,66,52),("f","fig1_panel_f",69,148,56,52),("g","fig1_panel_g",128,148,55,52)]),
 2: dict(size=(183, 172), panels=[("a","fig2_panel_a",0,0,88,52),("b","fig2_panel_b",91,0,92,52),("c","fig2_panel_c",0,56,56,62),("d","fig2_panel_d",59,56,66,66),
                                  ("e","fig2_panel_e",128,56,55,66),("f","fig2_panel_f",0,126,100,46),("g","fig2_panel_g",103,126,38,46),("h","fig2_panel_h",144,126,39,46)]),
}
def build(n, tag="draft"):
    ST.house_style(); W, H = FIG[n]["size"]; fig = ST.new_figure(W, H); nums = {}
    for L, mod, l, t, w, h in FIG[n]["panels"]:
        try:
            m = importlib.import_module(mod); nums[L] = m.draw(fig, l, t, w, h, L)
        except ModuleNotFoundError as e:
            if e.name != mod: raise
            ST.letter(fig, l, t, L); ax = ST.ax_mm(fig, l, t + 4, w, h - 4); ax.set_xticks([]); ax.set_yticks([]); ax.text(.5, .5, f"{mod}\nmissing", ha="center", va="center", transform=ax.transAxes, color=ST.MUTED)
        except Exception:
            print("PANEL", L, "FAILED"); traceback.print_exc(); ST.letter(fig, l, t, L)
    out = ROOT / "figures"; out.mkdir(exist_ok=True); ST.save(fig, str(out / f"Figure{n}_{tag}")); return nums
if __name__ == "__main__":
    a = sys.argv[1] if len(sys.argv) > 1 else "both"
    for n in ([1, 2] if a == "both" else [int(a)]): build(n); print("figure", n, "done")
