"""Small helpers shared by the Figure 1 panels e, f, g (paper figures 2026-10-27).  Everything is drawn in millimetres inside a panel axes."""
import sys, csv, pathlib
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import satellite_style as ST
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

PANELS = HERE.parent / "panels"
PALE = "#E3E5E8"            # pale grey: the rest of an interruption (other content)
SLATE = "#8E9AA3"           # darker grey: the 6-kb insert of type II units
S5 = ST.RTE                # the 5S rDNA (gold, as asked for panel f); NB: the same gold is RTE in panel e - see the caption stub
GOLD_EDGE = "#9a7410"
FS = 6.5                    # small text (pt); never below 6 pt at final size


def mm_axes(fig, left, top, w, h):
    """one axes covering the panel; data units = mm, y downwards, origin = panel top-left"""
    ax = ST.ax_mm(fig, left, top, w, h)
    ax.set_xlim(0, w); ax.set_ylim(h, 0); ax.axis("off")
    ax.set_zorder(5)
    return ax


def txt(ax, x, y, s, size=FS, ha="left", va="center", weight="normal", color=None, **kw):
    return ax.text(x, y, s, fontsize=size, ha=ha, va=va, fontweight=weight, color=color or ST.NAVY, clip_on=False, **kw)


def rect(ax, x, y, w, h, fc, ec="none", lw=0.0, z=2, **kw):
    ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z, **kw))


def arrow(ax, x0, x1, yc, h, fc, ec, forward=True, head=0.7, lw=0.25, z=4):
    """arrow-shaped block from x0 to x1 (mm) pointing right (forward) or left"""
    L = x1 - x0; hl = min(head, 0.55 * L)
    if forward:
        pts = [(x0, yc - h / 2), (x1 - hl, yc - h / 2), (x1, yc), (x1 - hl, yc + h / 2), (x0, yc + h / 2)]
    else:
        pts = [(x1, yc - h / 2), (x0 + hl, yc - h / 2), (x0, yc), (x0 + hl, yc + h / 2), (x1, yc + h / 2)]
    ax.add_patch(Polygon(pts, closed=True, facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z, joinstyle="miter"))


def hline(ax, x0, x1, y, color=None, lw=0.5, z=3):
    ax.plot([x0, x1], [y, y], color=color or ST.NAVY, lw=lw, solid_capstyle="butt", zorder=z)


def vline(ax, x, y0, y1, color=None, lw=0.5, z=3):
    ax.plot([x, x], [y0, y1], color=color or ST.NAVY, lw=lw, solid_capstyle="butt", zorder=z)


def write_numbers(path, rows, header=("item", "value", "note")):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(header)
        for r in rows: w.writerow(r)


def standalone(modname, drawfn, w, h, letter):
    ST.house_style()
    fig = plt.figure(figsize=(w * ST.MM, h * ST.MM), dpi=500)
    nums = drawfn(fig, 0, 0, w, h, letter)
    out = PANELS / "_standalone"; out.mkdir(exist_ok=True)
    base = str(out / modname)
    fig.savefig(base + ".png", dpi=500); fig.savefig(base + ".pdf")
    return nums
