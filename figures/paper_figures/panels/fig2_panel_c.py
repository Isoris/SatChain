"""Figure 2, panel c (house style module): RTE core copies inside / outside the centromere, by species and family group of the centromere (56 x 62 mm).
   draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None) -> dict of the numbers shown.   Standalone:  python fig2_panel_c.py  -> panels/_standalone/fig2_panel_c.{png,pdf} + fig2_panel_c_numbers.tsv
Data code = code/rte_core_copies_v47.py (chaining of the core BLAST hits into copies as code/06b_core_copies.py; SatChain v4.7 centromere intervals and runs of record, 2026-11-01).
Groups of a centromere (by the v4.7 sequence cluster of its dominant family): cen = the GarSat03 / MacSat03 group (RTE-carrying centromeric families, C2.1); 5S = GarSat04 (5S-derived satellite,
C2.2); other = everything else (gene-array-dominated such as the MacSat06 5S rDNA array on LG22, unsupported, or no v4.7 rows).  The word "lineage" is dropped.
Reference setting: copy covers >= 70 % of the 1,841-bp core (span <= 1.5 x core); inside = copy overlaps the centromere interval +- 1 Mb.
Sweep: window +-0.5 / 1 / 2 Mb x core coverage >= 50 / 70 / 90 % (9 settings)."""
import sys, pathlib, functools
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "code"))
import satellite_style as ST
import rte_core_copies_v47 as RC
import matplotlib.colors as mc
from matplotlib.patches import Rectangle

OLIVE = "#8a9a5b"                                     # other family group (palette_v46.OTHER_SAT; same colour in panel d)
ORDER0 = [("Mac", "cen"), ("Gar", "cen"), ("Gar", "5S"), ("Mac", "5S"), ("Mac", "other"), ("Gar", "other")]
LCOL = {"cen": ST.CEN, "5S": ST.FIVE_S, "other": OLIVE}
NAME = {"cen": "GarSat03 / MacSat03 group", "5S": "GarSat04 (5S-derived)", "other": "other"}


def pale(c, f=0.72):
    r = np.array(mc.to_rgb(c)); return tuple(r + (1 - r) * f)


@functools.lru_cache(maxsize=1)
def table():
    """all rows of the sweep (haplotype x group x coverage x window) = fig2_panel_c_sweep_table.tsv"""
    T, _ = RC.table(); return T


def draw(fig, left_mm, top_mm, w_mm=56.0, h_mm=62.0, letter=None):
    T = table()
    ref = T[(T.cov_min == RC.REF_COV) & (T.window_Mb == RC.REF_PAD / 1e6)].set_index(["haplotype", "centromere_group"])
    ORDER = [k for k in ORDER0 if ref.loc[k, "core_copies"] > 0]            # groups without a chromosome / copy are not drawn (reported in the numbers)
    if letter: ST.letter(fig, left_mm, top_mm, letter)
    ov = ST.ax_mm(fig, left_mm, top_mm, w_mm, h_mm); ov.set_xlim(0, w_mm); ov.set_ylim(h_mm, 0); ov.axis("off"); ov.patch.set_alpha(0)
    ov.text(4.2, 0.2, "RTE core copies in the centromere", fontsize=7, fontweight="bold", ha="left", va="top")
    # key: inside / outside
    for i, (lab, col) in enumerate((("inside ±1 Mb", LCOL["cen"]), ("outside", pale(LCOL["cen"])))):
        x0 = 4.2 + i * 19.5
        ov.add_patch(Rectangle((x0, 4.15), 2.2, 2.2, fc=col, ec="none")); ov.text(x0 + 3.0, 5.25, lab, fontsize=7, ha="left", va="center")
    # ---- stacked bars
    AX_L, AX_W, AX_T, AX_H = 11.5, 44.5, 7.8, 19.0
    ax = ST.ax_mm(fig, left_mm + AX_L, top_mm + AX_T, AX_W, AX_H); ST.hairline_axes(ax)
    pos = {}; cur = 0.0                                                      # bars side by side, a small gap before the "other" groups
    for k in ORDER:
        if k[1] == "other" and not any(q[1] == "other" for q in pos): cur += 0.1
        pos[k] = cur; cur += 1
    x = np.array([pos[k] for k in ORDER]); XMAX = cur - 0.9; out = {}
    for k in ORDER0:
        r = ref.loc[k]; kk = f"{k[0]}_{k[1]}"
        out.update({f"{kk}_core_copies": int(r.core_copies), f"{kk}_inside": int(r.inside), f"{kk}_outside": int(r.outside), f"{kk}_pct_inside": round(float(r.pct_inside), 1) if r.core_copies else "",
                    f"{kk}_n_chromosomes": int(r.n_chromosomes), f"{kk}_chromosomes": r.chromosomes, f"{kk}_pct_next_to_cen": round(float(r.pct_next_to_cen), 1) if r.core_copies else ""})
    for xi, (h, l) in zip(x, ORDER):
        r = ref.loc[(h, l)]; n = int(r.core_copies); a = int(r.inside); b = int(r.outside); pa, pb = 100 * a / n, 100 * b / n; c = LCOL[l]
        ax.bar(xi, pa, 0.78, color=c, lw=0, zorder=2); ax.bar(xi, pb, 0.78, bottom=pa, color=pale(c), lw=0, zorder=2)
        if pa >= 30: ax.text(xi, pa / 2, f"{a}\n{pa:.0f}%", ha="center", va="center", color="white", fontsize=6.5, linespacing=0.95, zorder=3)
        else: ax.text(xi, pa + 3, f"{a}\n{pa:.0f}%", ha="center", va="bottom", color=ST.NAVY, fontsize=6.5, linespacing=0.95, zorder=3)
        ax.text(xi, 98.5, f"{b}\n{pb:.0f}%", ha="center", va="top", color=ST.NAVY, fontsize=6.5, linespacing=0.95, zorder=3)
    ax.set_xlim(-0.5, XMAX + 0.5); ax.set_ylim(0, 100); ax.set_yticks([0, 50, 100]); ax.set_ylabel("% of core copies", labelpad=1)
    ax.set_xticks(x); ax.set_xticklabels([f"{h}\n{l}" for h, l in ORDER], linespacing=0.95, fontsize=7); ax.tick_params(axis="x", length=0, pad=1.2)
    # ---- table rows under the bars (x in mm in the overlay)
    xm = lambda xi: AX_L + AX_W * (xi + 0.5) / (XMAX + 1.0)
    rows = [("n chr.", lambda r: f"{int(r.n_chromosomes)}"), ("copies", lambda r: f"{int(r.core_copies):,}"), ("% near cen", lambda r: f"{r.pct_next_to_cen:.0f}")]
    for k, (lab, fn) in enumerate(rows):
        yy = 35.2 + 3.1 * k
        ov.text(0, yy, lab, ha="left", va="center", fontsize=7)
        for xi, key in zip(x, ORDER): ov.text(xm(xi), yy, fn(ref.loc[key]), ha="center", va="center", fontsize=7)
    # group key (two lines; every colour named)
    ov.text(0, 44.3, "cen = GarSat03 / MacSat03 group (RTE-carrying)", fontsize=6.5, ha="left", va="center")
    ov.text(0, 46.4, "5S = GarSat04 (5S-derived); other = other family", fontsize=6.5, ha="left", va="center")
    # ---- sweep small multiples (window)
    SW_T, SW_H, SW_W, GAP = 51.2, 6.0, 13.6, 1.65
    for kk, pad in enumerate((0.5, 1.0, 2.0)):
        a2 = ST.ax_mm(fig, left_mm + AX_L + kk * (SW_W + GAP), top_mm + SW_T, SW_W, SW_H); ST.hairline_axes(a2)
        for (h, l), xi in pos.items():
            v = T[(T.haplotype == h) & (T.centromere_group == l) & (T.window_Mb == pad)].sort_values("cov_min"); c = LCOL[l]
            if not (v.core_copies > 0).all(): continue
            lo, hi = v.pct_inside.min(), v.pct_inside.max(); r = v[v.cov_min == 0.7].iloc[0]
            a2.plot([xi, xi], [lo, hi], color=c, lw=0.9, solid_capstyle="butt", zorder=2)
            a2.plot(xi, r.pct_inside, marker="o" if h == "Mac" else "s", ms=2.4, mfc=c, mec="white", mew=0.25, ls="", zorder=3)
            k = f"{h}_{l}_pm{pad:g}Mb"; out.update({f"sweep_{k}_ref_pct": round(float(r.pct_inside), 1), f"sweep_{k}_min_pct": round(float(lo), 1), f"sweep_{k}_max_pct": round(float(hi), 1)})
        a2.set_xlim(-0.7, XMAX + 0.7); a2.set_ylim(-8, 108); a2.set_xticks([]); a2.set_yticks([0, 50, 100]); a2.tick_params(length=1.8, pad=1)
        if kk: a2.set_yticklabels([])
        a2.axhline(50, color=ST.GRID, lw=0.4, zorder=0)
        ov.text(AX_L + kk * (SW_W + GAP) + SW_W / 2, SW_T - 0.9, f"±{pad:g} Mb", fontsize=7, ha="center", va="bottom")
        if kk == 0: ov.text(0, SW_T - 0.9, "% inside", fontsize=7, ha="left", va="bottom")
    # key of the sweep markers, drawn as markers
    ky = 60.4
    ov.plot([0.9], [ky], marker="o", ms=2.4, mfc=ST.NAVY, mec="none", ls=""); ov.text(2.3, ky, "Mac", fontsize=7, va="center")
    ov.plot([9.0], [ky], marker="s", ms=2.4, mfc=ST.NAVY, mec="none", ls=""); ov.text(10.4, ky, "Gar", fontsize=7, va="center")
    ov.text(16.5, ky, "bar = 50–90 % core coverage", fontsize=7, va="center")
    out["groups_drawn"] = ",".join(f"{h}_{l}" for h, l in ORDER); out["centromere_source"] = f"{ST.RUNS_VERSION}: runs_of_record/{{Mac,Gar}}/centromere_check.tsv (cen_start0, cen_end, dominant_family)"
    for h in ("Mac", "Gar"):
        out[f"{h}_total_core_copies"] = int(sum(ref.loc[(h, l), "core_copies"] for l in RC.GROUPS)); out[f"{h}_total_inside"] = int(sum(ref.loc[(h, l), "inside"] for l in RC.GROUPS))
    return out


if __name__ == "__main__":
    ST.house_style(); W, H = 56.0, 62.0
    fig = ST.new_figure(W, H); res = draw(fig, 0, 0, W, H, "c")
    ST.save(fig, str(HERE / "_standalone" / "fig2_panel_c"))
    pd.DataFrame({"metric": list(res), "value": list(res.values())}).to_csv(HERE / "fig2_panel_c_numbers.tsv", sep="\t", index=False)
    table().to_csv(HERE / "fig2_panel_c_sweep_table.tsv", sep="\t", index=False, float_format="%.2f")
    print({k: v for k, v in res.items() if "total" in k or "groups" in k or k.endswith("_core_copies") or k.endswith("_inside") or k.endswith("_n_chromosomes")})
