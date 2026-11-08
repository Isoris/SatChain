"""Figure 2, panel b, v3 (2026-11-08: one colour ramp for every row = the RTE colour, so "how much of the RTE each centromere keeps" reads as gold intensity; type I / II is a side column, not a hue; v2 file untouched)
   Figure 2, panel b (house style module): what each MacSat03 centromere keeps of the RTE LINE, and the unit anatomy (92 x 52 mm).
   draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None) -> dict of the numbers shown.
   Standalone:  python fig2_panel_b.py  -> panels/_standalone/fig2_panel_b.{png,pdf,svg}, fig2_panel_b_numbers.tsv, fig2_panel_b_sweep.tsv
DATA (LTR_MacSat03_2026-10-22): split/v46_fragments.tsv (RTE fragment of every v4.6 interruption between two MacSat03 runs: aligned segments in consensus coordinates q, orientation flag),
  split/units_satellite_side.tsv (unit = MacSat03 run + the interruption that follows it, v4.6 catalogue / interruptions_between_runs), insert/Mac.insert_loci.tsv (type II insert loci).
UNIT = satellite side (run start -> next run start). It replaces the old unit set anchored on the transposon-annotation entry (293 units on 12 centromeres).
HEAT STRIP: x = LINE position p = 3,118 - q (5' left), 100-bp bins (the last bin is 17 bp); a bin counts as carried by a unit when >= BIN_COV (50 %) of its bp lie in an aligned RTE segment of that
unit's interruption; colour = fraction of the centromere's units carrying the bin (units without an RTE fragment count as not carrying).  'Kept' bins = fraction >= THR (50 %); 30 / 50 / 70 % and
bin coverage 25 / 50 / 75 % are swept (fig2_panel_b_sweep.tsv).
UNIT ANATOMY: consistent orientation (RTE 3' end toward the preceding satellite block): rows with the RTE written 3'->5' left to right (te_strand '-') are mirrored; medians over units."""
import sys, pathlib, functools
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "code"))
import satellite_style as ST
from matplotlib.patches import Rectangle, FancyArrow, Polygon
from matplotlib.colors import to_rgb

W = pathlib.Path(r"F:\data\satellite_FISH\LTR_MacSat03_2026-10-22")
L = 3117; NB = 32                                          # 31 bins of 100 bp + one of 17 bp
SLATE = "#3d4f66"; PALE = "#d9dcdd"
T2 = {"LG08", "LG16", "LG25"}
BIN_COV, THR = 0.5, 0.5
KEEP = "#8a5a00"                                           # v3: dark ochre, the RTE gold deepened; the ramp runs paper -> KEEP for every row


def ramp(base, f):
    b = np.array(to_rgb(base)); w = np.array(to_rgb("#f3f1ec")); return tuple(w + (b - w) * f)


def _segments(s):
    return [tuple(int(v) for v in x.split("-")) for x in s.split(";")] if isinstance(s, str) and s else []


@functools.lru_cache(maxsize=1)
def data():
    F = pd.read_csv(W / "split" / "v46_fragments.tsv", sep="\t"); U = pd.read_csv(W / "split" / "units_satellite_side.tsv", sep="\t")
    assert (F.chrom.values == U.chrom.values).all() and ((F.start0 > U.unit_lo) & (F.start0 < U.unit_hi)).all()
    F["sat_bp"] = U.sat_bp.values; F["has"] = F.aligned_bp >= 100
    keep = F.groupby("chrom").has.sum(); chroms = list(keep[keep >= 0.5 * F.groupby("chrom").size()].index)            # centromeres whose units carry the RTE
    F = F[F.chrom.isin(chroms)].copy()
    cov = np.zeros((len(F), NB))                                                       # fraction of every LINE-position bin covered, per unit
    for i, r in enumerate(F.itertuples()):
        if not r.has: continue
        m = np.zeros(L + 2, bool)
        for q0, q1 in _segments(r.segments): m[3118 - q1:3118 - q0 + 1] = True           # p = 3118 - q (LINE position 1..3117)
        for k in range(NB): lo, hi = 100 * k + 1, min(100 * (k + 1), L); cov[i, k] = m[lo:hi + 1].mean()
    return F.reset_index(drop=True), cov


def presence(F, cov, bin_cov=BIN_COV):
    out = {}
    for c in sorted(F.chrom.unique()): out[c] = (cov[(F.chrom == c).values] >= bin_cov).mean(0)
    return pd.DataFrame(out).T


def sweep():
    F, cov = data(); ref = presence(F, cov) >= THR; rows = []
    for bc in (0.25, 0.5, 0.75):
        P = presence(F, cov, bc)
        for th in (0.3, 0.5, 0.7):
            K = P >= th
            for c in K.index:
                a, b = K.loc[c].values, ref.loc[c].values
                rows.append(dict(bin_coverage=bc, threshold=th, chrom=c, type="II" if c in T2 else "I", bins_kept=int(a.sum()), bins_kept_ref=int(b.sum()), bins_changed=int((a != b).sum()),
                                 jaccard_vs_ref=round(float((a & b).sum() / max((a | b).sum(), 1)), 3)))
    return pd.DataFrame(rows)


@functools.lru_cache(maxsize=1)
def anatomy():
    F, _ = data(); G = F[F.has].copy(); G["flip"] = G.src_end_left > G.src_end_right
    G["Lf"] = np.where(G.flip, G.right_flank_bp, G.left_flank_bp); G["Rf"] = np.where(G.flip, G.left_flank_bp, G.right_flank_bp); G["t"] = np.where(G.chrom.isin(T2), "II", "I")
    I = pd.read_csv(W / "insert" / "Mac.insert_loci.tsv", sep="\t"); I = I[I.qcov_bp >= 3000]; rows = []
    for r in G[G.t == "II"].itertuples():                                                   # position of the type II insert inside its interruption (consistent orientation)
        for x in I[(I.chrom == "C_mac_" + r.chrom) & (I.hi > r.start0) & (I.lo < r.end)].itertuples():
            a, b = max(x.lo, r.start0) - r.start0, min(x.hi, r.end) - r.start0
            if r.flip: a, b = r.spacer_bp - b, r.spacer_bp - a
            rows.append((r.chrom, a, b))
    D = pd.DataFrame(rows, columns=["chrom", "a", "b"])
    out = {}
    for t in ("I", "II"):
        g = G[G.t == t]; out[t] = dict(n=len(g), sat=float(g.sat_bp.median()), left=float(g.Lf.median()), rte=float(g.frag_bp.median()), right=float(g.Rf.median()), spacer=float(g.spacer_bp.median()))
    out["ins"] = dict(n=len(D), start=float(D.a.median()), end=float(D.b.median()), units_with_insert=int(len(D)), units_II=int((G.t == "II").sum()))
    return out


def draw(fig, left_mm, top_mm, w_mm=92.0, h_mm=52.0, letter=None):
    F, cov = data(); P = presence(F, cov); An = anatomy()
    if letter: ST.letter(fig, left_mm, top_mm, letter)
    ov = ST.ax_mm(fig, left_mm, top_mm, w_mm, h_mm); ov.set_xlim(0, w_mm); ov.set_ylim(h_mm, 0); ov.axis("off"); ov.patch.set_alpha(0)
    ov.text(4.2, 0.2, "What the centromeres keep of the RTE", fontsize=7, fontweight="bold", ha="left", va="top")
    X0, X1 = 13.0, 77.5; sc = (X1 - X0) / L; X = lambda p: X0 + p * sc
    # ---- key (heat) in the title row
    kx, ky = 53.5, 0.9
    for i in range(24):
        ov.add_patch(Rectangle((kx + 0.9 * i * 1.0, ky), 0.9, 2.2, fc=ramp(KEEP, i / 23), ec="none"))
    ov.text(kx - 0.8, ky + 1.1, "0", fontsize=7, ha="right", va="center"); ov.text(kx + 21.6 + 0.8, ky + 1.1, "100 % of units\ncarry this part", fontsize=6, ha="left", va="center", linespacing=0.95)
    # ---- reference row + 15 centromere rows
    PITCH, BH, Y0 = 1.9, 1.55, 5.2
    ov.add_patch(Rectangle((X(0), Y0), X(L) - X(0), BH, fc=ST.RTE, ec="none")); ov.text(X0 - 1.0, Y0 + BH / 2, "RTE", fontsize=6, ha="right", va="center")
    order = sorted(P.index, key=lambda c: (c in T2, -P.loc[c].mean(), c)); y = Y0 + PITCH + 0.5; ytop_rows = y
    nunits = F.groupby("chrom").size()
    ov.text(w_mm - 0.3, Y0 + BH / 2, "units", fontsize=6, ha="right", va="center"); ov.text(w_mm - 7.4, Y0 + BH / 2, "type", fontsize=6, ha="right", va="center")
    first_t2 = True
    for c in order:
        base = KEEP
        if c in T2 and first_t2:
            ov.plot([X0 - 6.0, w_mm - 0.3], [y - 0.25, y - 0.25], color=ST.NAVY, lw=0.4, ls=(0, (1.5, 1.5))); first_t2 = False
        ov.add_patch(Rectangle((X(0), y), X(L) - X(0), BH, fc="#efede8", ec="none", lw=0))
        for k in range(NB):
            f = P.loc[c, k]
            if f > 0: ov.add_patch(Rectangle((X(100 * k), y), X(min(100 * (k + 1), L)) - X(100 * k) + 0.02, BH, fc=ramp(base, f), ec="none", lw=0))
        ov.text(X0 - 1.0, y + BH / 2, c, fontsize=6, ha="right", va="center", color=ST.NAVY)
        ov.text(w_mm - 0.3, y + BH / 2, str(nunits[c]), fontsize=6, ha="right", va="center", color=ST.NAVY)
        ov.text(w_mm - 7.4, y + BH / 2, "II" if c in T2 else "I", fontsize=6, ha="right", va="center", color=ST.NAVY)
        y += PITCH
    yb = y + 0.1
    ov.plot([X(0), X(L)], [yb, yb], color=ST.NAVY, lw=0.5)
    for p in range(0, L + 1, 500):
        ov.plot([X(p)] * 2, [yb, yb + 0.9], color=ST.NAVY, lw=0.5, solid_capstyle="butt")
        ov.text(X(p), yb + 1.2, f"{p:,}" + (" bp" if p == 3000 else ""), fontsize=7, ha="center", va="top")
    # ---- unit anatomy to scale
    ks = 4.7 / 1000.0; UX = X0                                                          # mm per bp
    def bar(x, ln, yy, col, txt=None, tc="white", h=3.1):
        ov.add_patch(Rectangle((x, yy), ln * ks, h, fc=col, ec="none"))
        if txt and ln * ks > 5.0: ov.text(x + ln * ks / 2, yy + h / 2, txt, fontsize=7, ha="center", va="center", color=tc)
        return x + ln * ks
    def rte(x, ln, yy, h=3.1):                                                          # arrow: 5' -> 3' points left (3' end toward the preceding block)
        hl = min(1.6, ln * ks * 0.4)
        ov.add_patch(Polygon([[x + ln * ks, yy], [x + hl, yy], [x, yy + h / 2], [x + hl, yy + h], [x + ln * ks, yy + h]], closed=True, fc=ST.RTE, ec="none"))
        if ln * ks > 8: ov.text(x + ln * ks / 2 + hl / 2, yy + h / 2, "RTE", fontsize=7, ha="center", va="center", color=ST.NAVY)
        elif ln * ks > 5: ov.text(x + ln * ks / 2 + hl / 2, yy + h / 2, "RTE", fontsize=7, ha="center", va="center", color=ST.NAVY)
        return x + ln * ks
    y1 = yb + 5.1; y2 = y1 + 4.3
    a = An["I"]; x = UX
    ov.text(X0 - 1.0, y1 + 1.55, "type I", fontsize=7, ha="right", va="center", color=ST.NAVY)
    x = bar(x, a["sat"], y1, ST.CEN, "block"); x = bar(x, a["left"], y1, PALE, "flank", ST.NAVY); x = rte(x, a["rte"], y1); bar(x, a["right"], y1, PALE, "flank", ST.NAVY)
    b = An["II"]; ins = An["ins"]; x = UX
    ov.text(X0 - 1.0, y2 + 1.55, "type II", fontsize=7, ha="right", va="center", color=ST.NAVY)
    x = bar(x, b["sat"], y2, ST.CEN, "block")
    pre = max(ins["start"], 0.0); ilen = ins["end"] - ins["start"]; post = b["left"] - ins["end"]
    if pre > 50: x = bar(x, pre, y2, PALE)
    x = bar(x, ilen, y2, SLATE, "insert"); x = bar(x, max(post, 0), y2, PALE); x = rte(x, b["rte"], y2)
    if b["right"] > 50: bar(x, b["right"], y2, PALE)
    # scale bar 5 kb
    ys = y2 + 3.1 + 2.0
    ov.plot([UX, UX + 5000 * ks], [ys, ys], color=ST.NAVY, lw=0.5, solid_capstyle="butt")
    for p in (0, 5000): ov.plot([UX + p * ks] * 2, [ys - 0.6, ys + 0.6], color=ST.NAVY, lw=0.5, solid_capstyle="butt")
    ov.text(UX + 5000 * ks + 1.2, ys, "5 kb", fontsize=7, ha="left", va="center")
    # ---- numbers
    K = P >= THR; sw = sweep()
    out = {"centromeres_with_RTE": len(P), "type_I_centromeres": int(sum(c not in T2 for c in P.index)), "type_II_centromeres": int(sum(c in T2 for c in P.index)),
           "units_total": int(len(F)), "units_with_RTE_fragment": int(F.has.sum()), "unit_def": "MacSat03 run + following v4.6 interruption (satellite side)",
           "typeI_unit_n": An["I"]["n"], "typeI_block_bp": An["I"]["sat"], "typeI_left_flank_bp": An["I"]["left"], "typeI_RTE_span_bp": An["I"]["rte"], "typeI_right_flank_bp": An["I"]["right"],
           "typeII_unit_n": An["II"]["n"], "typeII_block_bp": An["II"]["sat"], "typeII_insert_start_in_interruption_bp": An["ins"]["start"], "typeII_insert_end_bp": An["ins"]["end"], "typeII_left_flank_bp": An["II"]["left"],
           "typeII_RTE_span_bp": An["II"]["rte"], "typeII_right_flank_bp": An["II"]["right"], "typeII_units_with_insert_locus": An["ins"]["units_with_insert"], "typeII_units": An["ins"]["units_II"]}
    for c in order: out[f"{c}_units"] = int(nunits[c]); out[f"{c}_bins_kept_50pct"] = int(K.loc[c].sum()); out[f"{c}_bp_kept_50pct"] = int(sum(min(100 * (k + 1), L) - 100 * k for k in np.where(K.loc[c].values)[0]))
    s = sw[(sw.bin_coverage == 0.5)]
    for th in (0.3, 0.5, 0.7): out[f"sweep_thr{int(th * 100)}_max_bins_changed_vs_50"] = int(s[s.threshold == th].bins_changed.max()); out[f"sweep_thr{int(th * 100)}_min_jaccard"] = float(s[s.threshold == th].jaccard_vs_ref.min())
    out["sweep_all9_max_bins_changed"] = int(sw.bins_changed.max()); out["sweep_all9_min_jaccard"] = float(sw.jaccard_vs_ref.min())
    return out


if __name__ == "__main__":
    ST.house_style(); Wm, Hm = 92.0, 52.0
    fig = ST.new_figure(Wm, Hm); res = draw(fig, 0, 0, Wm, Hm, "b")
    ST.save(fig, str(HERE / "_standalone" / "fig2_panel_b_v3"))
    pd.DataFrame({"metric": list(res), "value": list(res.values())}).to_csv(HERE / "fig2_panel_b_v3_numbers.tsv", sep="\t", index=False)
    sweep().to_csv(HERE / "fig2_panel_b_v3_sweep.tsv", sep="\t", index=False)
    P = presence(*data()); P.round(3).to_csv(HERE / "fig2_panel_b_v3_bin_fractions.tsv", sep="\t")
    print({k: v for k, v in res.items() if not k.startswith("LG")})
