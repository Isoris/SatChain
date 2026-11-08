"""Figure 1, panel e (house style module): THE MacSat03 UNIT CARRIES THE RTE - the composite unit (one MacSat03 block + the interruption that follows it, carrying the RTE LINE fragment)
repeating along C. macrocephalus centromeres (66 x 52 mm).  MacSat03 = the C. macrocephalus member of the GarSat03 / MacSat03 group (the RTE-carrying centromeric families; SatChain v4.7
sequence cluster C2.1; the word "lineage" is dropped, 2026-11-01).
   draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None) -> dict of the numbers shown.
   Standalone:  python fig1_panel_e.py  -> panels/_standalone/fig1_panel_e.{png,pdf}, fig1_panel_e_numbers.tsv
DATA (LTR_MacSat03_2026-10-22/split): units_satellite_side.tsv (unit = MacSat03 run of record + the interruption that follows, v4.6 catalogue; period = unit length, sat_bp = block, spacer_bp = interruption),
  v46_fragments.tsv (RTE fragment of every interruption: position inside the interruption, aligned_bp, te_strand / src_end_left, src_end_right for the orientation), insert/Mac.insert_loci.tsv (type II 6-kb insert loci, qcov >= 3 kb).
  The rows of the two split tables are parallel (asserted).  Centromeres shown = those where >= 50 % of the units carry an RTE fragment (aligned_bp >= 100); selection swept (aligned_bp 50/100/300 x fraction 30/50/70 %).
TOP: two real stretches to scale (type I LG01, type II LG16): first 4 consecutive units whose length is within 10-16 kb and that carry an RTE fragment; arrow = RTE 5'->3' direction taken from the alignment
  (src_end_left > src_end_right = 3' end to the right); unit boundaries = start of each MacSat03 block.  BOTTOM: every unit of the 15 centromeres as a dot (unit length), median tick per centromere, n above."""
import sys, pathlib, functools
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "code"))
import fig1_efg_common as C
import satellite_style as ST

W = pathlib.Path(r"F:\data\satellite_FISH\LTR_MacSat03_2026-10-22")
BLOCK = ST.CEN; PALE = C.PALE; SLATE = "#3d4f66"; TXT = 6.5
AB, FR = 100, 0.5


@functools.lru_cache(maxsize=1)
def tables():
    U = pd.read_csv(W / "split" / "units_satellite_side.tsv", sep="\t"); F = pd.read_csv(W / "split" / "v46_fragments.tsv", sep="\t")
    assert len(U) == len(F) and (F.chrom.values == U.chrom.values).all() and ((F.start0 > U.unit_lo) & (F.start0 < U.unit_hi)).all()
    D = U.copy()
    for c in ("type", "start0", "end", "aligned_bp", "frag_spacer_from", "frag_spacer_to", "src_end_left", "src_end_right", "te_strand", "frag_bp"): D[c] = F[c].values
    D["has"] = D.aligned_bp >= AB
    return D


def select(D, ab=AB, fr=FR):
    h = (D.aligned_bp >= ab).groupby(D.chrom).mean(); return list(h[h >= fr].index)


def sweeps(D):
    rows = []
    for ab in (50, 100, 300):
        for fr in (0.3, 0.5, 0.7):
            ch = select(D, ab, fr); g = D[D.chrom.isin(ch)]
            rows.append(("sweep_selection", f"aligned_bp>={ab}; fraction>={int(fr*100)}%", f"{len(ch)} centromeres, {len(g)} units, type II {int(g.groupby('chrom').type.first().eq('II').sum())}", ""))
    g = D[D.chrom.isin(select(D))]
    for mn in (0, 500, 1000, 2000):
        s = g[g.sat_bp >= mn]
        rows.append(("sweep_min_block_bp", mn, f"{len(s)} units; median unit {s.period.median()/1000:.2f} kb; IQR {s.period.quantile(.25)/1000:.2f}-{s.period.quantile(.75)/1000:.2f} kb", ""))
    return rows


def stretch(g, I, chrom, n=4):
    g = g[g.chrom == chrom].reset_index(drop=True)
    for i in range(len(g) - n + 1):
        s = g.iloc[i:i + n]
        if s.period.between(10000, 16000).all() and s.has.all() and (s.frag_bp < 4000).all(): return s
    raise RuntimeError("no stretch")


def draw(fig, left_mm, top_mm, w_mm=66.0, h_mm=52.0, letter=None):
    D = tables(); chroms = select(D); G = D[D.chrom.isin(chroms)].copy(); G["flip"] = G.src_end_left > G.src_end_right
    I = pd.read_csv(W / "insert" / "Mac.insert_loci.tsv", sep="\t"); I = I[I.qcov_bp >= 3000]
    if letter: ST.letter(fig, left_mm, top_mm, letter)
    ax = C.mm_axes(fig, left_mm, top_mm, w_mm, h_mm)
    C.txt(ax, 4.2, 1.6, "The MacSat03 unit carries the RTE", size=7, weight="bold")
    typ = G.groupby("chrom").type.first()
    X0 = 11.5; SC = 0.80 / 1000.0                       # mm per bp
    BH = 3.4; out = {}
    ex = {"LG01": 9.6, "LG16": 17.8}
    shown = {}
    for c, y in ex.items():
        s = stretch(G, I, c); col = BLOCK if typ[c] == "I" else SLATE; lo = s.unit_lo.iloc[0]
        shown[c] = s
        C.txt(ax, X0 - 1.2, y + BH / 2 - 0.9, c, size=TXT, ha="right", color=BLOCK if typ[c] == "I" else SLATE, weight="bold")
        C.txt(ax, X0 - 1.2, y + BH / 2 + 1.5, "type " + typ[c], size=TXT, ha="right", color=BLOCK if typ[c] == "I" else SLATE)
        for r in s.itertuples():
            xb = X0 + (r.unit_lo - lo) * SC; xs = X0 + (r.start0 - lo) * SC; xe = X0 + (r.unit_hi - lo) * SC
            C.rect(ax, xb, y, xs - xb, BH, BLOCK, z=3)
            C.rect(ax, xs, y, xe - xs, BH, PALE, z=2)
            if typ[c] == "II":
                for x in I[(I.chrom == "C_mac_" + c) & (I.hi > r.start0) & (I.lo < r.end)].itertuples():
                    a, b = max(x.lo, r.start0), min(x.hi, r.end); C.rect(ax, X0 + (a - lo) * SC, y, (b - a) * SC, BH, SLATE, z=3)
            a = r.start0 + r.frag_spacer_from; b = r.start0 + r.frag_spacer_to
            C.arrow(ax, X0 + (a - lo) * SC, X0 + (b - lo) * SC, y + BH / 2, BH * 0.62, ST.RTE, C.GOLD_EDGE, forward=bool(r.flip), head=1.4, lw=0.3, z=5)
            C.vline(ax, xb, y - 1.3, y + BH + 1.3, lw=0.6, z=6)
        xe = X0 + (s.unit_hi.iloc[-1] - lo) * SC; C.vline(ax, xe, y - 1.3, y + BH + 1.3, lw=0.6, z=6)
        out[f"example_{c}_units"] = len(s); out[f"example_{c}_span_kb"] = round((s.unit_hi.iloc[-1] - lo) / 1000, 1)
        out[f"example_{c}_coords"] = f"C_mac_{c}:{int(lo)}-{int(s.unit_hi.iloc[-1])}"
    # unit bracket over the first LG01 unit
    s = shown["LG01"]; lo = s.unit_lo.iloc[0]; xa = X0; xb2 = X0 + (s.unit_hi.iloc[0] - lo) * SC; yb = 6.0
    C.hline(ax, xa, xb2, yb, lw=0.5); C.vline(ax, xa, yb - 0.5, yb + 0.5, lw=0.5); C.vline(ax, xb2, yb - 0.5, yb + 0.5, lw=0.5)
    C.txt(ax, (xa + xb2) / 2, yb - 1.7, "one unit", size=TXT, ha="center")
    # scale bar
    ys = 24.2; C.hline(ax, X0, X0 + 10000 * SC, ys, lw=0.5); C.vline(ax, X0, ys - 0.5, ys + 0.5); C.vline(ax, X0 + 10000 * SC, ys - 0.5, ys + 0.5)
    C.txt(ax, X0 + 10000 * SC + 1.0, ys, "10 kb", size=TXT)
    # key (one row, drawn symbols)
    ky = 24.2; kx = 29.0
    C.rect(ax, kx, ky - 0.9, 1.8, 1.8, BLOCK); C.txt(ax, kx + 2.4, ky, "block", size=TXT)
    C.arrow(ax, kx + 10.2, kx + 12.4, ky, 1.8, ST.RTE, C.GOLD_EDGE, forward=True, head=0.8, lw=0.3); C.txt(ax, kx + 13.0, ky, "RTE", size=TXT)
    C.rect(ax, kx + 19.4, ky - 0.9, 1.8, 1.8, SLATE); C.txt(ax, kx + 21.8, ky, "6-kb insert", size=TXT)
    # ---- summary: every unit as a dot
    order = sorted(chroms, key=lambda c: (typ[c] == "II", c)); n1 = sum(typ[c] == "I" for c in order)
    PX0, PX1 = 11.5, 64.5; pitch = (PX1 - PX0) / (len(order) + 0.6)
    xc = {c: PX0 + pitch * (i + 0.5 + (0.6 if typ[c] == "II" else 0)) for i, c in enumerate(order)}
    YT, YB = 33.6, 45.0; YMAX = 42.0; yy = lambda kb: YB - kb / YMAX * (YB - YT)
    med_all = G.period.median() / 1000
    C.txt(ax, 1.0, 29.0, f"{len(G)} units on {len(order)} candidate centromeres (n above)", size=TXT, ha="left", weight="normal")
    C.hline(ax, PX0 - 0.5, PX1, yy(med_all), color=ST.GRID, lw=0.5, z=1)
    C.hline(ax, PX0 - 0.5, PX1, YB, lw=0.5); C.vline(ax, PX0 - 0.5, YT, YB, lw=0.5)
    for k in (0, 20, 40):
        C.hline(ax, PX0 - 1.4, PX0 - 0.5, yy(k), lw=0.5); C.txt(ax, PX0 - 1.8, yy(k), str(k), size=TXT, ha="right")
    C.txt(ax, 1.6, (YT + YB) / 2 + 3.0, "unit length (kb)", size=TXT, ha="center", rotation=90)
    rng = np.random.RandomState(3)
    for c in order:
        g = G[G.chrom == c]; col = BLOCK if typ[c] == "I" else SLATE
        xs = xc[c] + rng.uniform(-0.8, 0.8, len(g))
        ax.scatter(xs, [yy(v / 1000) for v in g.period], s=1.1, color=col, linewidths=0, zorder=3, alpha=0.8, clip_on=False)
        m = g.period.median() / 1000; C.hline(ax, xc[c] - 1.15, xc[c] + 1.15, yy(m), color=ST.NAVY, lw=0.6, z=4)
        C.txt(ax, xc[c], YB + 1.7, c[2:], size=TXT, ha="center", color=col)
        C.txt(ax, xc[c], YT - 1.3, str(len(g)), size=TXT, ha="center", color=col)
        out[f"{c}_units"] = int(len(g)); out[f"{c}_median_unit_kb"] = round(m, 2)
    C.txt(ax, PX0 - 1.8, YB + 1.7, "LG", size=TXT, ha="right")
    # ---- numbers
    P = G.period / 1000
    out.update({"centromeres": len(order), "type_I_centromeres": int(n1), "type_II_centromeres": len(order) - n1, "units_total": int(len(G)), "units_type_I": int((G.type == "I").sum()), "units_type_II": int((G.type == "II").sum()),
                "units_with_RTE_fragment_aligned_ge100bp": int(G.has.sum()), "unit_length_median_kb": round(float(P.median()), 2), "unit_length_Q1_kb": round(float(P.quantile(.25)), 2), "unit_length_Q3_kb": round(float(P.quantile(.75)), 2),
                "unit_length_min_kb": round(float(P.min()), 2), "unit_length_max_kb": round(float(P.max()), 2), "units_per_centromere_min": int(G.groupby("chrom").size().min()), "units_per_centromere_max": int(G.groupby("chrom").size().max()),
                "units_per_centromere_median": float(G.groupby("chrom").size().median()), "block_median_bp_typeI": float(G[G.type == "I"].sat_bp.median()), "block_median_bp_typeII": float(G[G.type == "II"].sat_bp.median()),
                "interruption_median_bp_typeI": float(G[G.type == "I"].spacer_bp.median()), "interruption_median_bp_typeII": float(G[G.type == "II"].spacer_bp.median()),
                "median_unit_kb_typeI": round(float(G[G.type == "I"].period.median() / 1000), 2), "median_unit_kb_typeII": round(float(G[G.type == "II"].period.median() / 1000), 2),
                "fraction_arrows_pointing_right": round(float(G[G.has].flip.mean()), 3), "centromere_median_unit_range_kb": f"{min(out[f'{c}_median_unit_kb'] for c in order)}-{max(out[f'{c}_median_unit_kb'] for c in order)}",
                "units_over_30kb": int((P > 30).sum()), "units_on_other_centromeres_without_RTE": int((~D.chrom.isin(chroms)).sum())})
    return out


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    Wm, Hm = 66.0, 52.0
    res = C.standalone("fig1_panel_e", draw, Wm, Hm, "e")
    rows = [(k, v, "") for k, v in res.items()] + sweeps(tables())
    C.write_numbers(HERE / "fig1_panel_e_numbers.tsv", rows)
    print({k: v for k, v in res.items() if not k.startswith("LG")})
