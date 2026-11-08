"""Figure 1, panel g (house style module): the 5S rRNA gene mapped onto the 5S-related monomers (55 x 52 mm): GarSat04 (the 5S-derived satellite, SatChain v4.7 sequence cluster C2.2),
MacSat06 (the C. macrocephalus 5S rDNA array, a gene array) and GarSat49 (the C. gariepinus 5S rDNA array, a gene array); GarSat03 / MacSat03 (centromeric group) are the controls.
   draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None) -> dict of the numbers shown.
   Standalone:  python fig1_panel_g.py -> panels/_standalone/fig1_panel_g.{png,pdf}, fig1_panel_g_numbers.tsv
DATA (all computed by code/fig1_g_alignment.py, read here from panels/fig1_panel_g_data/, nothing is re-aligned):
  summary.json (115-bp medoid 5S gene vs the family consensus monomers, local alignment 2/-3/-5/-2, both strands, circular monomer),
  sweep.tsv (72 scoring sets), null.tsv (1,000 mononucleotide shuffles of each monomer), copies_check.tsv (200 genomic monomers per family), best_alignments.tsv (16 genes x consensus).
Orientation: every monomer is drawn in the strand in which the 5S gene reads 5'->3' (GarSat49 is stored as the opposite strand and is flipped here)."""
import sys, json, pathlib
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "code"))
import satellite_style as ST
import fig1_efg_common as C

D = HERE / "fig1_panel_g_data"
ORDER_B = ["GarSat04", "MacSat06", "GarSat49"]; ORDER_A = ["GarSat03", "MacSat03"]        # 5S-related monomers; centromeric controls
WHAT = {"GarSat04": "5S-derived", "MacSat06": "5S rDNA array", "GarSat49": "5S rDNA array"}     # second label line: GarSat04 = 5S-derived satellite; the other two = gene arrays


def _load():
    S = json.load(open(D / "summary.json")); return S, S["detail"]


def _sweep_stats():
    """length of the aligned segment over the 72 scoring sets; 'plausible' = break-even identity mismatch/(match+mismatch) >= 50 % (below that a local alignment extends into random sequence)"""
    s = pd.read_csv(D / "sweep.tsv", sep="\t"); s["breakeven"] = -s.mismatch / (s.match - s.mismatch)
    out = {}
    for lab, g0 in (("all72", s), ("plausible54", s[s.breakeven >= 0.5])):
        for f, g in g0.groupby("family"):
            out[(lab, f)] = dict(n=len(g), len_min=int(g.aln_len.min()), len_med=float(g.aln_len.median()), len_max=int(g.aln_len.max()), id_med=float(g.identity.median()),
                                 n_60_80=int(((g.aln_len >= 60) & (g.aln_len <= 80)).sum()), n_ge_85=int((g.aln_len >= 85).sum()), n_lt_60=int((g.aln_len < 60).sum()),
                                 g_start_med=float(g.gene_start.median()), g_end_med=float(g.gene_end.median()))
    return out


def _copies():
    c = pd.read_csv(D / "copies_check.tsv", sep="\t"); out = {}
    for f, g in c.groupby("family"):
        out[f] = dict(n=len(g), len_med=float(g.aln_len.median()), len_q25=float(g.aln_len.quantile(.25)), len_q75=float(g.aln_len.quantile(.75)), id_med=float(g.identity.median()),
                      frac_60_80=float((g.aln_len.between(60, 80)).mean()), frac_ge_60=float((g.aln_len >= 60).mean()), score_min=float(g.score.min()), score_max=float(g.score.max()))
    return out


def _segments(d):
    """list of (x0, x1) display segments of the aligned region on the monomer (5S sense), splitting at the seam; and the display-frame values"""
    L = d["mono_len"]; m0, m1 = d["m0"], d["m1"]
    if d["strand"] == "-": m0, m1 = L - m1, L - m0
    if m1 > L: return [(m0, L), (0, m1 - L)], True
    return [(m0, m1)], False


def draw(fig, left_mm, top_mm, w_mm=55.0, h_mm=52.0, letter=None):
    S, det = _load(); sw = _sweep_stats(); cp = _copies(); nul = pd.read_csv(D / "null.tsv", sep="\t").set_index("family")
    if letter: ST.letter(fig, left_mm, top_mm, letter)
    ax = C.mm_axes(fig, left_mm, top_mm, w_mm, h_mm)
    C.txt(ax, 4.2, 0.2, "5S gene on the 5S-related monomers", size=7, weight="bold", va="top")
    X0, SC = 15.8, 0.1715                                             # mm, mm per bp (common scale for monomers and gene)
    out = {}; fam_col = {f: ST.family_colour(f) for f in ORDER_B + ORDER_A}
    # ---- monomers (to scale, circular) with the 5S-related segment in gold
    Y0, P = 6.5, 5.15
    for i, f in enumerate(ORDER_B):
        d = det[f]; L = d["mono_len"]; yc = Y0 + i * P
        C.txt(ax, X0 - 1.2, yc - 0.9, f, size=C.FS, ha="right", color=None)
        C.txt(ax, X0 - 1.2, yc + 1.2, WHAT[f], size=6, ha="right", color=ST.MUTED)
        C.rect(ax, X0, yc - 0.75, L * SC, 1.5, fam_col[f], z=2, alpha=0.55)
        segs, wrap = _segments(d)
        for a, b in segs: C.rect(ax, X0 + a * SC, yc - 1.15, (b - a) * SC, 2.3, C.S5, ec=C.GOLD_EDGE, lw=0.35, z=4)
        # circular: tiny seam ticks at both ends of the bar
        for xx in (X0, X0 + L * SC): C.vline(ax, xx, yc - 1.15, yc + 1.15, lw=0.4, z=5)
        n = d["aln_len"]; idp = round(100 * d["ident"]); g0, g1 = d["g0"] + 1, d["g1"]
        full = (f == "GarSat49")
        txt = (f"whole gene, {n} bp, {idp}%" if full else f"{n} bp, {idp}% identity" + (", wraps seam" if wrap else ""))
        C.txt(ax, X0, yc + 2.45, txt, size=C.FS)
        out.update({f"{f}_monomer_bp": L, f"{f}_aligned_bp": n, f"{f}_identity_pct": round(100 * d["ident"], 1), f"{f}_5S_nt_start": g0, f"{f}_5S_nt_end": g1, f"{f}_wraps_seam": wrap,
                    f"{f}_score": d["score"], f"{f}_strand_in_consensus_file": d["strand"]})
    # ---- the gene (115 bp) with the covered interval of each family
    yg = 23.2; Lg = S["ref_len"]
    C.txt(ax, X0 - 1.2, yg, "5S gene", size=C.FS, ha="right")
    C.rect(ax, X0, yg - 0.9, Lg * SC, 1.8, "#F6E7B4", ec=C.GOLD_EDGE, lw=0.35, z=2)
    sx0 = X0 + Lg * SC + 3.0                                              # 50-bp scale bar (common to monomers and gene)
    C.hline(ax, sx0, sx0 + 50 * SC, yg, lw=0.5); C.vline(ax, sx0, yg - 0.5, yg + 0.5, lw=0.5); C.vline(ax, sx0 + 50 * SC, yg - 0.5, yg + 0.5, lw=0.5)
    C.txt(ax, sx0 + 50 * SC + 0.9, yg, "50 bp", size=C.FS)
    for j, f in enumerate(ORDER_B):
        d = det[f]; yy = yg + 2.7 + j * 2.0
        C.rect(ax, X0 + d["g0"] * SC, yy - 0.5, (d["g1"] - d["g0"]) * SC, 1.0, C.S5, ec="none", z=3)
        C.txt(ax, X0 + Lg * SC + 1.0, yy, f, size=6)
    # ---- null comparison
    yn = 32.9
    C.txt(ax, 0.6, yn, "Score vs 1,000 shuffled monomers", size=C.FS, weight="bold")
    NX0, NX1, NMAX = X0, X0 + 36.5, 100.0
    scx = lambda v: NX0 + v / NMAX * (NX1 - NX0)
    rows = ORDER_B + ORDER_A
    for k, f in enumerate(rows):
        yy = yn + 3.2 + k * 2.4; r = nul.loc[f]
        C.txt(ax, X0 - 1.2, yy, f, size=C.FS, ha="right")
        C.rect(ax, NX0, yy - 0.7, scx(r.null_max) - NX0, 1.4, "#D9DCDF", z=1)                         # shuffles: 0 .. maximum of 1,000
        C.vline(ax, scx(r.null_mean), yy - 0.7, yy + 0.7, color=ST.MUTED, lw=0.4, z=2)
        col = ST.FIVE_S if f in ORDER_B else ST.CEN
        over = r.real_score > NMAX
        ax.plot([scx(min(r.real_score, NMAX))], [yy], marker=">" if over else "o", ms=2.8 if over else 2.6, mfc=col, mec="none", zorder=5)
        lab = f"{int(r.real_score)}" + ("" if f in ORDER_B else "  centromeric control")
        if over: C.txt(ax, scx(NMAX) - 1.4, yy, lab, size=C.FS, color=col, ha="right")
        else: C.txt(ax, scx(r.real_score) + 1.2, yy, lab, size=C.FS, color=col)
        out.update({f"{f}_null_score": float(r.real_score), f"{f}_null_mean": float(r.null_mean), f"{f}_null_max": float(r.null_max), f"{f}_null_z": float(r.z), f"{f}_null_shuffles_ge_real": int(r.n_ge_real)})
    ya = yn + 3.2 + len(rows) * 2.4 - 0.7
    C.hline(ax, NX0, NX1, ya, lw=0.5)
    for v in (0, 50, 100):
        C.vline(ax, scx(v), ya, ya + 0.7, lw=0.5); C.txt(ax, scx(v), ya + 1.9, str(v), size=C.FS, ha="center")
    C.txt(ax, X0 - 1.2, ya + 1.9, "score", size=C.FS, ha="right")
    out["sweep"] = {f"{k[0]}|{k[1]}": v for k, v in sw.items()}; out["copies"] = cp
    out["ref_len"] = Lg; out["n_refs"] = S["n_refs"]; out["ref_medoid"] = S["ref_medoid"]
    return out


def numbers_rows(res):
    rows = []
    for k, v in res.items():
        if k in ("sweep", "copies"): continue
        rows.append((k, v, "summary.json / null.tsv (consensus monomer vs 115-bp medoid 5S gene, scoring 2/-3/-5/-2)"))
    for k, v in res["sweep"].items():
        for kk, vv in v.items(): rows.append((f"sweep_{k}_{kk}", round(vv, 3) if isinstance(vv, float) else vv, "sweep.tsv; all72 = 72 scoring sets, plausible54 = break-even identity >= 50%"))
    for k, v in res["copies"].items():
        for kk, vv in v.items(): rows.append((f"copies_{k}_{kk}", round(vv, 3) if isinstance(vv, float) else vv, "copies_check.tsv (200 random genomic monomers per family)"))
    return rows


if __name__ == "__main__":
    Wm, Hm = 55.0, 52.0
    res = C.standalone("fig1_panel_g", draw, Wm, Hm, "g")
    C.write_numbers(C.PANELS / "fig1_panel_g_numbers.tsv", numbers_rows(res))
    for k, v in res.items():
        if k not in ("sweep", "copies"): print(k, v)
    print(json.dumps(res["sweep"], indent=0)[:3000]); print(res["copies"])
