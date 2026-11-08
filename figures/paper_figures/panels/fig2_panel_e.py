"""Figure 2, panel e (house style module): positions of the ~6-kb type II insert along the chromosomes of both haplotypes (55 x 66 mm).
   draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None) -> dict of the numbers shown.
   Standalone:  python fig2_panel_e.py  -> panels/_standalone/fig2_panel_e.{png,pdf,svg}, fig2_panel_e_numbers.tsv, fig2_panel_e_sweep.tsv
DATA: LTR_MacSat03_2026-10-22/insert/{Mac,Gar}.insert_hits.tsv (dc-megablast of one representative insert against each haplotype, e <= 1e-20) RECHAINED HERE into loci exactly as code/11_insert.sh
  (hits within 2 kb merged; locus size = union of the covered insert positions, qcov_bp) and compared with insert/{Mac,Gar}.insert_loci.tsv (identical). Large = qcov >= 3 kb of the 6-kb insert,
  fragment = 1-3 kb. Chromosome sizes: fClaHyb_func_annotation_v1.0/fClaHyb_{Mac,Gar}/..._LG.chrom.sizes; centromere intervals: runs_of_record_v4.5.3_2026-10-18/{Mac,Gar}/centromere_check.tsv.
Sweep: large-locus threshold 2 / 3 / 4 / 5 kb; fragment lower bound 0.5 / 1 / 2 kb."""
import sys, pathlib, functools
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "code"))
import satellite_style as ST
from matplotlib.patches import Rectangle

W = pathlib.Path(r"F:\data\satellite_FISH\LTR_MacSat03_2026-10-22")
R = pathlib.Path(r"F:\data\satellite_FISH\runs_of_record_v4.5.3_2026-10-18"); A = pathlib.Path(r"F:\data\fClaHyb_func_annotation_v1.0")
SLATE = "#3d4f66"; SLATE_LIGHT = "#98a6b6"; CHROM = "#cfd3d6"
SP = {"Mac": "C. macrocephalus", "Gar": "C. gariepinus"}


def _qcov(q):
    q = sorted(q); tot = 0; s = e = None
    for a, b in q:
        if e is None or a > e:
            if e is not None: tot += e - s
            s, e = a, b
        else: e = max(e, b)
    return tot + (e - s)


@functools.lru_cache(maxsize=2)
def loci(h):
    t = pd.read_csv(W / "insert" / f"{h}.insert_hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
    t["lo"] = t[["s0", "s1"]].min(axis=1); t["hi"] = t[["s0", "s1"]].max(axis=1); rows = []
    for c, g in t.sort_values("lo").groupby("chrom"):
        cur = None
        for r in g.itertuples():
            if cur and r.lo - cur["hi"] <= 2000: cur["hi"] = max(cur["hi"], r.hi); cur["q"].append((r.q0, r.q1)); cur["pid"].append(r.pid)
            else:
                if cur: rows.append(cur)
                cur = dict(chrom=c, lo=r.lo, hi=r.hi, q=[(r.q0, r.q1)], pid=[r.pid])
        if cur: rows.append(cur)
    L = pd.DataFrame(rows); L["qcov_bp"] = L.q.apply(_qcov); L["pid"] = L.pid.apply(np.mean); L = L.drop(columns="q")
    old = pd.read_csv(W / "insert" / f"{h}.insert_loci.tsv", sep="\t")
    same = len(L) == len(old) and (L.qcov_bp.values == old.qcov_bp.values).all() and (L.lo.values == old.lo.values).all()
    return L, bool(same)


@functools.lru_cache(maxsize=2)
def genome(h):
    sz = pd.read_csv(A / f"fClaHyb_{h}" / f"fClaHyb_{h}_genome" / f"fClaHyb_{h}_LG.chrom.sizes", sep="\t", header=None, names=["chrom", "len"])
    sz = sz[sz.chrom.str.contains("_LG")].sort_values("chrom").reset_index(drop=True)
    cen = pd.read_csv(R / h / "centromere_check.tsv", sep="\t").set_index("chrom")
    return sz, cen


def counts(h, big=3000, lo=1000):
    L, _ = loci(h); return int((L.qcov_bp >= big).sum()), int(((L.qcov_bp >= lo) & (L.qcov_bp < big)).sum())


def draw(fig, left_mm, top_mm, w_mm=55.0, h_mm=66.0, letter=None):
    if letter: ST.letter(fig, left_mm, top_mm, letter)
    ov = ST.ax_mm(fig, left_mm, top_mm, w_mm, h_mm); ov.set_xlim(0, w_mm); ov.set_ylim(h_mm, 0); ov.axis("off"); ov.patch.set_alpha(0)
    ov.text(4.2, 0.2, "Type II insert along the chromosomes", fontsize=7, fontweight="bold", ha="left", va="top")
    sc = 0.38                                                                              # mm per Mb, common to both haplotypes
    COLX = {"Mac": 6.6, "Gar": 32.0}; Y0, PITCH = 8.6, 1.70
    nL = {}; out = {}
    for h in ("Mac", "Gar"):
        sz, cen = genome(h); L, same = loci(h); n = len(sz); nL[h] = n
        out[f"{h}_loci_rechained_equal_repo_table"] = same; out[f"{h}_chromosomes"] = n
        ov.text(COLX[h] + 11.0, 4.4, SP[h], fontsize=7, fontstyle="italic", ha="center", va="top")
        big = L[L.qcov_bp >= 3000]; frag = L[(L.qcov_bp >= 1000) & (L.qcov_bp < 3000)]
        for i, r in sz.iterrows():
            y = Y0 + i * PITCH + PITCH / 2; x0 = COLX[h]; num = int(r.chrom.split("LG")[-1])
            ov.add_patch(Rectangle((x0, y - 0.28), r.len / 1e6 * sc, 0.56, fc=CHROM, ec="none"))
            if r.chrom in cen.index:
                c0, c1 = cen.loc[r.chrom, "cen_start"], cen.loc[r.chrom, "cen_end"]
                ov.add_patch(Rectangle((x0 + c0 / 1e6 * sc, y - 0.62), max((c1 - c0) / 1e6 * sc, 0.5), 1.24, fc=ST.CENTROMERE, ec=ST.CEN_EDGE, lw=0.3, zorder=2))
            for l in frag[frag.chrom == r.chrom].itertuples():
                ov.add_patch(Rectangle((x0 + (l.lo + l.hi) / 2e6 * sc - 0.15, y - 0.62), 0.3, 1.24, fc=SLATE_LIGHT, ec="none", zorder=3))
            b = big[big.chrom == r.chrom]
            for l in b.itertuples():
                ov.add_patch(Rectangle((x0 + l.lo / 1e6 * sc, y - 0.75), max((l.hi - l.lo) / 1e6 * sc, 0.4), 1.5, fc=SLATE, ec="none", zorder=4))
            if len(b): ov.text(x0 + r.len / 1e6 * sc + 0.7, y, str(len(b)), fontsize=7, fontweight="bold", color=SLATE, ha="left", va="center")
            if h == "Mac" and num in (1, 5, 8, 12, 16, 20, 25):
                ov.text(COLX["Mac"] - 0.8, y, f"LG{num:02d}", fontsize=6, ha="right", va="center", color=SLATE if len(b) else ST.NAVY)
        # scale axis
        yb = Y0 + n * PITCH + 0.6 if h == "Mac" else Y0 + 28 * PITCH + 0.6
    yb = Y0 + 28 * PITCH + 0.6
    for h in ("Mac", "Gar"):
        x0 = COLX[h]; ov.plot([x0, x0 + 52 * sc], [yb, yb], color=ST.NAVY, lw=0.5)
        for m in (0, 20, 40):
            ov.plot([x0 + m * sc] * 2, [yb, yb + 0.9], color=ST.NAVY, lw=0.5, solid_capstyle="butt")
            ov.text(x0 + m * sc, yb + 1.2, f"{m}" + (" Mb" if m == 40 else ""), fontsize=7, ha="center", va="top")
    # legend (two rows, drawn symbols)
    nb, nf = counts("Mac"); gb, gf = counts("Gar"); ly = yb + 5.3
    ov.add_patch(Rectangle((4.2, ly - 0.75), 0.6, 1.5, fc=SLATE, ec="none")); ov.text(5.8, ly, f"insert \u2265 3 kb: Mac {nb}, Gar {gb}", fontsize=7, ha="left", va="center")
    ov.add_patch(Rectangle((4.2, ly + 2.7 - 0.75), 0.6, 1.5, fc=SLATE_LIGHT, ec="none")); ov.text(5.8, ly + 2.7, f"fragment 1\u20133 kb: Mac {nf}, Gar {gf}", fontsize=7, ha="left", va="center")
    ov.add_patch(Rectangle((38.4, ly - 0.62), 1.4, 1.24, fc=ST.CENTROMERE, ec=ST.CEN_EDGE, lw=0.3)); ov.text(40.4, ly, "centromere", fontsize=7, ha="left", va="center")
    # numbers / sweep
    sw = sweep()
    for h in ("Mac", "Gar"):
        L, _ = loci(h); big = L[L.qcov_bp >= 3000]
        out.update({f"{h}_large_loci_ge3kb": int(len(big)), f"{h}_fragments_1to3kb": int(((L.qcov_bp >= 1000) & (L.qcov_bp < 3000)).sum()), f"{h}_loci_total": int(len(L)), f"{h}_largest_locus_bp": int(L.qcov_bp.max()),
                    f"{h}_large_in_centromere_pm1Mb": int(pd.read_csv(W / "insert" / f"{h}.insert_loci.tsv", sep="	").query("qcov_bp >= 3000").in_cen.sum()), f"{h}_large_median_identity_pct": round(float(big.pid.median()), 1) if len(big) else ""})
        for c, v in big.chrom.str.split("_").str[-1].value_counts().sort_index().items(): out[f"{h}_large_{c}"] = int(v)
    for r in sw.itertuples():
        if r.chromosomes.startswith("fragment"): out[f"sweep_{r.haplotype}_fragments_lower_bound_{r.chromosomes.split()[-2]}bp"] = r.n_fragment_1kb_to_large
        else: out[f"sweep_{r.haplotype}_large_ge{r.large_kb}kb"] = r.n_large
    return out


def sweep():
    rows = []
    for h in ("Mac", "Gar"):
        L, _ = loci(h)
        for kb in (2, 3, 4, 5):
            b = L[L.qcov_bp >= kb * 1000]
            rows.append(dict(haplotype=h, large_kb=kb, n_large=int(len(b)), chromosomes=";".join(f"{c.split('_')[-1]}:{v}" for c, v in b.chrom.value_counts().sort_index().items()),
                             n_fragment_1kb_to_large=int(((L.qcov_bp >= 1000) & (L.qcov_bp < kb * 1000)).sum())))
        for lo_ in (500, 1000, 2000):
            rows.append(dict(haplotype=h, large_kb=3, n_large=int((L.qcov_bp >= 3000).sum()), chromosomes=f"fragment lower bound {lo_} bp", n_fragment_1kb_to_large=int(((L.qcov_bp >= lo_) & (L.qcov_bp < 3000)).sum())))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    ST.house_style(); Wm, Hm = 55.0, 66.0
    fig = ST.new_figure(Wm, Hm); res = draw(fig, 0, 0, Wm, Hm, "e")
    ST.save(fig, str(HERE / "_standalone" / "fig2_panel_e"))
    pd.DataFrame({"metric": list(res), "value": list(res.values())}).to_csv(HERE / "fig2_panel_e_numbers.tsv", sep="\t", index=False)
    sweep().to_csv(HERE / "fig2_panel_e_sweep.tsv", sep="\t", index=False)
    print(res)
