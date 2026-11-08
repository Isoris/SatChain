"""Figure 1, panel f (56 x 52 mm): 5S rDNA clusters along the chromosomes of both haplotypes, drawn from coordinates, against the 5S-related arrays
(GarSat04 = the 5S-derived satellite; GarSat49 and MacSat06 = the 5S rDNA gene arrays of the two subgenomes, gene arrays and not satellites; SatChain v4.7 wording, 2026-11-01).
   draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None) -> dict of the numbers shown.
   Standalone:  python fig1_panel_f.py  -> panels/_standalone/fig1_panel_f.{png,pdf}, fig1_panel_f_numbers.tsv
DATA: annotation_overlays/{Mac,Gar}/fClaHyb_{Mac,Gar}.rRNA.sorted.gff3.gz  (barrnap 5S calls, Name=5S_rRNA)
      fClaHyb_func_annotation_v1.0/fClaHyb_{Mac,Gar}/fClaHyb_{Mac,Gar}_genome/*_LG.fa.fai (chromosome sizes)
      satchain_v4.7_2026-11-01/runs_of_record/{Mac,Gar}/runs_of_record.tsv  (5S-related families MacSat06, GarSat04, GarSat49; NB v4.7 has no rows for Mac LG14/17/18 and Gar LG23/27/28)
CLUSTER = 5S calls merged when the gap between neighbours is <= 50 kb, kept when >= 5 calls (the rule of the old overview figure); the gap is swept 10 kb - 1 Mb.
AT A 5S-RELATED ARRAY = the cluster, widened by 10 kb each side, overlaps a GarSat04 / GarSat49 / MacSat06 run of record; swept 0 / 10 / 50 / 100 kb.
Encoding: gold mark = cluster (mark HEIGHT = log10 of the number of 5S calls, width = true extent, min 0.3 mm); solid gold = at a 5S-related array, open gold = elsewhere;
          navy segment on a chromosome = 5S-related array (runs merged at <= 50 kb)."""
import sys, pathlib
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "code"))
import satellite_style as ST
import fig1_efg_common as C
from matplotlib.patches import Rectangle

OV = pathlib.Path(r"F:\data\satellite_FISH\annotation_overlays")
RR = ST.RUNS
AN = pathlib.Path(r"F:\data\fClaHyb_func_annotation_v1.0")
OLD = pathlib.Path(r"F:\data\satellite_FISH\LTR_MacSat03_2026-10-22\figure\diagnostics\markers_rdna_telomere\5S_rDNA\5S_rDNA_sites.tsv")
LB = {"MacSat06", "GarSat04", "GarSat49"}            # the 5S-related arrays: 5S-derived satellite + the two 5S rDNA gene arrays
GAP, MINN, TOL = 50_000, 5, 10_000
CHROM = "#cfd3d6"
SP = {"Mac": "C. macrocephalus", "Gar": "C. gariepinus"}
BLUE = ST.FIVE_S
FS = 6.5


def genes(h):
    g = pd.read_csv(OV / h / f"fClaHyb_{h}.rRNA.sorted.gff3.gz", sep="\t", comment="#", header=None, names="chrom src typ s e sc st ph attr".split())
    g = g[g.attr.str.contains("Name=5S_rRNA")].copy(); g["s"] -= 1
    return g.sort_values(["chrom", "s"]).reset_index(drop=True)


def clusters(g, gap=GAP, minn=MINN):
    rows = []
    for c, x in g.groupby("chrom"):
        cur = None
        for r in x.itertuples():
            if cur and r.s - cur[1] <= gap: cur[1] = max(cur[1], r.e); cur[2] += 1
            else:
                if cur: rows.append([c] + cur)
                cur = [r.s, r.e, 1]
        if cur: rows.append([c] + cur)
    d = pd.DataFrame(rows, columns=["chrom", "s", "e", "n"])
    return d[d.n >= minn].reset_index(drop=True)


def sizes(h):
    d = AN / f"fClaHyb_{h}" / f"fClaHyb_{h}_genome"
    f = next(d.glob("*_LG.fa.fai"))
    z = pd.read_csv(f, sep="\t", header=None, usecols=[0, 1], names=["chrom", "len"])
    z["lg"] = z.chrom.str.extract(r"LG(\d+)")[0].astype(int)
    return z.sort_values("lg").reset_index(drop=True)


def runs(h):
    return pd.read_csv(RR / h / "runs_of_record.tsv", sep="\t")


def lb_arrays(h, merge=50_000):
    t = runs(h); t = t[t.family.isin(LB)].sort_values(["chrom", "start0"]); rows = []
    for c, x in t.groupby("chrom"):
        cur = None
        for r in x.itertuples():
            if cur and r.start0 - cur[1] <= merge: cur[1] = max(cur[1], r.end)
            else:
                if cur: rows.append([c] + cur)
                cur = [r.start0, r.end]
        if cur: rows.append([c] + cur)
    return pd.DataFrame(rows, columns=["chrom", "s", "e"])


def classify(h, d, tol=TOL):
    """B = overlaps a 5S-related run (GarSat04 / GarSat49 / MacSat06, +-tol); S = overlaps another run of record; - = none"""
    t = runs(h); out = []
    for r in d.itertuples():
        x = t[(t.chrom == r.chrom) & (t.end >= r.s - tol) & (t.start0 <= r.e + tol)]
        out.append("B" if x.family.isin(LB).any() else ("S" if len(x) else "-"))
    return out


def kfmt(n):
    return f"{n / 1000:.1f}k" if n >= 1000 else str(int(n))


def draw(fig, left_mm, top_mm, w_mm=56.0, h_mm=52.0, letter=None):
    if letter: ST.letter(fig, left_mm, top_mm, letter)
    ax = C.mm_axes(fig, left_mm, top_mm, w_mm, h_mm)
    C.txt(ax, 4.2, 0.2, "5S rDNA clusters along the chromosomes", size=7, weight="bold", va="top")
    SC = 0.33                                   # mm per Mb, common to both haplotypes
    X0 = {"Mac": 4.2, "Gar": 30.6}; Y0, PITCH = 8.4, 1.2
    out = {}; tab = []
    hh = lambda n: 0.5 + 0.222 * np.log10(max(n, 5) / 5.0)       # mark height (mm)
    for h in ("Mac", "Gar"):
        sz = sizes(h); g = genes(h); d = clusters(g); d["cls"] = classify(h, d); arr = lb_arrays(h)
        n = len(sz); x0 = X0[h]
        C.txt(ax, x0 + 9.0, 5.3, SP[h], size=7, ha="center", color=ST.NAVY, fontstyle="italic")
        top = d.sort_values("n", ascending=False).head(3)
        labrow = {r: [] for r in top.chrom.unique()}
        for r in top.itertuples(): labrow[r.chrom].append(r.n)
        req = {int(c.split("LG")[-1]) for c in labrow}
        shown = set(req)
        for lg in (1, 5, 10, 15, 20, 25):
            if all(abs(lg - q) >= 2 for q in shown): shown.add(lg)
        for i, r in sz.iterrows():
            y = Y0 + (i + 0.5) * PITCH
            C.rect(ax, x0, y - 0.2, r.len / 1e6 * SC, 0.4, CHROM, z=1)
            for a in arr[arr.chrom == r.chrom].itertuples():
                C.rect(ax, x0 + a.s / 1e6 * SC, y - 0.2, max((a.e - a.s) / 1e6 * SC, 0.3), 0.4, BLUE, z=2)
            cl = d[d.chrom == r.chrom]
            for c in cl.itertuples():
                hgt = hh(c.n); xx = x0 + c.s / 1e6 * SC; ww = max((c.e - c.s) / 1e6 * SC, 0.3)
                if c.cls == "B": C.rect(ax, xx, y - hgt / 2, ww, hgt, C.S5, ec=C.GOLD_EDGE, lw=0.2, z=4)
                else: C.rect(ax, xx, y - hgt / 2, ww, hgt, "white", ec=C.GOLD_EDGE, lw=0.35, z=4)
            if r.lg in shown:
                C.txt(ax, x0 - 0.7, y, f"{r.lg}", size=6, ha="right", color=ST.NAVY, weight="bold" if r.lg in req else "normal")
            if r.chrom in labrow:
                C.txt(ax, x0 + r.len / 1e6 * SC + 0.9, y, ", ".join(kfmt(v) for v in sorted(labrow[r.chrom], reverse=True)), size=6, color=C.GOLD_EDGE, weight="bold")
        C.txt(ax, x0 - 0.7, Y0 - 1.0, "LG", size=6, ha="right", va="center", color=ST.MUTED)
        # numbers
        d["kb"] = (d.e - d.s) / 1e3
        nB = int((d.cls == "B").sum()); nS = int((d.cls == "S").sum()); nN = int((d.cls == "-").sum())
        out.update({f"{h}_chromosomes": n, f"{h}_5S_calls": len(g), f"{h}_clusters_ge5_gap50kb": len(d), f"{h}_clusters_at_5S_related_array": nB,
                    f"{h}_clusters_at_other_satellite_run": nS, f"{h}_clusters_no_satellite_run": nN, f"{h}_calls_in_clusters": int(d.n.sum()),
                    f"{h}_chromosomes_with_cluster": int(d.chrom.nunique()), f"{h}_largest_cluster_calls": int(d.n.max()),
                    f"{h}_largest_cluster": "%s:%d-%d" % (d.loc[d.n.idxmax(), "chrom"], d.loc[d.n.idxmax(), "s"], d.loc[d.n.idxmax(), "e"]),
                    f"{h}_largest_cluster_class": d.loc[d.n.idxmax(), "cls"]})
        for r in d.itertuples(): tab.append((h, r.chrom, r.s, r.e, round(r.kb, 1), r.n, r.cls))
    # scale bar + legend (four lines)
    ly = 43.4; LS = 2.25
    xs = X0["Mac"]
    ax.add_patch(Rectangle((xs, ly - 0.55), 1.0, 1.1, fc=C.S5, ec=C.GOLD_EDGE, lw=0.2, zorder=4))
    C.txt(ax, xs + 1.7, ly, f"5S cluster at a 5S-related array: Mac {out['Mac_clusters_at_5S_related_array']}, Gar {out['Gar_clusters_at_5S_related_array']}", size=FS)
    ly2 = ly + LS
    ax.add_patch(Rectangle((xs, ly2 - 0.55), 1.0, 1.1, fc="white", ec=C.GOLD_EDGE, lw=0.35, zorder=4))
    C.txt(ax, xs + 1.7, ly2, f"elsewhere: Mac {out['Mac_clusters_at_other_satellite_run'] + out['Mac_clusters_no_satellite_run']}, Gar {out['Gar_clusters_at_other_satellite_run'] + out['Gar_clusters_no_satellite_run']}", size=FS)
    C.rect(ax, xs + 31.5, ly2 - 0.2, 1.6, 0.4, BLUE, z=2); C.txt(ax, xs + 33.6, ly2, "5S-related array", size=FS)
    ly2b = ly2 + LS
    C.txt(ax, xs, ly2b, "5S-related = GarSat04, GarSat49, MacSat06 (5S rDNA)", size=6)
    ly3 = ly2b + LS
    C.txt(ax, xs, ly3, "height = 5S calls:", size=FS)
    xk = xs + 19.0
    for n_, lab in ((10, "10"), (100, "100"), (1000, "1000")):
        C.rect(ax, xk, ly3 - hh(n_) / 2, 0.5, hh(n_), C.S5, ec=C.GOLD_EDGE, lw=0.2, z=4)
        C.txt(ax, xk + 1.0, ly3, lab, size=FS); xk += 1.0 + 1.25 * len(lab) + 1.3
    xb = X0["Gar"] + 12.0
    C.hline(ax, xb, xb + 10 * SC, ly3, lw=0.5); C.txt(ax, xb + 10 * SC + 0.8, ly3, "10 Mb", size=FS)
    # sweeps and comparison with the old table
    for h in ("Mac", "Gar"):
        g = genes(h); t = runs(h)
        for gap in (10_000, 50_000, 100_000, 500_000, 1_000_000):
            d = clusters(g, gap); top = sorted(d.n.tolist(), reverse=True)[:3]
            out[f"sweep_{h}_gap{gap // 1000}kb_clusters"] = len(d); out[f"sweep_{h}_gap{gap // 1000}kb_calls_in_clusters"] = int(d.n.sum())
            out[f"sweep_{h}_gap{gap // 1000}kb_top3_calls"] = "/".join(map(str, top)); out[f"sweep_{h}_gap{gap // 1000}kb_chromosomes"] = int(d.chrom.nunique())
        d = clusters(g)
        for tol in (0, 10_000, 50_000, 100_000):
            c = classify(h, d, tol); out[f"sweep_{h}_tol{tol // 1000}kb_at_5S_related"] = c.count("B"); out[f"sweep_{h}_tol{tol // 1000}kb_at_other_satellite"] = c.count("S"); out[f"sweep_{h}_tol{tol // 1000}kb_none"] = c.count("-")
        for minn in (3, 5, 10, 20):
            out[f"sweep_{h}_gap50kb_min{minn}_clusters"] = len(clusters(g, GAP, minn))
    # old table
    if OLD.exists():
        o = pd.read_csv(OLD, sep="\t"); nd = 0
        for h in ("Mac", "Gar"):
            d = clusters(genes(h)); oo = o[o.species == h]
            m = d.merge(oo, left_on=["chrom", "s"], right_on=["chrom", "start0"], how="outer", indicator=True)
            same = (m._merge == "both") & (m.n_x == m.copies) & (m.e == m.end)
            out[f"old_table_{h}_clusters"] = len(oo); out[f"old_table_{h}_identical_to_gff_clusters"] = int(same.sum()); nd += int((~same).sum())
        out["old_table_differences"] = nd
    # other satellite run families at 'elsewhere' clusters of Mac
    d = clusters(genes("Mac")); d["cls"] = classify("Mac", d); t = runs("Mac"); fam = {}
    for r in d[d.cls == "S"].itertuples():
        x = t[(t.chrom == r.chrom) & (t.end >= r.s - TOL) & (t.start0 <= r.e + TOL)]
        fam["%s:%d" % (r.chrom.split("_")[-1], r.s)] = "/".join(sorted(set(x.family)))
    out["Mac_other_satellite_families_at_clusters"] = "; ".join(f"{k}={v}" for k, v in fam.items())
    global LAST_TABLE; LAST_TABLE = tab
    return out


if __name__ == "__main__":
    Wm, Hm = 56.0, 52.0
    res = C.standalone("fig1_panel_f", draw, Wm, Hm, "f")
    tab = LAST_TABLE; rows = [(k, v, "") for k, v in res.items()]
    rows += [("cluster", f"{h} {c}:{s}-{e}", f"{kb} kb; {n} 5S calls; class {cl} (B = 5S-related array GarSat04 / GarSat49 / MacSat06, S = other run of record, - = none)") for h, c, s, e, kb, n, cl in tab]
    C.write_numbers(HERE / "fig1_panel_f_numbers.tsv", rows)
    for k, v in res.items(): print(k, v)
