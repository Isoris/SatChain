"""What sits inside the major centromere arrays besides the satellite (2026-10-27).
Per major centromere (centromeres/index.tsv, 27 Mac + 27 Gar): span = first run start .. last run end; non-satellite part = span minus the
catalogue runs of this array. Overlapped with (each class = union of its features, clipped to the non-satellite part):
  tRNA (tRNAscan), rRNA (barrnap), ncRNA (Rfam all classes; Rfam tRNA/rRNA counted under ncRNA), genes (LiftOn: exon / CDS bp, gene names),
  SINE (AnnoSINE), RTE (genome-wide hits of the RTE source, HSP >= 150 bp), TE whole elements (EDTA feature >= 500 bp with >= 80 % of its
  length inside the non-satellite part; EDTA also labels satellites as TEs, so partial overlaps are not counted), other satellite families
  (catalogue runs of other arrays nested in the span), unannotated = none of these.
Classes are not exclusive (a SINE inside the RTE is counted in both); 'unannotated' = covered by no class.
Outputs: residuals_per_array.tsv, residuals_features.tsv (every ncRNA / tRNA / rRNA / gene feature with its name), residuals.png.
usage: python -P residuals.py"""
import os, glob, gzip, collections, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "F:/data/satellite_FISH/runs_of_record_v4.6_2026-10-27/code"); import palette_v46 as PAL
B = os.path.dirname(os.path.abspath(__file__)); R0 = "F:/data/fClaHyb_func_annotation_v1.0"
def plot_residuals(P):
    cls = ["RTE", "TE whole element (EDTA)", "SINE", "other satellite family", "tRNA", "rRNA", "ncRNA (Rfam)", "gene (exon)"]
    fig, axs = plt.subplots(1, 2, figsize=(18, 8), sharey=False); fig.patch.set_facecolor("white")
    colr = {"RTE": PAL.TE["LINE"], "TE whole element (EDTA)": "#8f8f8f", "SINE": "#8c6a4f", "other satellite family": PAL.CLASS["satellite"], "tRNA": PAL.CLASS["gene array"], "rRNA": PAL.RDNA["45S rDNA"], "ncRNA (Rfam)": "#b5468a", "gene (exon)": "#3d4f66", "unannotated": "#d9d9d9"}
    HATCH = {"tRNA": "////"}   # tRNA gold is almost the RTE (LINE) amber: tRNA also carries a hatch glyph
    for ax, sp in zip(axs, ("Mac", "Gar")):
        d = P[P.species == sp].reset_index(drop=True); y = np.arange(len(d)); left = np.zeros(len(d))
        for k in cls + ["unannotated"]:
            v = d[k + " kb"].values; ax.barh(y, v, left=left, facecolor=colr[k], edgecolor="white" if k in HATCH else "none", hatch=HATCH.get(k), lw=0, label=k, height=0.8); left += v
        ax.set_yticks(y); ax.set_yticklabels([f"{c.split('_')[-1]} {f}" for c, f in zip(d.chrom, d.family)], fontsize=7.5); ax.invert_yaxis()
        ax.set_xlabel("kb inside the non-satellite part of the major centromere (classes can overlap)"); ax.set_title(f"{sp}: what the interruptions contain", loc="left", fontsize=12, fontweight="bold")
        ax.grid(axis="x", alpha=0.3)
    h_, l_ = axs[0].get_legend_handles_labels(); fig.legend(h_, l_, loc="lower center", ncol=5, fontsize=9, frameon=False); fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig(f"{B}/residuals.png", dpi=450); fig.savefig(f"{B}/residuals.png".replace(".png", ".pdf"))
if os.environ.get("PLOT_ONLY"):                                  # redraw from the saved table (colour changes only)
    plot_residuals(pd.read_csv(f"{B}/residuals_per_array.tsv", sep="\t")); raise SystemExit
D0 = "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/figure/diagnostics/centromeres"
I = pd.read_csv(f"{D0}/index.tsv", sep="\t")
RTE = pd.read_csv("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/split/rte_genome_hits_both.tsv", sep="\t"); RTE["a"] = RTE[["sstart", "send"]].min(1) - 1; RTE["b"] = RTE[["sstart", "send"]].max(1)
def find(sp, sub, pat):
    g = sorted(glob.glob(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_{sub}/{pat}")); return g[0] if g else None
def gff(path, keep_types=None):
    out = collections.defaultdict(list)
    if not path: return out
    op = gzip.open if path.endswith(".gz") else open
    for l in op(path, "rt"):
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or (keep_types and f[2] not in keep_types): continue
        at = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        nm = at.get("Name") or at.get("gene_name") or at.get("product") or at.get("rfam_name") or at.get("ID", "")
        out[f[0]].append((int(f[3]) - 1, int(f[4]), f[2], nm, at.get("biotype", at.get("gene_biotype", "")), at.get("classification", "")))
    return out
def union(iv):
    iv = sorted(iv); m = []
    for a, b in iv:
        if m and a <= m[-1][1]: m[-1][1] = max(m[-1][1], b)
        else: m.append([a, b])
    return m
def clip(iv, segs):
    out = []
    for a, b in iv:
        for s, e in segs:
            x, y = max(a, s), min(b, e)
            if y > x: out.append((x, y))
    return union(out)
bp = lambda iv: sum(b - a for a, b in iv)
rows, feats = [], []
for sp in ("Mac", "Gar"):
    C = pd.read_csv(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_satelliteDNA/fClaHyb_{sp}.satDNA_v46.catalogue.tsv", sep="\t")
    src = {"tRNA": gff(find(sp, "tRNA_tRNAscan", "*.sorted.gff3.gz"), {"tRNA", "pseudogene", "tRNA_pseudogene"}),
           "rRNA": gff(find(sp, "rRNA_barrnap", "*.rRNA.sorted.gff3.gz"), {"rRNA"}),
           "ncRNA (Rfam)": gff(find(sp, "ncRNA_rfam", "*.rfam.all_classes.sorted.gff3.gz")),
           "gene (exon)": gff(find(sp, "genes_lifton", "*.lifton.sorted.gff3.gz"), {"exon", "CDS", "gene"}),
           "SINE": gff(find(sp, "repeats_AnnoSINE", "*.SINE.AnnoSINE.annotations.gff3.gz"))}
    te = gff(find(sp, "repeats_EDTA", f"fClaHyb_{sp}.TEanno.sorted.gff3.gz"))
    print(sp, {k: (os.path.basename(find(sp, *{"tRNA": ("tRNA_tRNAscan", "*.sorted.gff3.gz"), "rRNA": ("rRNA_barrnap", "*.rRNA.sorted.gff3.gz"), "ncRNA (Rfam)": ("ncRNA_rfam", "*.rfam.all_classes.sorted.gff3.gz"),
                                                    "gene (exon)": ("genes_lifton", "*.lifton.sorted.gff3.gz"), "SINE": ("repeats_AnnoSINE", "*.SINE.AnnoSINE.annotations.gff3.gz")}[k]) or "MISSING"), sum(len(v) for v in d.values())) for k, d in src.items()})
    for r in I[I.species == sp].itertuples():
        x = C[C.array_id == r.array_id].sort_values("start0"); c = r.chrom; lo, hi = int(x.start0.min()), int(x.end.max())
        runs = union(list(zip(x.start0, x.end))); non = []; cur = lo
        for a, b in runs:
            if a > cur: non.append((cur, a))
            cur = max(cur, b)
        if hi > cur: non.append((cur, hi))
        rec = dict(species=sp, chrom=c, array_id=r.array_id, family=r.family, span_kb=round((hi - lo) / 1e3, 1), satellite_kb=round(bp(runs) / 1e3, 1), non_satellite_kb=round(bp(non) / 1e3, 1))
        cover = {}
        for k, d in src.items():
            ivs = [(a, b) for a, b, *_ in d.get(c, []) if b > lo and a < hi]
            if k == "gene (exon)": ivs = [(a, b) for a, b, t, *_ in d.get(c, []) if b > lo and a < hi and t in ("exon", "CDS")]
            cover[k] = clip(ivs, non)
            for a, b, t, nm, bt, cl in d.get(c, []):
                if b > lo and a < hi and (k != "gene (exon)" or t == "gene") and clip([(a, b)], non):
                    feats.append(dict(species=sp, chrom=c, array_id=r.array_id, cls=k, type=t, name=nm, biotype=bt, start=a, end=b))
        cover["RTE"] = clip([(a, b) for a, b in zip(RTE[RTE.chrom == c].a, RTE[RTE.chrom == c].b) if b > lo and a < hi], non)
        whole = [(a, b) for a, b, t, nm, bt, cl in te.get(c, []) if b > lo and a < hi and b - a >= 500 and bp(clip([(a, b)], non)) >= 0.8 * (b - a)
                 and t not in ("repeat_region", "target_site_duplication", "long_terminal_repeat")]
        cover["TE whole element (EDTA)"] = clip(whole, non)
        others = C[(C.chrom == c) & (C.array_id != r.array_id) & (C.end > lo) & (C.start0 < hi)]
        cover["other satellite family"] = clip(list(zip(others.start0, others.end)), non)
        anyc = union([iv for v in cover.values() for iv in v])
        for k, v in cover.items(): rec[k + " kb"] = round(bp(v) / 1e3, 2)
        rec["unannotated kb"] = round((bp(non) - bp(anyc)) / 1e3, 1); rec["unannotated share"] = round((bp(non) - bp(anyc)) / max(1, bp(non)), 2)
        rows.append(rec)
P = pd.DataFrame(rows); P.to_csv(f"{B}/residuals_per_array.tsv", sep="\t", index=False); F = pd.DataFrame(feats); F.to_csv(f"{B}/residuals_features.tsv", sep="\t", index=False)
pd.set_option("display.width", 300); pd.set_option("display.max_columns", 30)
print(P.drop(columns=["array_id"]).to_string(index=False))
print("\nfeatures inside the non-satellite part (counts by class / type):"); print(F.groupby(["species", "cls", "type"]).size().to_string() if len(F) else "none")
if len(F): print("\nnamed ncRNA / tRNA / rRNA / genes:"); print(F[F.cls != "SINE"].groupby(["species", "cls", "name"]).size().sort_values(ascending=False).head(40).to_string())
plot_residuals(P)
