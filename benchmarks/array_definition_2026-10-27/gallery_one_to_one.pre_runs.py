"""Gallery of the step-02 intervals that equal one run-based array (kind 1:1 in compare_<sp>_20k.tsv, >= 2 runs) (2026-10-27).
One row per array, largest first, 30 per page: barcode of its step-03 copies along the array span (x = position in the span), coloured by class
rank (C01 = most copies; U grey), interruptions white, assembly gaps black ticks. Right: chrom and Mb, family, runs, satellite / span kb, classes,
composition domains (domain_seg.py, W = 40 copies, 100 shuffles; '-' when < 100 copies), centromeric flag. NEW arrays have no step-03 copies in
the v4.5 output -> barcode empty. Also writes one_run_arrays_<sp>.tsv.
usage: python -P gallery_one_to_one.py [Mac,Gar]"""
import os, sys, csv, collections, pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
B = os.path.dirname(os.path.abspath(__file__)); R = "F:/data/fClaHyb_func_annotation_v1.0"
sys.path.insert(0, "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/figure/diagnostics"); from domain_seg import domains
PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#8a4fd1", "#e6b400", "#d1495b", "#17becf", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#393b79"]
for sp in (sys.argv[1] if len(sys.argv) > 1 else "Mac,Gar").split(","):
    A = pd.read_csv(f"{B}/compare_{sp}_20k.tsv", sep="\t"); A = A[(A.kind == "1:1") & (A.runs >= 2)].sort_values("bp", ascending=False)
    C = pd.read_csv(f"{R}/fClaHyb_{sp}/fClaHyb_{sp}_satelliteDNA/fClaHyb_{sp}.satDNA_v46.catalogue.tsv", sep="\t")
    G = pd.read_csv(f"{R}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.gaps.detgaps.bed", sep="\t", header=None, usecols=[0, 1, 2], names=["c", "a", "b"])
    need = set(A.array_id); runs = {a: list(zip(g.start0, g.end)) for a, g in C[C.array_id.isin(need)].groupby("array_id")}
    fam = C[C.array_id.isin(need)].groupby("array_id").agg(family=("family_v46", "first"), name=("name_v46", "first"))
    cop = collections.defaultdict(list)
    with open(f"F:/data/satellite_FISH/pipeline_v4_5_fixed/out/fClaHyb_{sp}_satellites/03_monomer_copies.tsv", newline="") as fh:
        rd = csv.reader(fh, delimiter="\t"); h = next(rd); ix = {k: i for i, k in enumerate(h)}
        for row in rd:
            a = row[ix["array_id"]]
            if a not in need: continue
            s, e = int(row[ix["genomic_start_0based"]]), int(row[ix["genomic_end"]])
            if not any(x <= s and e <= y for x, y in runs[a]): continue
            ok = row[ix["complete_or_truncated"]] == "complete" and float(row[ix["best_alignment_identity"]] or 0) >= 80
            cop[a].append((s, e, row[ix["monomer_class"]] if ok else "U"))
    rows = []
    for r in A.itertuples():
        v = sorted(cop.get(r.array_id, [])); ncls = len({c for *_, c in v if c != "U"})
        nd = "-"
        if len([x for x in v if x[2] != "U"]) >= 100:
            gp = [int(g.a) for g in G[(G.c == r.chrom) & (G.a > r.lo) & (G.a < r.hi)].itertuples()]
            D, _ = domains([x for x in v if x[2] != "U"], gp, Ws=(40,), nperm=100); nd = len(D[40]["domains"])
        rows.append(dict(array_id=r.array_id, chrom=r.chrom, start=r.lo, end=r.hi, runs=r.runs, sat_kb=round(r.bp / 1e3, 1), span_kb=round((r.hi - r.lo) / 1e3, 1),
                         family=fam.family.get(r.array_id, ""), name=fam.name.get(r.array_id, ""), copies=len(v), classes=ncls, U=sum(1 for x in v if x[2] == "U"),
                         domains=nd, centromeric=bool(r.cen)))
    T = pd.DataFrame(rows); T.to_csv(f"{B}/one_run_arrays_{sp}.tsv", sep="\t", index=False)
    per = 30
    for p in range(0, len(T), per):
        sub = T.iloc[p:p + per]; fig = plt.figure(figsize=(18, 0.42 * len(sub) + 1.4)); fig.patch.set_facecolor("white")
        ax = fig.add_axes([0.01, 0.4 / (0.42 * len(sub) + 1.4), 0.55, (0.42 * len(sub)) / (0.42 * len(sub) + 1.4)])
        for y, r in enumerate(sub.itertuples()):
            span = max(1, r.end - r.start); yy = len(sub) - 1 - y
            ax.add_patch(Rectangle((0, yy - 0.06), 1, 0.12, color="#efefef", lw=0))
            cnt = collections.Counter(c for *_, c in cop.get(r.array_id, []) if c != "U"); rank = {c: i for i, (c, _) in enumerate(cnt.most_common())}
            for s, e, c in cop.get(r.array_id, []):
                ax.add_patch(Rectangle(((s - r.start) / span, yy - 0.32), max((e - s) / span, 0.0008), 0.64, color="#bbbbbb" if c == "U" else PAL[rank[c] % len(PAL)], lw=0))
            for g in G[(G.c == r.chrom) & (G.a > r.start) & (G.a < r.end)].itertuples(): ax.plot([(g.a - r.start) / span] * 2, [yy - 0.45, yy + 0.45], color="black", lw=1.5)
            ax.text(1.01, yy, f"{r.array_id}  {r.chrom.split('_')[-1]}:{r.start / 1e6:.2f} Mb  {r.family} ({r.name})  runs {r.runs}  sat {r.sat_kb:.0f} / span {r.span_kb:.0f} kb  "
                    f"classes {r.classes}  U {r.U}  domains {r.domains}" + ("  CEN" if r.centromeric else ""), va="center", fontsize=8, color="#c0392b" if r.domains not in ("-", 1) else "#222")
        ax.set_xlim(0, 1); ax.set_ylim(-0.6, len(sub) - 0.4); ax.axis("off")
        fig.suptitle(f"fClaHyb {sp}: step-02 intervals equal to one run-based array (n = {len(T)}), largest first, page {p // per + 1}/{(len(T) - 1) // per + 1}. "
                     "Barcode = copies along the span (colour = class rank, C01 first; grey U; white = interruption; black = assembly gap); red text = > 1 composition domain",
                     x=0.01, ha="left", fontsize=10.5, fontweight="bold")
        fig.savefig(f"{B}/gallery_{sp}_p{p // per + 1:02d}.png", dpi=105); plt.close(fig)
    print(sp, len(T), "arrays;", "domains > 1:", int((T.domains.astype(str).isin(["-", "1"]) == False).sum()), "| no copies (NEW):", int((T.copies == 0).sum()), "| pages:", (len(T) - 1) // per + 1)