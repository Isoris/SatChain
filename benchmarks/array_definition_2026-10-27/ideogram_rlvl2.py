"""Resolution level 2 of the array-type ideogram (rlvl2; 2026-10-27): micro / mini-satellites resolved by motif and purity from the sequence.
Level 1 (ideogram_types.py) = type from catalogue unit_bp, relative span, k-mer homogeneity, gene flag. Level 2 keeps those types except
'micro / mini-satellite', which is replaced, for each multi-run interval with catalogue unit <= 60 bp, by:
  best motif  = the 1-6 bp motif whose tandem runs (>= 3 copies and >= 12 bp) cover most of the runs' sequence; canonical form = primitive motif
                (ATAT -> AT), smallest rotation of it and of its reverse complement ((CA)n = (TG)n = (GT)n)
  purity      = bp covered by tandem runs of the best motif / bp of the runs
  unit identity = share of positions i with s[i] == s[i + unit] (catalogue unit), over the runs
Labels: 'micro (MOTIF)n' if purity >= 0.5 (top motifs own colour, others 'micro, other motif'); else 'mini 7-15 / 16-30 / 31-60 bp' if unit
identity >= 0.7; else 'mixed / not a clean tandem'.
Outputs: array_types_rlvl2_<sp>.tsv, ideogram_types_<sp>_rlvl2.png (and ideogram_types_<sp>_rlvl1.png = the level-1 figure under the new name).
usage: python -P ideogram_rlvl2.py"""
import os, re, sys, shutil, collections, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
B = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, B); R0 = "F:/data/fClaHyb_func_annotation_v1.0"
from kmer_homogeneity import Fasta
BS = chr(92); RCT = str.maketrans("ACGT", "TGCA")
def canon(m):
    for p in range(1, len(m) + 1):                                              # primitive motif
        if len(m) % p == 0 and m[:p] * (len(m) // p) == m: m = m[:p]; break
    rc = m.translate(RCT)[::-1]
    return min([m[i:] + m[:i] for i in range(len(m))] + [rc[i:] + rc[:i] for i in range(len(rc))])
def micro_profile(s):
    cov = {}                                                                     # per canonical motif: union of covered positions (no double count across k)
    for k in range(1, 7):
        for x in re.finditer(("(([ACGT]{%d})" + BS + "2{2,})") % k, s):
            if x.end() - x.start() >= 12:
                m = canon(x.group(2)); v = cov.setdefault(m, np.zeros(len(s), dtype=bool)); v[x.start():x.end()] = True
    if not cov: return "", 0.0
    m, v = max(cov.items(), key=lambda t: t[1].sum()); return m, float(v.sum()) / max(1, len(s))
L1 = {"gene array": "#7f7f7f", "colinear, homogeneous": "#1f4e9c", "colinear, multi-variant": "#6aa1e0", "composite periodic": "#c0392b",
      "interrupted, irregular": "#e67e22", "scattered": "#b39ddb"}
for sp in ("Mac", "Gar"):
    if os.path.exists(f"{B}/ideogram_types_{sp}.png"): shutil.copyfile(f"{B}/ideogram_types_{sp}.png", f"{B}/ideogram_types_{sp}_rlvl1.png")
    T = pd.read_csv(f"{B}/array_types_{sp}.tsv", sep="\t"); C = pd.read_csv(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_satelliteDNA/fClaHyb_{sp}.satDNA_v46.catalogue.tsv", sep="\t")
    fa = Fasta(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.fa"); runs = {a: g for a, g in C.groupby("array_id")}
    motif, purity, uid, lab = [], [], [], []
    for r in T.itertuples():
        if r.type_R8 != "micro / mini-satellite": motif.append(""); purity.append(np.nan); uid.append(np.nan); lab.append(r.type_R8); continue
        seqs = [fa.get(x.chrom, int(x.start0), int(x.end)).decode().upper() for x in runs[r.array_id].itertuples()]   # genome is soft-masked
        s = "N".join(seqs); m, p = micro_profile(s); u = int(r.unit_bp)
        a = np.frombuffer(s.encode(), dtype=np.uint8); ui = float(np.mean(a[:-u] == a[u:])) if len(a) > u else np.nan
        motif.append(m); purity.append(round(p, 3)); uid.append(round(ui, 3))
        if p >= 0.5: lab.append(f"micro ({m})n")
        elif ui >= 0.7: lab.append("mini 7-15 bp" if u <= 15 else "mini 16-30 bp" if u <= 30 else "mini 31-60 bp")
        else: lab.append("mixed / not a clean tandem")
    T["best_motif"], T["motif_purity"], T["unit_identity"], T["rlvl2"] = motif, purity, uid, lab
    top = [l for l, _ in collections.Counter(x for x in lab if x.startswith("micro (")).most_common(8)]
    T["rlvl2_plot"] = [x if not x.startswith("micro (") or x in top else "micro, other motif" for x in lab]
    T.to_csv(f"{B}/array_types_rlvl2_{sp}.tsv", sep="\t", index=False)
    mm = T[T.type_R8 == "micro / mini-satellite"]
    print(f"== {sp}: {len(mm)} micro/mini intervals (level 1) ->"); print(mm.rlvl2.value_counts().head(20).to_string())
    print("   purity of the best motif (micro/mini intervals): median", mm.motif_purity.median(), "| share >= 0.9:", round(float((mm.motif_purity >= 0.9).mean()), 2),
          "| share < 0.5:", round(float((mm.motif_purity < 0.5).mean()), 2))
    pal = plt.get_cmap("tab10"); COL = dict(L1)
    for i, l in enumerate(top): COL[l] = pal(i % 10)
    COL.update({"micro, other motif": "#d9d9a0", "mini 7-15 bp": "#9edae5", "mini 16-30 bp": "#17becf", "mini 31-60 bp": "#0b6e78", "mixed / not a clean tandem": "#000000"})
    Lc = {l.split()[0]: int(l.split()[1]) for l in open(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.fa.fai")}
    chroms = sorted(Lc, key=lambda c: int(c.split("LG")[-1])); n = len(chroms)
    fig, ax = plt.subplots(figsize=(18, 0.55 * n + 2.5)); fig.patch.set_facecolor("white")
    for i, c in enumerate(chroms):
        y = n - 1 - i; ax.add_patch(Rectangle((0, y - 0.1), Lc[c] / 1e6, 0.2, color="#eceae4", ec="#9a9890", lw=0.6))
        for r in T[T.chrom == c].itertuples():
            ax.add_patch(Rectangle((r.lo / 1e6, y - 0.36), max((r.hi - r.lo) / 1e6, 0.07), 0.24, color=COL.get(r.rlvl2_plot, "#555"), lw=0))
            if r.centromeric: ax.plot([(r.lo + r.hi) / 2e6], [y + 0.22], marker="v", ms=5, color="black")
        ax.text(-0.6, y, c.split("_")[-1], ha="right", va="center", fontsize=9)
    ax.set_xlim(-2.5, max(Lc.values()) / 1e6 + 1); ax.set_ylim(-0.7, n - 0.3); ax.set_yticks([]); ax.set_xlabel("Mb")
    for s_ in ("left", "right", "top"): ax.spines[s_].set_visible(False)
    vc = T.rlvl2_plot.value_counts(); order = list(L1) + top + ["micro, other motif", "mini 7-15 bp", "mini 16-30 bp", "mini 31-60 bp", "mixed / not a clean tandem"]
    ax.legend(handles=[Patch(color=COL[k], label=f"{k} ({vc.get(k, 0)})") for k in order if vc.get(k, 0)] + [plt.Line2D([], [], marker="v", ls="", color="black", label="centromeric (L1)")],
              loc="lower right", fontsize=8.5, frameon=False, ncol=2)
    ax.set_title(f"fClaHyb {sp}: array types, resolution level 2 - micro / mini-satellites resolved by motif and purity (n = {len(T)} multi-run intervals)", loc="left", fontsize=12, fontweight="bold")
    fig.tight_layout(); fig.savefig(f"{B}/ideogram_types_{sp}_rlvl2.png", dpi=450); plt.close(fig)
# ---------------- resolution level 3: every interval by its own identity; composite arrays by spacer type and unit period; colours = hue of the
# level-1 type, shades for the identities inside it (identities < 5 kb pooled per type); micro transparency = motif purity
RTE = pd.read_csv("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/split/rte_genome_hits_both.tsv", sep="	"); RTE["a"] = RTE[["sstart", "send"]].min(1); RTE["b"] = RTE[["sstart", "send"]].max(1)
HUE = {"composite": "Reds", "colinear": "Blues", "micro": "Greens", "mini": "PuBu", "mixed": "Greys", "gene array": "Purples", "interrupted": "Oranges", "scattered": "RdPu"}
for sp in ("Mac", "Gar"):
    T = pd.read_csv(f"{B}/array_types_rlvl2_{sp}.tsv", sep="	"); C = pd.read_csv(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_satelliteDNA/fClaHyb_{sp}.satDNA_v46.catalogue.tsv", sep="	")
    runs = {a: g.sort_values("start0") for a, g in C.groupby("array_id")}; extra = []
    for r in T.itertuples():
        g = runs[r.array_id]; gaps = list(zip(g.end.values[:-1], g.start0.values[1:])); h = RTE[RTE.chrom == r.chrom]
        w = sum(((h.a < b) & (h.b > a) & ((h.b - h.a) >= 500)).any() for a, b in gaps) / max(1, len(gaps))
        period = float(np.median(np.diff(g.start0.values))) if len(g) > 2 else np.nan
        extra.append((round(w, 2), period))
    T["rte_share"], T["period_bp"] = [e[0] for e in extra], [e[1] for e in extra]
    def l3(r):
        if r.rlvl2.startswith("micro ("): return "micro: " + r.rlvl2.split(" ", 1)[1]
        if r.rlvl2.startswith("mini "): return f"mini: {int(r.unit_bp)} bp ({r.family})"
        if r.rlvl2 == "mixed / not a clean tandem": return f"mixed: {r.family} (unit {int(r.unit_bp)})"
        if r.rlvl2 == "composite periodic":
            spc = "RTE spacer" if r.rte_share >= 0.5 else ("no-RTE spacer" if r.rte_share == 0 else "mixed spacer")
            return f"composite: {r.family}, {spc}, unit ~{r.period_bp / 1e3:.0f} kb"
        return f"{r.rlvl2.split(',')[0]}: {r.family}"
    T["rlvl3"] = [l3(r) for r in T.itertuples()]
    T["rlvl3_group"] = T.rlvl3.str.split(":").str[0]
    bpw = T.groupby("rlvl3").sat_kb.sum(); T["rlvl3_plot"] = [l if bpw[l] >= 5 else f"{grp}: {int((T[(T.rlvl3_group == grp)].rlvl3.map(bpw) < 5).groupby(T.rlvl3).any().sum())} small others"
                                                                for l, grp in zip(T.rlvl3, T.rlvl3_group)]
    T.to_csv(f"{B}/array_types_rlvl2_{sp}.tsv", sep="	", index=False)
    COL = {}; order = []
    for grp, cm in HUE.items():
        ids = T[T.rlvl3_group == grp].groupby("rlvl3_plot").sat_kb.sum().sort_values(ascending=False)
        for j, l in enumerate(ids.index): COL[l] = plt.get_cmap(cm)(0.9 - 0.6 * j / max(1, len(ids) - 1) if len(ids) > 1 else 0.75); order.append((grp, l, ids[l], (T.rlvl3_plot == l).sum()))
    Lc = {l.split()[0]: int(l.split()[1]) for l in open(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.fa.fai")}
    chroms = sorted(Lc, key=lambda c: int(c.split("LG")[-1])); n = len(chroms)
    fig = plt.figure(figsize=(24, 0.55 * n + 2)); fig.patch.set_facecolor("white"); ax = fig.add_axes([0.03, 0.05, 0.6, 0.88])
    for i, c in enumerate(chroms):
        y = n - 1 - i; ax.add_patch(Rectangle((0, y - 0.1), Lc[c] / 1e6, 0.2, color="#eceae4", ec="#9a9890", lw=0.6))
        for r in T[T.chrom == c].itertuples():
            a = 0.35 + 0.65 * (r.motif_purity if r.rlvl3_group == "micro" and not np.isnan(r.motif_purity) else 1.0)
            ax.add_patch(Rectangle((r.lo / 1e6, y - 0.36), max((r.hi - r.lo) / 1e6, 0.07), 0.24, color=COL[r.rlvl3_plot], alpha=min(1, a), lw=0))
            if r.centromeric: ax.plot([(r.lo + r.hi) / 2e6], [y + 0.22], marker="v", ms=5, color="black")
        ax.text(-0.6, y, c.split("_")[-1], ha="right", va="center", fontsize=9)
    ax.set_xlim(-2.5, max(Lc.values()) / 1e6 + 1); ax.set_ylim(-0.7, n - 0.3); ax.set_yticks([]); ax.set_xlabel("Mb")
    for s_ in ("left", "right", "top"): ax.spines[s_].set_visible(False)
    ax.set_title(f"fClaHyb {sp}: array types, resolution level 3 - every identity coloured (hue = level-1 type); composite = satellite blocks + spacers, labelled by spacer type and unit period",
                 loc="left", fontsize=12, fontweight="bold")
    ax2 = fig.add_axes([0.64, 0.02, 0.35, 0.93]); ax2.axis("off"); y0 = 1.0; step = min(0.022, 0.98 / max(1, len(order) + len(HUE)))
    for grp in HUE:
        items = [o for o in order if o[0] == grp]
        if not items: continue
        ax2.text(0, y0, grp.upper() + (" (satellite 30-80 % of the span, >= 5 interruptions / 100 kb)" if grp == "composite" else ""), fontsize=8, fontweight="bold", transform=ax2.transAxes); y0 -= step
        for _, l, kb, cnt in items:
            ax2.add_patch(Rectangle((0.0, y0 - 0.004), 0.03, step * 0.7, color=COL[l], transform=ax2.transAxes)); ax2.text(0.04, y0, f"{l.split(': ', 1)[-1]}  ({cnt}; {kb:.0f} kb)", fontsize=7, transform=ax2.transAxes, va="bottom"); y0 -= step
    fig.savefig(f"{B}/ideogram_types_{sp}_rlvl3.png", dpi=450); plt.close(fig)
    comp = T[T.rlvl3_group == "composite"]
    print(f"== {sp}: rlvl3 identities {T.rlvl3.nunique()} (plotted {T.rlvl3_plot.nunique()}); composite periodic:")
    print(comp.groupby("rlvl3").agg(n=("array_id", "size"), kb=("sat_kb", "sum"), cen=("centromeric", "sum"), rel_span=("rel_span", "median")).sort_values("kb", ascending=False).round(2).to_string())
