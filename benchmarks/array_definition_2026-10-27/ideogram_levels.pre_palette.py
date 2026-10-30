"""Array types in three nested resolution levels, both species (2026-10-27; replaces ideogram_types.py / ideogram_rlvl2.py levels).
One classification table, each level refines the previous one; all from the sequence of the runs (upper-cased: the genome is soft-masked).
Measured per multi-run step-02 interval: unit u (catalogue mode unit_bp), unit identity UI = share of i with s[i] == s[i + u]; best 1-6 bp motif
(canonical: primitive, min rotation incl. reverse complement) and its purity (union of its tandem runs >= 12 bp / bp); relative span = run bp /
span; blocks = runs; period = median distance between consecutive run starts and its coefficient of variation (CV); RTE share = interruptions
containing an RTE HSP >= 500 bp; gene flag (name_v46 gene_array:...).
rlvl1  WHAT REPEAT: gene array | microsatellite (purity >= 0.5) | minisatellite (u 7-99, k-mer recurrence >= T) | satellite (u >= 100, recurrence >= T) |
       unresolved; recurrence = share of 15-mers occurring >= 2x in the runs (indel-robust; a fixed-lag unit identity collapsed on real satellites).
       T swept 0.3 / 0.5 / 0.7; 0.5 drawn
rlvl2  HOW ARRANGED (within rlvl1): satellite / minisatellite -> colinear (relative span >= 0.8) | composite periodic (>= 5 blocks, relative span
       0.3-0.8, period CV <= 0.5) | interrupted (other relative span 0.3-0.8) | scattered (< 0.3); microsatellite -> its motif; gene array -> gene
rlvl3  WHICH ONE: family (satellite / mini), exact motif (micro), gene family; composite arrays + spacer type (RTE / no-RTE / mixed) and unit period
Outputs: array_levels_<sp>.tsv, ideogram_<sp>_rlvl1.png / _rlvl2.png / _rlvl3.png (legend in its own column, never over the tracks).
usage: python -P ideogram_levels.py"""
import os, re, sys, collections, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
B = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, B); R0 = "F:/data/fClaHyb_func_annotation_v1.0"
from kmer_homogeneity import Fasta, kcodes
BS = chr(92); RCT = str.maketrans("ACGT", "TGCA")
RTE = pd.read_csv("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/split/rte_genome_hits_both.tsv", sep="\t"); RTE["a"] = RTE[["sstart", "send"]].min(1); RTE["b"] = RTE[["sstart", "send"]].max(1)
def canon(m):
    for p in range(1, len(m) + 1):
        if len(m) % p == 0 and m[:p] * (len(m) // p) == m: m = m[:p]; break
    rc = m.translate(RCT)[::-1]; return min([m[i:] + m[:i] for i in range(len(m))] + [rc[i:] + rc[:i] for i in range(len(rc))])
def micro(s):
    cov = {}
    for k in range(1, 7):
        for x in re.finditer(("(([ACGT]{%d})" + BS + "2{2,})") % k, s):
            if x.end() - x.start() >= 12:
                m = canon(x.group(2)); v = cov.setdefault(m, np.zeros(len(s), dtype=bool)); v[x.start():x.end()] = True
    if not cov: return "", 0.0
    m, v = max(cov.items(), key=lambda t: t[1].sum()); return m, float(v.sum()) / max(1, len(s))
def unit_identity(s, u):
    """median over windows of 3 units (>= 300 bp) of the best identity at lags within +-10 % of the unit (indels and monomer-length variation
    shift a fixed lag; the step-02 scan used the same best-lag-near-the-period idea)"""
    a = np.frombuffer(s.encode(), dtype=np.uint8); W = max(300, 3 * u); lags = range(max(1, int(0.9 * u)), int(1.1 * u) + 2); vals = []
    for w0 in range(0, len(a) - W - max(lags), W):
        x = a[w0:w0 + W]
        if (x == 78).any(): continue
        vals.append(max(float(np.mean(x == a[w0 + l:w0 + l + W])) for l in lags))
    return float(np.median(vals)) if vals else np.nan
def gene_kind(name): return name.split("|")[0].split(":")[1] if name.startswith("gene_array:") else ""
def level1(r, T):
    if r.gene: return "gene array"
    if r.motif_purity >= 0.5: return "microsatellite"
    if 7 <= r.unit_bp < 100 and r.kmer_recurrence >= T: return "minisatellite"
    if r.unit_bp >= 100 and r.kmer_recurrence >= T: return "satellite"
    return "unresolved"
def arrangement(r):
    if r.rel_span >= 0.8: return "colinear"
    if r.rel_span < 0.3: return "scattered"
    if r.runs >= 5 and r.period_cv <= 0.5: return "composite periodic"
    return "interrupted"
def lv2(r):
    if r.rlvl1 == "microsatellite": return f"({r.best_motif})n"
    if r.rlvl1 == "gene array": return gene_kind(r.name)
    if r.rlvl1 == "unresolved": return "unresolved"
    return arrangement(r)
def lv3(r):
    if r.rlvl1 == "microsatellite": return f"({r.best_motif})n"
    if r.rlvl1 == "gene array": return f"{gene_kind(r.name)} ({r.family})"
    if r.rlvl1 == "unresolved": return f"unresolved {r.family} (unit {r.unit_bp})"
    if r.rlvl2 == "composite periodic":
        spc = "RTE spacer" if r.rte_share >= 0.5 else ("no-RTE spacer" if r.rte_share == 0 else "mixed spacer")
        return f"{r.family} composite, {spc}, unit ~{r.period / 1e3:.0f} kb"
    return f"{r.family} {r.rlvl2}"
L1COL = {"satellite": "#c0392b", "minisatellite": "#2e86c1", "microsatellite": "#27ae60", "gene array": "#7f7f7f", "unresolved": "#d4ac0d"}
L1HUE = {"satellite": "Reds", "minisatellite": "Blues", "microsatellite": "Greens", "gene array": "Greys", "unresolved": "YlOrBr"}
def draw(T, sp, col_of, level, legend_rows, title):
    Lc = {l.split()[0]: int(l.split()[1]) for l in open(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.fa.fai")}
    chroms = sorted(Lc, key=lambda c: int(c.split("LG")[-1])); n = len(chroms); H = max(0.55 * n + 2, 0.2 * len(legend_rows) + 2)
    fig = plt.figure(figsize=(24, H)); fig.patch.set_facecolor("white"); ax = fig.add_axes([0.04, 0.4 / H, 0.6, 1 - 1.4 / H])
    for i, c in enumerate(chroms):
        y = n - 1 - i; ax.add_patch(Rectangle((0, y - 0.1), Lc[c] / 1e6, 0.2, color="#eceae4", ec="#9a9890", lw=0.6))
        for r in T[T.chrom == c].itertuples():
            ax.add_patch(Rectangle((r.lo / 1e6, y - 0.36), max((r.hi - r.lo) / 1e6, 0.07), 0.24, color=col_of(r), lw=0))
            if r.centromeric: ax.plot([(r.lo + r.hi) / 2e6], [y + 0.24], marker="v", ms=5, color="black")
        ax.text(-0.6, y, c.split("_")[-1], ha="right", va="center", fontsize=9)
    ax.set_xlim(-2.5, max(Lc.values()) / 1e6 + 1); ax.set_ylim(-0.7, n - 0.3); ax.set_yticks([]); ax.set_xlabel("Mb")
    for s_ in ("left", "right", "top"): ax.spines[s_].set_visible(False)
    ax.set_title(title, loc="left", fontsize=12, fontweight="bold")
    lg = fig.add_axes([0.66, 0.4 / H, 0.33, 1 - 1.4 / H]); lg.axis("off"); lg.set_xlim(0, 1); lg.set_ylim(0, 1)
    step = min(0.03, 0.98 / max(1, len(legend_rows))); y = 0.99
    for kind, text, color in legend_rows:
        if kind == "head": lg.text(0, y, text, fontsize=8.5, fontweight="bold", va="top")
        else: lg.add_patch(Rectangle((0.0, y - step * 0.75), 0.035, step * 0.6, color=color)); lg.text(0.05, y, text, fontsize=7.5, va="top")
        y -= step
    lg.plot([0.01], [y - step], marker="v", color="black", ms=5); lg.text(0.05, y - step + step * 0.4, "centromeric (L1)", fontsize=7.5, va="top")
    fig.savefig(f"{B}/ideogram_{sp}_rlvl{level}.png", dpi=110); plt.close(fig)
for sp in ("Mac", "Gar"):
    C = pd.read_csv(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_satelliteDNA/fClaHyb_{sp}.satDNA_v46.catalogue.tsv", sep="\t"); fa = Fasta(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.fa")
    rows = []
    for aid, g in C.groupby("array_id"):
        if len(g) < 2: continue
        g = g.sort_values("start0"); s = "N".join(fa.get(x.chrom, int(x.start0), int(x.end)).decode().upper() for x in g.itertuples())
        u = int(g.unit_bp.mode().iloc[0]); ui = unit_identity(s, u)
        cnt = np.unique(kcodes(s.encode(), 15), return_counts=True)[1]; rec = float(cnt[cnt >= 2].sum() / max(1, cnt.sum()))   # share of 15-mers occurring >= 2x
        m, p = micro(s); lo, hi = int(g.start0.min()), int(g.end.max()); d = np.diff(g.start0.values)
        gaps = list(zip(g.end.values[:-1], g.start0.values[1:])); h = RTE[RTE.chrom == g.chrom.iloc[0]]
        rte = sum(((h.a < b) & (h.b > a_) & ((h.b - h.a) >= 500)).any() for a_, b in gaps) / max(1, len(gaps))
        rows.append(dict(array_id=aid, chrom=g.chrom.iloc[0], lo=lo, hi=hi, runs=len(g), sat_kb=round(g.bp.sum() / 1e3, 2), rel_span=round(g.bp.sum() / max(1, hi - lo), 3),
                         unit_bp=u, unit_identity=round(ui, 3), kmer_recurrence=round(rec, 3), best_motif=m, motif_purity=round(p, 3), period=float(np.median(d)) if len(d) else np.nan,
                         period_cv=round(float(np.std(d) / np.mean(d)), 2) if len(d) > 1 and np.mean(d) > 0 else np.nan, rte_share=round(rte, 2),
                         gene=str(g.name_v46.iloc[0]).startswith("gene_array"), name=str(g.name_v46.iloc[0]), family=g.family_v46.iloc[0],
                         centromeric=(g.putative_centromeric_array_L1 == "yes").any()))
    T = pd.DataFrame(rows)
    for t in (0.3, 0.5, 0.7): T[f"rlvl1_T{int(t * 10)}"] = [level1(r, t) for r in T.itertuples()]
    T["rlvl1"] = T.rlvl1_T5; T["rlvl2"] = [lv2(r) for r in T.itertuples()]; T["rlvl3"] = [lv3(r) for r in T.itertuples()]
    T.to_csv(f"{B}/array_levels_{sp}.tsv", sep="\t", index=False)
    print(f"== {sp}: {len(T)} multi-run intervals; rlvl1 by k-mer recurrence cut T:")
    print(pd.DataFrame({f"T = {t}": T[f"rlvl1_T{int(t * 10)}"].value_counts() for t in (0.3, 0.5, 0.7)}).fillna(0).astype(int).to_string())
    print("   rlvl2 within satellite:", T[T.rlvl1 == "satellite"].rlvl2.value_counts().to_dict(), "| within minisatellite:", T[T.rlvl1 == "minisatellite"].rlvl2.value_counts().to_dict())
    print("   centromeric:", T[T.centromeric].groupby(["rlvl1", "rlvl2"]).size().to_dict())
    kb = lambda col, v: T[T[col] == v].sat_kb.sum()
    # level 1
    rows1 = [("item", f"{k}  ({(T.rlvl1 == k).sum()} intervals; {kb('rlvl1', k):.0f} kb)", c) for k, c in L1COL.items() if (T.rlvl1 == k).any()]
    draw(T, sp, lambda r: L1COL[r.rlvl1], 1, [("head", "rlvl1: what repeat (from the sequence)", None)] + rows1,
         f"fClaHyb {sp} - resolution level 1: what repeat (micro: 1-6 bp motif >= 50 %; mini: unit 7-99 bp; satellite: unit >= 100 bp; tandem = >= 50 % of 15-mers recur)")
    # level 2 and 3: hue of rlvl1, shades for the sub-levels; < 5 kb identities pooled per rlvl1 class at level 3
    for level, col in ((2, "rlvl2"), (3, "rlvl3")):
        T["lab"] = T[col]
        if level == 3:
            w = T.groupby("rlvl3").sat_kb.sum(); small = T.rlvl3.map(w) < 5
            T.loc[small, "lab"] = T.loc[small, "rlvl1"] + ": " + T.loc[small].groupby("rlvl1").rlvl3.transform("nunique").astype(str) + " small identities (< 5 kb each)"
        COL, leg = {}, []
        for c1, cm in L1HUE.items():
            ids = T[T.rlvl1 == c1].groupby("lab").sat_kb.sum().sort_values(ascending=False)
            if not len(ids): continue
            leg.append(("head", c1.upper(), None))
            for j, l in enumerate(ids.index):
                COL[(c1, l)] = plt.get_cmap(cm)(0.9 - 0.55 * j / max(1, len(ids) - 1)); leg.append(("item", f"{l}  ({(T[T.rlvl1 == c1]['lab'] == l).sum()}; {ids[l]:.0f} kb)", COL[(c1, l)]))
        draw(T, sp, lambda r: COL[(r.rlvl1, r.lab)], level, leg,
             f"fClaHyb {sp} - resolution level {level}: " + ("how each repeat is arranged (colinear / composite periodic / interrupted / scattered; motif; gene)" if level == 2 else
                                                             "which one exactly (family / motif / gene; composite + spacer type and unit period)"))
    print("   wrote rlvl1-3", sp)