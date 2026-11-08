"""Divergence, span and occupancy between homologous compartments of anchored units (2026-10-23).
User rules: fixed boundaries; compare only homologous segments; no pooling of monomers; span (aligned length) and occupancy (which reference
coordinates are present) are complementary metrics, reported separately from divergence.
Units: unit_anchored/units_A6500.tsv (anchor = element position 6,500; oriented so the element reads forward).
Compartments per unit:
  element windows: element consensus coordinates in 500-bp windows; the unit's bases mapped to window w through the element hits (strand-aware);
  satellite block: from the first to the last MacSat03 monomer inside the unit (36_monomers_final), oriented with the unit; edge monomers = first and last
    complete monomer of the block (next to fixed element landmarks);
  insert windows: insert coordinates in 500-bp windows (type II units only).
Metrics per homologous pair:
  divergence = K2P over aligned, gap-free columns of a global alignment (free end gaps); only when both segments occupy >= 90% of a window;
  span = aligned columns / length of the longer segment, and absolute length difference (bp; satellite block also in monomer equivalents of 224 bp);
  occupancy = covered consensus bp / window bp per unit; shared occupancy = Jaccard of the covered consensus coordinates of the two units.
Pair categories: same array at lag 1-5; different centromeres, same type (<= 20 pairs per centromere pair); type I vs type II (<= 20 per pair).
Summaries: median per centromere pair, then median across pairs; 95% CI by hierarchical bootstrap (centromeres, then units), 200 replicates.
Outputs: homologous/pairs.tsv, homologous/summary.tsv, figure/homologous_divergence.png/.pdf"""
import pathlib, random, itertools, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/code"); import ltr_colours as LC
from Bio import Align
import sys; sys.path.insert(0, str(pathlib.Path(__file__).parent)); from satblock import block_and_edges, window_hits

W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22"); A = pathlib.Path("F:/data/fClaHyb_func_annotation_v1.0")
O = pathlib.Path("F:/data/satellite_FISH/pipeline_v4_5_fixed/out/fClaHyb_Mac_satellites"); OUT = W / "homologous"; OUT.mkdir(exist_ok=True)
random.seed(20261023); WS = 500; EL = 8164; INS = 6068; MON = 224
T2 = {"C_mac_LG08", "C_mac_LG25", "C_mac_LG16"}
lab = lambda c: "LG" + c.split("LG")[-1]
class Fa:
    def __init__(self, p):
        self.f = open(p, "rb"); self.i = {}
        for l in open(str(p) + ".fai"):
            n, ln, off, lb, lby = l.split("\t")[:5]; self.i[n] = (int(ln), int(off), int(lb), int(lby))
    def get(self, c, a, b):
        ln, off, lb, lby = self.i[c]; a = max(0, int(a)); b = min(ln, int(b))
        s0 = off + (a // lb) * lby + a % lb; s1 = off + (b // lb) * lby + b % lb
        self.f.seek(s0); return self.f.read(s1 - s0).replace(b"\n", b"").replace(b"\r", b"").decode().upper()
fa = Fa(A / "fClaHyb_Mac/fClaHyb_Mac_genome/fClaHyb_Mac_LG.fa"); comp = str.maketrans("ACGTN", "TGCAN")
rc = lambda s: s.translate(comp)[::-1]
AL = Align.PairwiseAligner(mode="global", match_score=2, mismatch_score=-3, open_gap_score=-5, extend_gap_score=-2)
AL.target_end_gap_score = AL.query_end_gap_score = 0
TS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

def k2p_span(a, b):
    if not a or not b: return np.nan, np.nan, np.nan
    aln = AL.align(a, b)[0]; P = Q = n = 0
    for (s0, s1), (t0, t1) in zip(*aln.aligned):
        for x, y in zip(a[s0:s1], b[t0:t1]):
            if x in "ACGT" and y in "ACGT":
                n += 1
                if x != y: P += (x, y) in TS; Q += (x, y) not in TS
    if n < 50: return np.nan, n / max(len(a), len(b)), n
    p, q = P / n, Q / n
    try: d = -0.5 * np.log(1 - 2 * p - q) - 0.25 * np.log(1 - 2 * q)
    except Exception: d = np.nan
    return d, n / max(len(a), len(b)), n

U = pd.read_csv(W / "unit_anchored" / "units_A6500.tsv", sep="\t").reset_index(drop=True)
U["type"] = np.where(U.chrom.isin(T2), "II", "I")
def hits(path):
    h = pd.read_csv(path, sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
    h["lo"] = h[["s0", "s1"]].min(1); h["hi"] = h[["s0", "s1"]].max(1); return h
eh = hits(W / "copies" / "Mac.hits.tsv"); ih = hits(W / "insert" / "Mac.insert_hits.tsv")
mono = pd.read_csv(O / "36_monomers_final.tsv", sep="\t", usecols=["chrom", "genomic_start_0based", "genomic_end", "family_id", "complete_or_truncated"])
mono = mono[mono.family_id == "MacSat03"]

def windows(h, L, unit_strand):
    """per consensus window: covered consensus coordinates (set of 10-bp bins) and the unit sequence mapped to it, in element orientation"""
    out = {}
    for w0 in range(0, L, WS):
        w1 = min(L, w0 + WS); seg = []; cov = set()
        for r in window_hits(h, w0, w1, sat_filter=(L == EL)):
            a, b = max(r.q0, w0), min(r.q1, w1)
            if b - a < 20: continue
            f = lambda q: r.s0 + (q - r.q0) * (r.s1 - r.s0) / max(1, r.q1 - r.q0)
            g0, g1 = sorted((int(f(a)), int(f(b)))); s = fa.get(r.chrom, g0, g1)
            seg.append((a, s if r.s1 >= r.s0 else rc(s))); cov |= set(range(a // 10, b // 10))
        seq = "".join(s for _, s in sorted(seg))
        out[w0] = dict(seq=seq, cov=cov, occ=len(cov) * 10 / (w1 - w0))
    return out

feat = []
for r in U.itertuples():
    lo, hi = int(r.lo), int(r.hi)
    he = eh[(eh.chrom == r.chrom) & (eh.hi > lo) & (eh.lo < hi)]
    hn = ih[(ih.chrom == r.chrom) & (ih.hi > lo) & (ih.lo < hi) & (ih.len >= 300)]
    m = mono[(mono.chrom == r.chrom) & (mono.genomic_start_0based >= lo) & (mono.genomic_end <= hi)]
    block, first_, last_, nmono = block_and_edges(m, he, r.strand, fa.get, r.chrom, rc)
    feat.append(dict(el=windows(he, EL, r.strand), ins=windows(hn, INS, r.strand) if r.type == "II" else {}, block=block, n_mono=nmono,
                     first=first_, last=last_))
print("units featured:", len(feat), flush=True)

# pairs
pairs = []
for c, g in U.groupby("chrom"):
    idx = g.sort_values("lo").index.tolist()
    for lag in range(1, 6):
        for i in range(len(idx) - lag):
            if U.loc[idx[i + lag], "lo"] - U.loc[idx[i], "hi"] <= 2000 * lag:             # truly consecutive along the array
                pairs.append((idx[i], idx[i + lag], f"same array, lag {lag}"))
cs = sorted(U.chrom.unique())
for c1, c2 in itertools.combinations(cs, 2):
    i1, i2 = U.index[U.chrom == c1].tolist(), U.index[U.chrom == c2].tolist()
    cat = "different centromeres, same type" if U.loc[i1[0], "type"] == U.loc[i2[0], "type"] else "type I vs type II"
    for _ in range(20): pairs.append((random.choice(i1), random.choice(i2), cat))
print("pairs:", len(pairs), flush=True)

rows = []
for n_, (i, j, cat) in enumerate(pairs):
    fi, fj = feat[i], feat[j]; base = dict(category=cat, c1=lab(U.loc[i, "chrom"]), c2=lab(U.loc[j, "chrom"]), u1=i, u2=j)
    for w0 in fi["el"]:
        a, b = fi["el"][w0], fj["el"][w0]
        shared = len(a["cov"] & b["cov"]) / max(1, len(a["cov"] | b["cov"]))
        row = dict(base, compartment="element", window=w0, occ1=a["occ"], occ2=b["occ"], shared_occupancy=shared)
        if a["occ"] >= 0.9 and b["occ"] >= 0.9:
            d, sp, n = k2p_span(a["seq"], b["seq"]); row.update(divergence=d, span=sp, aligned=n, len_diff=abs(len(a["seq"]) - len(b["seq"])))
        rows.append(row)
    if fi["ins"] and fj["ins"]:
        for w0 in fi["ins"]:
            a, b = fi["ins"][w0], fj["ins"][w0]
            row = dict(base, compartment="insert", window=w0, occ1=a["occ"], occ2=b["occ"], shared_occupancy=len(a["cov"] & b["cov"]) / max(1, len(a["cov"] | b["cov"])))
            if a["occ"] >= 0.9 and b["occ"] >= 0.9:
                d, sp, n = k2p_span(a["seq"], b["seq"]); row.update(divergence=d, span=sp, aligned=n, len_diff=abs(len(a["seq"]) - len(b["seq"])))
            rows.append(row)
    for part in ("first", "last"):
        d, sp, n = k2p_span(fi[part], fj[part])
        rows.append(dict(base, compartment=f"satellite edge, {part} 448 bp", window=-1, divergence=d, span=sp, aligned=n, len_diff=abs(len(fi[part]) - len(fj[part]))))
    if fi["block"] and fj["block"] and max(len(fi["block"]), len(fj["block"])) <= 12000:
        d, sp, n = k2p_span(fi["block"], fj["block"])
        rows.append(dict(base, compartment="satellite block", window=-1, divergence=d, span=sp, aligned=n, len_diff=abs(len(fi["block"]) - len(fj["block"])),
                         monomer_diff=abs(fi["n_mono"] - fj["n_mono"]), indel_monomer_equiv=abs(len(fi["block"]) - len(fj["block"])) / MON))
    if n_ % 200 == 0: print("pair", n_, flush=True)
P = pd.DataFrame(rows); P.to_csv(OUT / "pairs.tsv", sep="\t", index=False)

# balanced summaries + hierarchical bootstrap
def balanced(df, col):
    v = df.dropna(subset=[col]).groupby(["c1", "c2"])[col].median()
    return v.median() if len(v) else np.nan
CATS = ["same array, lag 1", "same array, lag 2", "same array, lag 5", "different centromeres, same type", "type I vs type II"]
summ = []
for cat in CATS:
    d0 = P[P.category == cat]
    for comp_ in ["element", "satellite edge, first 448 bp", "satellite edge, last 448 bp", "satellite block", "insert"]:
        d = d0[d0.compartment == comp_]
        if d.empty: continue
        rec = dict(category=cat, compartment=comp_, pairs=len(d))
        for col in ("divergence", "span", "shared_occupancy", "indel_monomer_equiv"):
            if col not in d or d[col].notna().sum() == 0: continue
            est = balanced(d, col); boots = []
            cents = sorted(set(d.c1) | set(d.c2))
            for _ in range(200):
                pick = random.choices(cents, k=len(cents)); sub = d[d.c1.isin(pick) & d.c2.isin(pick)]
                if len(sub): sub = sub.sample(frac=1, replace=True, random_state=random.randint(0, 10 ** 6)); boots.append(balanced(sub, col))
            rec[col] = est; rec[col + "_lo"] = np.nanpercentile(boots, 2.5) if boots else np.nan; rec[col + "_hi"] = np.nanpercentile(boots, 97.5) if boots else np.nan
        summ.append(rec)
S = pd.DataFrame(summ); S.round(4).to_csv(OUT / "summary.tsv", sep="\t", index=False)
print(S[["category", "compartment", "pairs", "divergence", "span", "shared_occupancy", "indel_monomer_equiv"]].round(3).to_string(index=False))

# figure: three metrics x compartments x categories
fig, axs = plt.subplots(1, 3, figsize=(18, 6), facecolor="white")
COMPS = ["element", "satellite edge, first 448 bp", "satellite edge, last 448 bp", "satellite block", "insert"]
COL = dict(zip(CATS, [LC.REL[0], LC.REL[1], LC.REL[2], LC.REL[3], LC.ELEM]))
for ax, (m, t) in zip(axs, (("divergence", "Divergence (K2P, aligned bases)"), ("span", "Span (aligned / longer length)"), ("shared_occupancy", "Shared occupancy (Jaccard of covered coordinates)"))):
    x = np.arange(len(COMPS)); w = 0.16
    for k, cat in enumerate(CATS):
        s = S[S.category == cat].set_index("compartment").reindex(COMPS)
        if m not in s: continue
        y = s[m].values; lo = s.get(m + "_lo", s[m]).values; hi = s.get(m + "_hi", s[m]).values
        ax.bar(x + (k - 2) * w, y, width=w, color=COL[cat], label=cat)
        ax.errorbar(x + (k - 2) * w, y, yerr=[np.nan_to_num(y - lo), np.nan_to_num(hi - y)], fmt="none", ecolor="black", lw=0.8, capsize=2)
    ax.set_xticks(x); ax.set_xticklabels(["element\nwindows", "first\nmonomer", "last\nmonomer", "satellite\nblock", "insert\nwindows"], fontsize=9)
    ax.set_title(t, loc="left", fontsize=12, fontweight="bold")
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
axs[0].legend(fontsize=8, frameon=False)
fig.suptitle("Homologous compartments of anchored units: divergence, span and occupancy (balanced by centromere pair; 95% hierarchical bootstrap)",
             fontsize=13, fontweight="bold", x=0.02, ha="left")
fig.tight_layout()
for ext in ("png", "pdf"): fig.savefig(W / "figure" / f"homologous_divergence.{ext}", dpi=450, bbox_inches="tight")
