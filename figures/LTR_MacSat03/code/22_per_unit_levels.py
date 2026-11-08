"""Per-unit divergence at four levels, alignment + k-mer grid (2026-10-23; DESIGN_BOARD.md sections 2-5).
Homology: element windows (500 bp of consensus coordinates) used only when BOTH members of a pair occupy >= 90% of the window; metrics are summed
over these shared windows (alignment: P, Q, n summed -> K2P; k-mers: intersections/unions/min/max/dot summed per window, k-mers counted per window).
Satellite: whole satellite block, k-mer metrics only (length-robust ones: containment, weighted Jaccard, cosine; Mash shown for comparison);
edge monomers (first, last) by alignment K2P.
Arrays: maximal runs of consecutive anchored units (<= 5 kb between units).
Levels per unit: L-a lag-1 neighbours in its array; L-b units of other arrays of the same chromosome (<= 10); L-c units of other chromosomes, same type
(<= 10) and other type (<= 10, shared windows only); L-d Gar element copies (<= 20; element windows only).
k-mer grid: k = 11, 15, 21, 31 x measure (Mash, containment = max of both directions, weighted-Jaccard distance, cosine distance).
Outputs: homologous/per_unit.tsv (one row per unit x level x method), homologous/per_unit_wide.tsv, homologous/arrays.tsv,
figure/per_unit_divergence.png/.pdf, figure/per_unit_kmer_grid.png/.pdf"""
import pathlib, random, collections, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/code"); import ltr_colours as LC
from Bio import Align
import sys; sys.path.insert(0, str(pathlib.Path(__file__).parent)); from satblock import block_and_edges, window_hits

W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22"); A = pathlib.Path("F:/data/fClaHyb_func_annotation_v1.0")
O = pathlib.Path("F:/data/satellite_FISH/pipeline_v4_5_fixed/out/fClaHyb_Mac_satellites"); OUT = W / "homologous"; OUT.mkdir(exist_ok=True)
random.seed(20261023); WS = 500; EL = 8164; KS = [11, 15, 21, 31]
T2 = {"C_mac_LG08", "C_mac_LG25", "C_mac_LG16"}; lab = lambda c: "LG" + c.split("LG")[-1]
class Fa:
    def __init__(self, p):
        self.f = open(p, "rb"); self.i = {}
        for l in open(str(p) + ".fai"):
            n, ln, off, lb, lby = l.split("\t")[:5]; self.i[n] = (int(ln), int(off), int(lb), int(lby))
    def get(self, c, a, b):
        ln, off, lb, lby = self.i[c]; a = max(0, int(a)); b = min(ln, int(b))
        s0 = off + (a // lb) * lby + a % lb; s1 = off + (b // lb) * lby + b % lb
        self.f.seek(s0); return self.f.read(s1 - s0).replace(b"\n", b"").replace(b"\r", b"").decode().upper()
mfa = Fa(A / "fClaHyb_Mac/fClaHyb_Mac_genome/fClaHyb_Mac_LG.fa"); gfa = Fa(A / "fClaHyb_Gar/fClaHyb_Gar_genome/fClaHyb_Gar_LG.fa")
comp = str.maketrans("ACGTN", "TGCAN"); rc = lambda s: s.translate(comp)[::-1]
AL = Align.PairwiseAligner(mode="global", match_score=2, mismatch_score=-3, open_gap_score=-5, extend_gap_score=-2)
AL.target_end_gap_score = AL.query_end_gap_score = 0
TS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
def pq(a, b):
    if not a or not b: return 0, 0, 0
    aln = AL.align(a, b)[0]; P = Q = n = 0
    for (s0, s1), (t0, t1) in zip(*aln.aligned):
        for x, y in zip(a[s0:s1], b[t0:t1]):
            if x in "ACGT" and y in "ACGT":
                n += 1
                if x != y: P += (x, y) in TS; Q += (x, y) not in TS
    return P, Q, n
def k2p(P, Q, n):
    if n < 100: return np.nan
    p, q = P / n, Q / n
    with np.errstate(all="ignore"):
        v = -0.5 * np.log(1 - 2 * p - q) - 0.25 * np.log(1 - 2 * q)
    return float(v) if np.isfinite(v) else np.nan
def kc(s, k):
    c = collections.Counter()
    for i in range(len(s) - k + 1):
        x = s[i:i + k]
        if "N" not in x: c[min(x, rc(x))] += 1
    return c

def hits(path):
    h = pd.read_csv(path, sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
    h["lo"] = h[["s0", "s1"]].min(1); h["hi"] = h[["s0", "s1"]].max(1); return h
def windows(h, fa):
    out = {}
    for w0 in range(0, EL, WS):
        w1 = min(EL, w0 + WS); seg = []; cov = 0
        for r in window_hits(h, w0, w1):
            a, b = max(r.q0, w0), min(r.q1, w1)
            if b - a < 20: continue
            f = lambda q: r.s0 + (q - r.q0) * (r.s1 - r.s0) / max(1, r.q1 - r.q0)
            g0, g1 = sorted((int(f(a)), int(f(b)))); s = fa.get(r.chrom, g0, g1)
            seg.append((a, s if r.s1 >= r.s0 else rc(s))); cov += b - a
        if cov / (w1 - w0) >= 0.9:
            seq = "".join(s for _, s in sorted(seg)); out[w0] = dict(seq=seq, k={k: kc(seq, k) for k in KS})
    return out

U = pd.read_csv(W / "unit_anchored" / "units_A6500.tsv", sep="\t").sort_values(["chrom", "lo"]).reset_index(drop=True)
U["type"] = np.where(U.chrom.isin(T2), "II", "I")
arr = []; aid = 0
for c, g in U.groupby("chrom", sort=False):
    prev = None
    for i, r in g.iterrows():
        if prev is None or r.lo - prev > 5000: aid += 1
        arr.append((i, f"{lab(c)}_A{aid}")); prev = r.hi
U["array"] = pd.Series(dict(arr))
U.groupby(["chrom", "array"]).size().rename("units").to_csv(OUT / "arrays.tsv", sep="\t")
eh = hits(W / "copies" / "Mac.hits.tsv"); gh = hits(W / "copies" / "Gar.hits.tsv")
mono = pd.read_csv(O / "36_monomers_final.tsv", sep="\t", usecols=["chrom", "genomic_start_0based", "genomic_end", "family_id", "complete_or_truncated"])
mono = mono[(mono.family_id == "MacSat03")]
F = []
for r in U.itertuples():
    he = eh[(eh.chrom == r.chrom) & (eh.hi > r.lo) & (eh.lo < r.hi)]
    m = mono[(mono.chrom == r.chrom) & (mono.genomic_start_0based >= r.lo) & (mono.genomic_end <= r.hi)]
    block, first_, last_, nmono = block_and_edges(m, he, r.strand, mfa.get, r.chrom, rc)
    F.append(dict(win=windows(he, mfa), block={k: kc(block, k) for k in KS}, first=first_, last=last_))
print("units featured", len(F), "arrays", U.array.nunique(), flush=True)
gc = pd.read_csv(W / "core" / "Gar.core_copies.tsv", sep="\t").sample(n=60, random_state=1)
G = [dict(win=windows(gh[(gh.chrom == r.chrom) & (gh.hi > r.lo - 8000) & (gh.lo < r.hi + 8000)], gfa)) for r in gc.itertuples()]
print("Gar copies featured", len(G), flush=True)

def kmer_metrics(pairs_of_counters):
    out = {}
    for k in KS:
        inter = uni = mn = mx = dot = na = nb = sa = sb = 0
        for a, b in pairs_of_counters:
            ca, cb = a[k], b[k]; ka, kb = set(ca), set(cb); I = ka & kb
            inter += len(I); uni += len(ka | kb); sa += len(ka); sb += len(kb)
            mn += sum(min(ca[x], cb[x]) for x in I); mx += sum(ca.values()) + sum(cb.values()) - sum(min(ca[x], cb[x]) for x in I)
            dot += sum(ca[x] * cb[x] for x in I); na += sum(v * v for v in ca.values()); nb += sum(v * v for v in cb.values())
        if uni == 0: continue
        j = inter / uni
        out[f"k{k}_mash"] = 1.0 if j == 0 else max(0.0, -1 / k * np.log(2 * j / (1 + j)))
        out[f"k{k}_containment"] = 1 - max(inter / max(sa, 1), inter / max(sb, 1))
        out[f"k{k}_wjaccard"] = 1 - mn / max(mx, 1)
        out[f"k{k}_cosine"] = 1 - dot / np.sqrt(max(na * nb, 1))
    return out

def compare(fa_, fb_, sat=True):
    shared = sorted(set(fa_["win"]) & set(fb_["win"]))
    P = Q = n = 0
    for w in shared:
        p, q, m = pq(fa_["win"][w]["seq"], fb_["win"][w]["seq"]); P += p; Q += q; n += m
    rec = dict(shared_windows=len(shared), el_align_k2p=k2p(P, Q, n))
    rec.update({f"el_{k}": v for k, v in kmer_metrics([(fa_["win"][w]["k"], fb_["win"][w]["k"]) for w in shared]).items()})
    if sat and fa_.get("block") and fb_.get("block") and sum(fa_["block"][11].values()) and sum(fb_["block"][11].values()):
        rec.update({f"sat_{k}": v for k, v in kmer_metrics([(fa_["block"], fb_["block"])]).items()})
        for part in ("first", "last"):
            rec[f"sat_{part}_k2p"] = k2p(*pq(fa_[part], fb_[part])) if fa_[part] and fb_[part] else np.nan
    return rec

rows = []
for i, r in U.iterrows():
    r_arr, r_chr, r_typ, r_lo = r["array"], r["chrom"], r["type"], r["lo"]
    same_arr = U[(U.array == r_arr)]; pos = list(same_arr.index); k_ = pos.index(i)
    targets = {"L-a within array (lag 1)": [p for p in (pos[k_ - 1] if k_ > 0 else None, pos[k_ + 1] if k_ + 1 < len(pos) else None) if p is not None],
               "L-b other arrays, same chromosome": random.sample(list(U.index[(U.chrom == r_chr) & (U.array != r_arr)]), min(10, ((U.chrom == r_chr) & (U.array != r_arr)).sum())),
               "L-c other chromosomes, same type": random.sample(list(U.index[(U.chrom != r_chr) & (U.type == r_typ)]), min(10, ((U.chrom != r_chr) & (U.type == r_typ)).sum())),
               "L-c other type": random.sample(list(U.index[U.type != r_typ]), min(10, (U.type != r_typ).sum()))}
    for lvl, tg in targets.items():
        if not tg: continue
        recs = [compare(F[i], F[j]) for j in tg]
        df = pd.DataFrame(recs).median(numeric_only=True).to_dict()
        rows.append(dict(unit=i, chrom=lab(r_chr), array=r_arr, type=r_typ, lo=r_lo, level=lvl, n=len(tg), **df))
    recs = [compare(F[i], g, sat=False) for g in random.sample(G, 20)]
    rows.append(dict(unit=i, chrom=lab(r_chr), array=r_arr, type=r_typ, lo=r_lo, level="L-d Gar element copies", n=20,
                     **pd.DataFrame(recs).median(numeric_only=True).to_dict()))
    if i % 25 == 0: print("unit", i, flush=True)
R = pd.DataFrame(rows); R.to_csv(OUT / "per_unit.tsv", sep="\t", index=False)
LV = ["L-a within array (lag 1)", "L-b other arrays, same chromosome", "L-c other chromosomes, same type", "L-c other type", "L-d Gar element copies"]
summ = R.groupby(["level"])[["el_align_k2p", "el_k15_mash", "el_k21_containment", "sat_k15_wjaccard", "sat_first_k2p"]].median().reindex(LV)
print(summ.round(4).to_string())
# claim stability across the k-mer grid: within array < other arrays < other chromosomes < Gar (element), per measure
order_ok = {}
for col in [c for c in R.columns if c.startswith("el_k")]:
    v = R.groupby("level")[col].median().reindex(LV)
    order_ok[col] = bool(v[LV[0]] < v[LV[2]] < v[LV[4]])
pd.Series(order_ok, name="within<between_chrom<Gar").to_csv(OUT / "kmer_grid_order.tsv", sep="\t")
print("ordering within < between chromosomes < Gar holds in", sum(order_ok.values()), "of", len(order_ok), "element k-mer variants")

# figure 1: per unit along the arrays (alignment K2P), one row per centromere
cs = [c for c in ["LG07", "LG04", "LG01", "LG23", "LG02", "LG15", "LG19", "LG21", "LG26", "LG08", "LG25", "LG16"] if c in set(R.chrom)]
fig, axs = plt.subplots(len(cs), 1, figsize=(15, 2.0 * len(cs)), sharex=False, facecolor="white")
COL = dict(zip(LV, [LC.REL[0], LC.REL[1], LC.REL[2], LC.REL[3], LC.ELEM]))
for ax, c in zip(axs, cs):
    d = R[R.chrom == c]; x0 = d.lo.min()
    for lvl in LV:
        e = d[d.level == lvl].sort_values("lo")
        if len(e): ax.plot((e.lo - x0) / 1000, e.el_align_k2p, "-o", ms=3, lw=1, color=COL[lvl], label=lvl)
    for a_ in d.array.unique():
        e = d[d.array == a_]; ax.axvspan((e.lo.min() - x0) / 1000, (e.lo.max() - x0) / 1000 + 12, color="#eceae4", zorder=0)
    ax.set_ylabel(f"{c}\n({'II' if 'C_mac_' + c in T2 else 'I'})", fontsize=10, rotation=0, ha="right", va="center"); ax.set_ylim(0, None)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axs[0].legend(fontsize=8, ncol=5, frameon=False, loc="upper left", bbox_to_anchor=(0, 1.6))
axs[-1].set_xlabel("position along the centromere (kb from the first unit); shaded = arrays")
fig.suptitle("Divergence of each unit (element, alignment K2P over shared homologous windows) at four levels", fontsize=13, fontweight="bold", x=0.02, ha="left")
fig.tight_layout()
for ext in ("png", "pdf"): fig.savefig(W / "figure" / f"per_unit_divergence.{ext}", dpi=450, bbox_inches="tight")
plt.close(fig)
# figure 2: k-mer grid
cols = ["el_align_k2p"] + [f"el_k{k}_{m}" for m in ("mash", "containment", "wjaccard", "cosine") for k in KS]
M = R.groupby("level")[cols].median().reindex(LV)
fig, ax = plt.subplots(figsize=(15, 4.5), facecolor="white")
Mn = M / M.loc[LV[4]]
im = ax.imshow(Mn.values, aspect="auto", cmap="viridis"); ax.set_yticks(range(len(LV))); ax.set_yticklabels(LV, fontsize=9)
ax.set_xticks(range(len(cols))); ax.set_xticklabels([c.replace("el_", "") for c in cols], rotation=60, ha="right", fontsize=8)
for i_ in range(len(LV)):
    for j_ in range(len(cols)): ax.text(j_, i_, f"{M.values[i_, j_]:.3f}", ha="center", va="center", fontsize=6.5, color="white" if Mn.values[i_, j_] < 0.6 else "black")
fig.colorbar(im, ax=ax, label="relative to Gar level (= 1)")
ax.set_title("Element divergence per level across methods (median over units; colour scaled to the between-species level)", loc="left", fontsize=12, fontweight="bold")
for ext in ("png", "pdf"): fig.savefig(W / "figure" / f"per_unit_kmer_grid.{ext}", dpi=450, bbox_inches="tight")
