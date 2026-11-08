"""Divergence along the unit, part by part (2026-10-23; user: establish divergence in some parts and less in centromeres).
Parts: MacSat03 satellite block of the unit (satellite runs, concatenated), element in 500-bp windows of the element consensus (unit sequence mapped
through the element hits), type II insert (insert hits). For each part, k-mer distance (k = 11, Mash formula, full sets) for up to 300 random pairs in:
  within centromere | between centromeres, same type | between type I and II | dispersed arm copies (element windows only; copies outside centromeres).
Outputs: kmer_tree/divergence_by_part.tsv, figure/divergence_by_part.png/.pdf"""
import pathlib, itertools, random, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/code"); import ltr_colours as LC

W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22"); A = pathlib.Path("F:/data/fClaHyb_func_annotation_v1.0")
R = pathlib.Path("F:/data/satellite_FISH/runs_of_record_v4.5.3_2026-10-18/Mac"); K = 11; NP = 300; random.seed(20261023)
T1 = ["C_mac_LG07", "C_mac_LG04", "C_mac_LG01", "C_mac_LG23", "C_mac_LG02"]; T2 = ["C_mac_LG08", "C_mac_LG25", "C_mac_LG16"]
class Fa:
    def __init__(self, p):
        self.f = open(p, "rb"); self.i = {}
        for l in open(str(p) + ".fai"):
            n, ln, off, lb, lby = l.split("\t")[:5]; self.i[n] = (int(ln), int(off), int(lb), int(lby))
    def get(self, c, a, b):
        ln, off, lb, lby = self.i[c]; a = max(0, a); b = min(ln, b)
        s0 = off + (a // lb) * lby + a % lb; s1 = off + (b // lb) * lby + b % lb
        self.f.seek(s0); return self.f.read(s1 - s0).replace(b"\n", b"").replace(b"\r", b"").decode().upper()
fa = Fa(A / "fClaHyb_Mac/fClaHyb_Mac_genome/fClaHyb_Mac_LG.fa"); comp = str.maketrans("ACGT", "TGCA")
def kset(s):
    return {min(s[i:i + K], s[i:i + K].translate(comp)[::-1]) for i in range(len(s) - K + 1) if "N" not in s[i:i + K]}
def mash(a, b):
    if not a or not b: return np.nan
    j = len(a & b) / len(a | b); return 1.0 if j == 0 else max(0.0, -1 / K * np.log(2 * j / (1 + j)))

U = pd.read_csv(W / "unit" / "units.tsv", sep="\t"); U = U[U.chrom.isin(T1 + T2)].reset_index(drop=True)
runs = pd.read_csv(R / "runs_of_record.tsv", sep="\t"); runs = runs[(runs.unit_class == "satellite") & runs.tier.isin(["S1", "S2"]) & (runs.family == "MacSat03")]
eh = pd.read_csv(W / "copies" / "Mac.hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
eh = eh[~((eh.q0 >= 2900) & (eh.q1 <= 3420))]; eh["lo"] = eh[["s0", "s1"]].min(1); eh["hi"] = eh[["s0", "s1"]].max(1)
ih = pd.read_csv(W / "insert" / "Mac.insert_hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
ih["lo"] = ih[["s0", "s1"]].min(1); ih["hi"] = ih[["s0", "s1"]].max(1)
WIN = [(w, w + 500) for w in range(0, 8164, 500)]

def window_seq(h, w0, w1):
    """genome sequence of the element hits h mapped onto element window [w0, w1) (linear interpolation inside each hit)"""
    parts = []
    for r in h.itertuples():
        a, b = max(r.q0, w0), min(r.q1, w1)
        if b - a < 60: continue
        f = lambda q: r.s0 + (q - r.q0) * (r.s1 - r.s0) / max(1, r.q1 - r.q0)
        x0, x1 = sorted((int(f(a)), int(f(b)))); parts.append(fa.get(r.chrom, x0, x1))
    return "".join(parts)

# per unit: parts
items = []
for r in U.itertuples():
    a, b = int(r.start), int(r.start + r.period)
    sat = runs[(runs.chrom == r.chrom) & (runs.end > a) & (runs.start0 < b)]
    rec = dict(group=r.chrom, type="I" if r.chrom in T1 else "II", kind="cen",
               sat=kset("".join(fa.get(r.chrom, max(s, a), min(e, b)) for s, e in zip(sat.start0, sat.end))))
    h = eh[(eh.chrom == r.chrom) & (eh.hi > a) & (eh.lo < b)]
    for w0, w1 in WIN: rec[f"el_{w0}"] = kset(window_seq(h, w0, w1))
    hi = ih[(ih.chrom == r.chrom) & (ih.hi > a) & (ih.lo < b) & (ih.len >= 300)]
    rec["insert"] = kset("".join(fa.get(r.chrom, max(s, a), min(e, b)) for s, e in zip(hi.lo, hi.hi)))
    items.append(rec)
# dispersed arm copies: element hits within +-8 kb of each arm core copy
cp = pd.read_csv(W / "core" / "Mac.core_copies.tsv", sep="\t"); arm = cp[~cp.in_cen_1Mb]
for r in arm.itertuples():
    h = eh[(eh.chrom == r.chrom) & (eh.hi > r.lo - 8000) & (eh.lo < r.hi + 8000)]
    rec = dict(group=f"arm_{r.chrom}_{r.lo}", type="-", kind="arm", sat=set(), insert=set())
    for w0, w1 in WIN: rec[f"el_{w0}"] = kset(window_seq(h, w0, w1))
    items.append(rec)

cen = [i for i, x in enumerate(items) if x["kind"] == "cen"]; armi = [i for i, x in enumerate(items) if x["kind"] == "arm"]
def pairs(cat):
    if cat == "within centromere": P = [(i, j) for i, j in itertools.combinations(cen, 2) if items[i]["group"] == items[j]["group"]]
    elif cat == "between centromeres, same type": P = [(i, j) for i, j in itertools.combinations(cen, 2) if items[i]["group"] != items[j]["group"] and items[i]["type"] == items[j]["type"]]
    elif cat == "between type I and II": P = [(i, j) for i, j in itertools.combinations(cen, 2) if items[i]["type"] != items[j]["type"]]
    else: P = list(itertools.combinations(armi, 2))
    return random.sample(P, min(NP, len(P)))
CATS = ["within centromere", "between centromeres, same type", "between type I and II", "dispersed copies (outside centromeres)"]
PARTS = ["sat"] + [f"el_{w0}" for w0, _ in WIN] + ["insert"]
rows = []
for cat in CATS:
    P = pairs(cat)
    for p in PARTS:
        if p in ("sat", "insert") and cat.startswith("dispersed"): continue
        d = [mash(items[i][p], items[j][p]) for i, j in P if len(items[i][p]) > 50 and len(items[j][p]) > 50]
        if p == "insert": d = [x for x, (i, j) in zip(d, P) if items[i]["type"] == "II" and items[j]["type"] == "II"] if d else d
        rows.append(dict(category=cat, part=p, pairs=len(d), median=np.nanmedian(d) if d else np.nan,
                         q25=np.nanpercentile(d, 25) if d else np.nan, q75=np.nanpercentile(d, 75) if d else np.nan))
T = pd.DataFrame(rows); T.round(4).to_csv(W / "kmer_tree" / "divergence_by_part.tsv", sep="\t", index=False)

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(17, 6.5), gridspec_kw=dict(width_ratios=[4, 1.1]), facecolor="white")
COL = dict(zip(CATS, [LC.REL[0], LC.REL[2], LC.REL[3], LC.ELEM]))
x = [w0 / 1000 + 0.25 for w0, _ in WIN]
for cat in CATS:
    t = T[(T.category == cat) & T.part.str.startswith("el_")].set_index("part").reindex([f"el_{w0}" for w0, _ in WIN])
    ok = t.pairs.fillna(0).values >= 10
    m = np.where(ok, t["median"].values, np.nan); lo = np.where(ok, t.q25.values, np.nan); hi = np.where(ok, t.q75.values, np.nan)   # NaN breaks the line where data are missing
    ax0.plot(x, m, "-o", color=COL[cat], label=cat, lw=2, ms=5)
    ax0.fill_between(x, lo, hi, color=LC.P.shade(COL[cat], 0.8), zorder=0)
ax0.axvspan(5.56, 7.4, color=LC.P.shade(LC.ELEM_CORE, 0.85), zorder=0); ax0.text(6.48, 0.325, "RT core", ha="center", fontsize=9, color=LC.ELEM_CORE, va="top"); ax0.text(3.16, 0.325, "built-in satellite", ha="center", fontsize=9, color=LC.SAT, va="top")
ax0.axvspan(2.93, 3.39, color=LC.P.shade(LC.SAT, 0.85), zorder=0)
ax0.set_xlabel("element position (kb, 500-bp windows)", fontsize=11); ax0.set_ylabel("k-mer divergence (median, IQR)", fontsize=11)
ax0.set_title("A  Divergence along the element, by comparison", loc="left", fontsize=13, fontweight="bold"); ax0.legend(fontsize=9, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=4)
for s in ("top", "right"): ax0.spines[s].set_visible(False)
xs = np.arange(2); w = 0.25
for k, cat in enumerate(CATS[:3]):
    v = [T[(T.category == cat) & (T.part == p)]["median"].values for p in ("sat", "insert")]
    v = [x_[0] if len(x_) else np.nan for x_ in v]
    ax1.bar(xs + (k - 1) * w, v, width=w, color=COL[cat], label=cat)
ax1.set_xticks(xs); ax1.set_xticklabels(["MacSat03\nsatellite", "type II\ninsert"], fontsize=10)
ax1.set_title("B  Satellite and insert", loc="left", fontsize=13, fontweight="bold")
for s in ("top", "right"): ax1.spines[s].set_visible(False)
fig.suptitle("Where the unit diverges: part by part, within vs between centromeres, and vs dispersed copies", fontsize=14, fontweight="bold", x=0.02, ha="left")
for ext in ("png", "pdf"): fig.savefig(W / "figure" / f"divergence_by_part.{ext}", dpi=450, bbox_inches="tight")
print(T.pivot_table(index="part", columns="category", values="median").reindex(PARTS).round(3).to_string())
