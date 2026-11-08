"""Normalized homogenization and element age (2026-10-23; user: B/C mix things; normalize by number of units).
Figure 1 homogenization (whole units only): rarefy every centromere to N_MIN units (smallest centromere), 200 replicates; per centromere
  within = median Mash distance among its sampled units; between = median distance to the sampled units of all other centromeres;
  H = 1 - within/between (1 = units of a centromere identical, 0 = no more alike than units of other centromeres). Medians and 95% intervals.
  Plus the unit tree (from kmer_tree/units.nwk).
Figure 2 age of the element (element core copies): Mac arm copies, Gar copies, Mac centromeric copies pooled; equal sample size, 200 replicates;
  diversity = median pairwise Mash distance (higher = older relicts; centromeric = homogenized).
Outputs: figure/homogenization.png/.pdf, figure/element_age.png/.pdf, kmer_tree/homogenization_rarefied.tsv, kmer_tree/element_age_rarefied.tsv"""
import pathlib, itertools, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/code"); import ltr_colours as LC
from Bio import Phylo

W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22"); A = pathlib.Path("F:/data/fClaHyb_func_annotation_v1.0"); OUT = W / "kmer_tree"
K = 15; rng = np.random.default_rng(20261023); REP = 200
T1 = ["C_mac_LG07", "C_mac_LG04", "C_mac_LG01", "C_mac_LG23", "C_mac_LG02"]; T2 = ["C_mac_LG08", "C_mac_LG25", "C_mac_LG16"]
CH_COL = dict(zip(T1 + T2, LC.CHROM_SEQ))
lab = lambda c: "LG" + c.split("LG")[-1]
class Fa:
    def __init__(self, p):
        self.f = open(p, "rb"); self.i = {}
        for l in open(str(p) + ".fai"):
            n, ln, off, lb, lby = l.split("\t")[:5]; self.i[n] = (int(ln), int(off), int(lb), int(lby))
    def get(self, c, a, b):
        ln, off, lb, lby = self.i[c]; s0 = off + (a // lb) * lby + a % lb; s1 = off + (b // lb) * lby + b % lb
        self.f.seek(s0); return self.f.read(s1 - s0).replace(b"\n", b"").replace(b"\r", b"").decode().upper()
comp = str.maketrans("ACGT", "TGCA")
def kset(s):
    out = set()
    for i in range(len(s) - K + 1):
        k = s[i:i + K]
        if "N" not in k: r = k.translate(comp)[::-1]; out.add(h64(min(k, r)))
    return out
import hashlib
def h64(k): return int.from_bytes(hashlib.blake2b(k.encode(), digest_size=8).digest(), "little")   # stable hash (Python hash() is salted per run)
def sketch(ks, size):
    """fixed-size bottom-s MinHash sketch: every sequence contributes the same number of k-mers (normalizes for length; 2026-10-23)"""
    return frozenset(sorted(ks)[:size])
def mash(a, b):
    s = max(len(a), len(b)); u = sorted(a | b)[:s]; inter = len(set(u) & a & b)     # Mash estimator on bottom-s of the union
    j = inter / s if s else 0; return 1.0 if j == 0 else max(0.0, -1 / K * np.log(2 * j / (1 + j)))
macfa = Fa(A / "fClaHyb_Mac/fClaHyb_Mac_genome/fClaHyb_Mac_LG.fa"); garfa = Fa(A / "fClaHyb_Gar/fClaHyb_Gar_genome/fClaHyb_Gar_LG.fa")

# ---------- Figure 1: homogenization, units only, rarefied
U = pd.read_csv(W / "unit" / "units.tsv", sep="\t"); U = U[U.chrom.isin(T1 + T2)].reset_index(drop=True)
S = [sketch(kset(macfa.get(r.chrom, int(r.start), int(r.start + r.period))), 2000) for r in U.itertuples()]
n = len(U); D = np.zeros((n, n))
for i, j in itertools.combinations(range(n), 2): D[i, j] = D[j, i] = mash(S[i], S[j])
idx = {c: np.where(U.chrom.values == c)[0] for c in T1 + T2}; N_MIN = min(len(v) for v in idx.values())
res = {c: {"within": [], "between": [], "H": []} for c in T1 + T2}
for _ in range(REP):
    samp = {c: rng.choice(v, N_MIN, replace=False) for c, v in idx.items()}
    for c in T1 + T2:
        s = samp[c]; w = np.median([D[i, j] for i, j in itertools.combinations(s, 2)])
        o = np.concatenate([samp[x] for x in T1 + T2 if x != c]); b = np.median(D[np.ix_(s, o)])
        res[c]["within"].append(w); res[c]["between"].append(b); res[c]["H"].append(1 - w / b)
rows = []
for c in T1 + T2:
    r = {k: np.array(v) for k, v in res[c].items()}
    rows.append(dict(chrom=lab(c), type="I" if c in T1 else "II", units=len(idx[c]), sampled=N_MIN,
                     within=np.median(r["within"]), between=np.median(r["between"]), H=np.median(r["H"]),
                     H_lo=np.percentile(r["H"], 2.5), H_hi=np.percentile(r["H"], 97.5),
                     within_lo=np.percentile(r["within"], 2.5), within_hi=np.percentile(r["within"], 97.5)))
HT = pd.DataFrame(rows); HT.round(4).to_csv(OUT / "homogenization_rarefied.tsv", sep="\t", index=False)

fig = plt.figure(figsize=(16, 9), facecolor="white"); gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1.0], wspace=0.25)
ax = fig.add_subplot(gs[0, 0]); t = Phylo.read(str(OUT / "units.nwk"), "newick")
depths = t.depths(); tips = t.get_terminals(); y = {tp: i for i, tp in enumerate(tips)}
def place(cl):
    if cl.is_terminal(): return y[cl]
    ys = [place(c) for c in cl.clades]; y[cl] = (min(ys) + max(ys)) / 2; return y[cl]
place(t.root)
for cl in t.find_clades():
    for ch in cl.clades:
        ax.plot([depths[cl]] * 2, [y[cl], y[ch]], color="#999", lw=0.5); ax.plot([depths[cl], depths[ch]], [y[ch]] * 2, color="#999", lw=0.5)
for tp in tips: ax.scatter(depths[tp], y[tp], s=10, color=CH_COL["C_mac_" + tp.name.split("_")[0]], zorder=3, edgecolor="none")
for c in T1 + T2:
    ys_ = [y[tp] for tp in tips if tp.name.split("_")[0] == lab(c)]
    ax.text(max(depths.values()) + 0.003, np.mean(ys_), f"{lab(c)} ({'I' if c in T1 else 'II'})", va="center", fontsize=10, color=CH_COL[c], fontweight="bold")
ax.set_yticks([]); ax.set_xlabel("k-mer (Mash) distance from the root", fontsize=10)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.set_title("A  Tree of the 177 whole units (each tip = one unit)", loc="left", fontsize=13, fontweight="bold")
ax = fig.add_subplot(gs[0, 1]); x = np.arange(len(HT))
ax.bar(x, HT.H, color=[CH_COL["C_mac_" + c] for c in HT.chrom], width=0.65)
ax.errorbar(x, HT.H, yerr=[HT.H - HT.H_lo, HT.H_hi - HT.H], fmt="none", ecolor="black", capsize=3, lw=1)
for i, r in enumerate(HT.itertuples()): ax.text(i, 0.02, f"{r.units}\nunits", ha="center", fontsize=8, color="white", fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels([f"{r.chrom}\n{r.type}" for r in HT.itertuples()], fontsize=10); ax.set_ylim(0, 1.05)
ax.set_ylabel("homogenization index H = 1 − within / between", fontsize=10)
ax.set_title(f"B  Homogenization per centromere\n({N_MIN} units per centromere, 2,000-k-mer sketch per unit, {REP} replicates, 95% interval)", loc="left", fontsize=13, fontweight="bold")
ax.axhline(1, color="#888", lw=0.8, ls="--")
fig.suptitle("Homogenization of the composite units within each centromere", fontsize=15, fontweight="bold", x=0.02, ha="left")
for ext in ("png", "pdf"): fig.savefig(W / "figure" / f"homogenization.{ext}", dpi=450, bbox_inches="tight")
plt.close(fig)

# ---------- Figure 2: age of the element, equal sample sizes
mc = pd.read_csv(W / "core" / "Mac.core_copies.tsv", sep="\t"); gc = pd.read_csv(W / "core" / "Gar.core_copies.tsv", sep="\t")
groups = {"Mac, outside centromeres": [macfa.get(r.chrom, int(r.lo), int(r.hi)) for r in mc[~mc.in_cen_1Mb].itertuples()],
          "Gar, all copies": [garfa.get(r.chrom, int(r.lo), int(r.hi)) for r in gc.sample(n=min(300, len(gc)), random_state=1).itertuples()],
          "Mac, in centromeres (pooled)": [macfa.get(r.chrom, int(r.lo), int(r.hi)) for r in mc[mc.in_cen_1Mb].itertuples()]}
gs_ = {g: [sketch(kset(s), 500) for s in v] for g, v in groups.items()}; NS = min(60, min(len(v) for v in gs_.values()))
age = {g: [] for g in gs_}
for _ in range(REP):
    for g, v in gs_.items():
        s = rng.choice(len(v), NS, replace=False); age[g].append(np.median([mash(v[i], v[j]) for i, j in itertools.combinations(s, 2)]))
AT = pd.DataFrame([dict(group=g, copies=len(gs_[g]), sampled=NS, diversity=np.median(v), lo=np.percentile(v, 2.5), hi=np.percentile(v, 97.5)) for g, v in age.items()])
AT.round(4).to_csv(OUT / "element_age_rarefied.tsv", sep="\t", index=False)
fig, ax = plt.subplots(figsize=(10, 5), facecolor="white")
ax.barh(range(len(AT)), AT.diversity, color=[LC.REL[2], LC.REL[3], LC.ELEM], height=0.6)
ax.errorbar(AT.diversity, range(len(AT)), xerr=[AT.diversity - AT.lo, AT.hi - AT.diversity], fmt="none", ecolor="black", capsize=4)
ax.set_yticks(range(len(AT))); ax.set_yticklabels([f"{r.group}\n({r.copies} copies)" for r in AT.itertuples()], fontsize=11); ax.invert_yaxis()
ax.set_xlabel(f"diversity among element copies (median pairwise k-mer distance; {NS} copies per group, 500-k-mer sketch per copy, {REP} replicates)", fontsize=10)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.set_title("Diversity of element copies: Mac dispersed copies are the oldest; centromeric pool mixes 8 centromeres", loc="left", fontsize=12, fontweight="bold")
for ext in ("png", "pdf"): fig.savefig(W / "figure" / f"element_age.{ext}", dpi=450, bbox_inches="tight")
print(HT.round(3).to_string(index=False)); print(); print(AT.round(4).to_string(index=False))
