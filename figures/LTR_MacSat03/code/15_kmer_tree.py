"""k-mer divergence trees (2026-10-23; user: k-mer tree to see older/younger copies across chromosomes, a divergence proposal).
1 Unit tree: every Mac centromeric unit (units.tsv; 8-25 kb) -> canonical 15-mer sets -> Mash distance D = -1/k ln(2J/(1+J)) -> NJ, midpoint root.
2 Element-core tree: Mac core copies (centromeric + arm) + 120 random Gar core copies -> 15-mer Mash -> NJ rooted on the Gar copies.
  Relative age of each Mac centromere's copies = median root-to-tip depth of its copies below the node where Mac separates from Gar
  (age 1 = as deep as that node; 0 = identical copies). Within-group diversity = median pairwise Mash distance.
Outputs: kmer_tree/units.nwk, units_by_chrom.tsv, core.nwk, core_ages.tsv, figure/kmer_divergence.png/.pdf"""
import pathlib, random, itertools, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/code"); import ltr_colours as LC
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22"); A = pathlib.Path("F:/data/fClaHyb_func_annotation_v1.0")
OUT = W / "kmer_tree"; OUT.mkdir(exist_ok=True); K = 15; random.seed(20261023)
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
        if "N" in k: continue
        r = k.translate(comp)[::-1]; out.add(hash(min(k, r)))
    return out
def mash(a, b):
    j = len(a & b) / len(a | b) if a or b else 0
    return 1.0 if j == 0 else max(0.0, -1 / K * np.log(2 * j / (1 + j)))
def nj(names, sets):
    n = len(names); M = [[0.0] * (i + 1) for i in range(n)]
    for i in range(n):
        for j in range(i): M[i][j] = mash(sets[i], sets[j])
    return DistanceTreeConstructor().nj(DistanceMatrix(names, M)), M

macfa = Fa(A / "fClaHyb_Mac/fClaHyb_Mac_genome/fClaHyb_Mac_LG.fa"); garfa = Fa(A / "fClaHyb_Gar/fClaHyb_Gar_genome/fClaHyb_Gar_LG.fa")

# ---- 1 unit tree
U = pd.read_csv(W / "unit" / "units.tsv", sep="\t"); U = U[U.chrom.isin(T1 + T2)].reset_index(drop=True)
un = [f"{lab(r.chrom)}_{r.start}" for r in U.itertuples()]
us = [kset(macfa.get(r.chrom, int(r.start), int(r.start + r.period))) for r in U.itertuples()]
ut, UM = nj(un, us); ut.root_at_midpoint(); Phylo.write(ut, OUT / "units.nwk", "newick")
chrom_of = {n: "C_mac_" + n.split("_")[0] for n in un}
rows = []
for c in T1 + T2:
    idx = [i for i, n in enumerate(un) if chrom_of[n] == c]
    within = [UM[max(i, j)][min(i, j)] for i, j in itertools.combinations(idx, 2)]
    other = [UM[max(i, j)][min(i, j)] for i in idx for j in range(len(un)) if chrom_of[un[j]] != c]
    same_t = [UM[max(i, j)][min(i, j)] for i in idx for j in range(len(un)) if chrom_of[un[j]] != c and (chrom_of[un[j]] in T1) == (c in T1)]
    rows.append(dict(chrom=lab(c), type="I" if c in T1 else "II", units=len(idx), within_centromere=np.median(within),
                     to_same_type=np.median(same_t), to_all_others=np.median(other)))
UC = pd.DataFrame(rows); UC.round(4).to_csv(OUT / "units_by_chrom.tsv", sep="\t", index=False)

# ---- 2 element-core tree with Gar root
def core_seqs(h, fa, n_sample=None):
    c = pd.read_csv(W / "core" / f"{h}.core_copies.tsv", sep="\t")
    if n_sample: c = c.sample(n=min(n_sample, len(c)), random_state=1)
    return [(f"{h}|{lab(r.chrom)}|{r.lo}|{'cen' if r.in_cen_1Mb else 'arm'}", fa.get(r.chrom, int(r.lo), int(r.hi))) for r in c.itertuples()]
cs = core_seqs("Mac", macfa) + core_seqs("Gar", garfa, 120)
cn = [n for n, _ in cs]; cset = [kset(s) for _, s in cs]
ct, CM = nj(cn, cset)
gar = [t for t in ct.get_terminals() if t.name.startswith("Gar")]
ct.root_with_outgroup(ct.common_ancestor(gar)) if len(gar) > 1 else None
Phylo.write(ct, OUT / "core.nwk", "newick")
# relative age per group (2026-10-23 fix: MRCA-based depth failed because groups are not monophyletic):
# within-group diversity = median pairwise Mash distance among the group's core copies; distance to Gar = median distance to the Gar sample.
# Scaled by the median Mac-vs-Gar core distance (1 = as diverged as copies of the two species).
D = lambda i, j: CM[max(i, j)][min(i, j)]
gar_i = [i for i, n in enumerate(cn) if n.startswith("Gar")]
macgar = np.median([D(i, j) for i in range(len(cn)) if cn[i].startswith("Mac") for j in gar_i])
groups = [(lab(c), "I" if c in T1 else "II", [i for i, n in enumerate(cn) if n.startswith("Mac") and n.split("|")[1] == lab(c) and n.endswith("cen")]) for c in T1 + T2]
groups += [("Mac arm copies", "-", [i for i, n in enumerate(cn) if n.startswith("Mac") and n.endswith("arm")]), ("Gar copies (sample)", "-", gar_i)]
rows = []
for name, typ, idx in groups:
    if len(idx) < 3: continue
    w = np.median([D(i, j) for i, j in itertools.combinations(idx, 2)])
    g = np.median([D(i, j) for i in idx for j in gar_i if i != j]) if not name.startswith("Gar") else np.nan
    rows.append(dict(chrom=name, type=typ, copies=len(idx), within_group=w, to_Gar=g, relative_age_of_expansion=w / macgar))
CA = pd.DataFrame(rows); CA.round(4).to_csv(OUT / "core_ages.tsv", sep="\t", index=False)
print("median Mac-Gar core distance (scale = 1):", round(macgar, 4))

# ---- figure
fig = plt.figure(figsize=(17, 14), facecolor="white"); gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1.0], hspace=0.3, wspace=0.25)
def draw_tree(ax, tree, colour_of, title):
    depths = tree.depths(); tips = tree.get_terminals(); y = {t: i for i, t in enumerate(tips)}
    def place(cl):
        if cl.is_terminal(): return y[cl]
        ys = [place(c) for c in cl.clades]; y[cl] = (min(ys) + max(ys)) / 2; return y[cl]
    place(tree.root)
    for cl in tree.find_clades():
        for ch in cl.clades:
            ax.plot([depths[cl], depths[cl]], [y[cl], y[ch]], color="#999", lw=0.5); ax.plot([depths[cl], depths[ch]], [y[ch], y[ch]], color="#999", lw=0.5)
    for t in tips: ax.scatter(depths[t], y[t], s=10, color=colour_of(t.name), zorder=3, edgecolor="none")
    ax.set_yticks([]); ax.set_xlabel("k-mer (Mash) distance from the root", fontsize=10); ax.set_title(title, loc="left", fontsize=13, fontweight="bold")
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax = fig.add_subplot(gs[:, 0])
draw_tree(ax, ut, lambda n: CH_COL[chrom_of[n]], "A  Unit tree: 177 whole units from 8 Mac centromeres (15-mer Mash, NJ)")
for c in T1 + T2: ax.scatter([], [], color=CH_COL[c], s=40, label=f"{lab(c)} (type {'I' if c in T1 else 'II'})")
ax.legend(fontsize=9, frameon=False, loc="lower left")
ax = fig.add_subplot(gs[0, 1]); x = np.arange(len(UC))
for k, (col, nm) in enumerate((("within_centromere", "within its centromere"), ("to_same_type", "to other centromeres, same type"), ("to_all_others", "to all other centromeres"))):
    ax.bar(x + (k - 1) * 0.27, UC[col], width=0.27, label=nm, color=[LC.REL[0], LC.REL[2], LC.REL[3]][k])
ax.set_xticks(x); ax.set_xticklabels([f"{r.chrom}\n{r.type}" for r in UC.itertuples()], fontsize=9); ax.set_ylabel("median unit-to-unit Mash distance", fontsize=10)
ax.legend(fontsize=9, frameon=False); ax.set_title("B  Units: diversity within vs between centromeres", loc="left", fontsize=13, fontweight="bold")
ax = fig.add_subplot(gs[1, 1]); ca = CA.dropna(subset=["relative_age_of_expansion"])
ax.barh(range(len(ca)), ca.relative_age_of_expansion, color=[CH_COL.get("C_mac_" + c, "#999") for c in ca.chrom])
ax.set_yticks(range(len(ca))); ax.set_yticklabels([f"{r.chrom} ({r.type})" if r.type != "-" else r.chrom for r in ca.itertuples()], fontsize=10)
ax.axvline(1, color="black", lw=1, ls="--"); ax.text(1.0, -0.6, " Mac vs Gar copies", fontsize=8, va="bottom")
ax.set_xlabel("diversity among copies, relative to Mac-vs-Gar (higher = older)", fontsize=10); ax.invert_yaxis()
ax.set_title("C  Element core: diversity per group (age proxy)", loc="left", fontsize=13, fontweight="bold")
fig.suptitle("k-mer divergence of the composite units and their element copies", fontsize=16, fontweight="bold", x=0.02, ha="left")
for ext in ("png", "pdf"): fig.savefig(W / "figure" / f"kmer_divergence.{ext}", dpi=450, bbox_inches="tight")
print(UC.round(3).to_string(index=False)); print(); print(CA.round(3).to_string(index=False))
