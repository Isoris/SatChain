"""Bootstrap stability of the unit tree (2026-10-23; user: is each unit stable against the tree topology; bootstrapping).
Units -> canonical 15-mers -> presence matrix over a 1/50 hash subsample of all k-mers (columns). 200 replicates: resample columns with replacement,
weighted Jaccard -> Mash distance -> NJ (midpoint root). Per replicate: monophyly of each centromere, of type I and of type II, sister pairs
(LG07+LG01, LG08+LG25), and per unit: (a) its sister group holds only units of its own centromere, (b) nearest unit (tree distance)
from the same centromere. Outputs: kmer_tree/bootstrap_units.tsv (per unit), kmer_tree/bootstrap_clades.tsv, figure/bootstrap.png/.pdf"""
import pathlib, hashlib, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/code"); import ltr_colours as LC
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22"); A = pathlib.Path("F:/data/fClaHyb_func_annotation_v1.0")
K, SUB, REP = 15, 50, 200; rng = np.random.default_rng(20261023)
T1 = ["C_mac_LG07", "C_mac_LG04", "C_mac_LG01", "C_mac_LG23", "C_mac_LG02"]; T2 = ["C_mac_LG08", "C_mac_LG25", "C_mac_LG16"]
lab = lambda c: "LG" + c.split("LG")[-1]
class Fa:
    def __init__(self, p):
        self.f = open(p, "rb"); self.i = {}
        for l in open(str(p) + ".fai"):
            n, ln, off, lb, lby = l.split("\t")[:5]; self.i[n] = (int(ln), int(off), int(lb), int(lby))
    def get(self, c, a, b):
        ln, off, lb, lby = self.i[c]; s0 = off + (a // lb) * lby + a % lb; s1 = off + (b // lb) * lby + b % lb
        self.f.seek(s0); return self.f.read(s1 - s0).replace(b"\n", b"").replace(b"\r", b"").decode().upper()
fa = Fa(A / "fClaHyb_Mac/fClaHyb_Mac_genome/fClaHyb_Mac_LG.fa"); comp = str.maketrans("ACGT", "TGCA")
def h64(k): return int.from_bytes(hashlib.blake2b(k.encode(), digest_size=8).digest(), "little")
def ksub(s):
    out = set()
    for i in range(len(s) - K + 1):
        k = s[i:i + K]
        if "N" in k: continue
        v = h64(min(k, k.translate(comp)[::-1]))
        if v % SUB == 0: out.add(v)
    return out
U = pd.read_csv(W / "unit" / "units.tsv", sep="\t"); U = U[U.chrom.isin(T1 + T2)].reset_index(drop=True)
names = [f"{lab(r.chrom)}_{r.start}" for r in U.itertuples()]; chrom = [lab(c) for c in U.chrom]
sets = [ksub(fa.get(r.chrom, int(r.start), int(r.start + r.period))) for r in U.itertuples()]
cols = sorted(set().union(*sets)); ci = {v: i for i, v in enumerate(cols)}
X = np.zeros((len(sets), len(cols)), dtype=np.float32)
for i, s in enumerate(sets): X[i, [ci[v] for v in s]] = 1
print("units", len(names), "k-mer columns", len(cols), flush=True)

def tree_from(w):
    Xw = X * w; inter = Xw @ X.T; tot = Xw.sum(1); union = tot[:, None] + tot[None, :] - inter
    J = np.where(union > 0, inter / union, 0); D = np.where(J > 0, -1 / K * np.log(np.clip(2 * J / (1 + J), 1e-12, 1)), 1.0); np.fill_diagonal(D, 0)
    M = [[float(max(D[i, j], 0)) for j in range(i + 1)] for i in range(len(names))]
    t = DistanceTreeConstructor().nj(DistanceMatrix(names, M)); t.root_at_midpoint(); return t

GROUPS = {lab(c): {n for n, ch in zip(names, chrom) if ch == lab(c)} for c in T1 + T2}
GROUPS["type I"] = {n for n, ch in zip(names, chrom) if "C_mac_" + ch in T1}; GROUPS["type II"] = {n for n, ch in zip(names, chrom) if "C_mac_" + ch in T2}
GROUPS["LG07+LG01"] = GROUPS["LG07"] | GROUPS["LG01"]; GROUPS["LG08+LG25"] = GROUPS["LG08"] | GROUPS["LG25"]
clade_hits = {g: 0 for g in GROUPS}; in_own = np.zeros(len(names)); nn_same = np.zeros(len(names))
for r in range(REP + 1):
    w = np.ones(len(cols), dtype=np.float32) if r == 0 else np.bincount(rng.integers(0, len(cols), len(cols)), minlength=len(cols)).astype(np.float32)
    t = tree_from(w)
    if r == 0: Phylo.write(t, W / "kmer_tree" / "units_bootref.nwk", "newick"); continue
    splits = [frozenset(x.name for x in cl.get_terminals()) for cl in t.get_nonterminals()]
    sset = set(splits); alln = frozenset(names)
    for g, members in GROUPS.items():
        fm = frozenset(members); clade_hits[g] += (fm in sset) or ((alln - fm) in sset)      # monophyly (either side of the root)
    # per unit: its sister group (smallest clade holding it and >= 1 other unit) contains only units of its own centromere
    for cl in t.get_nonterminals():
        kids = cl.clades
        for kid in kids:
            if kid.is_terminal():
                sib = {x.name for k2 in kids if k2 is not kid for x in k2.get_terminals()}
                own = kid.name.split("_")[0]
                if sib and all(x.split("_")[0] == own for x in sib): in_own[names.index(kid.name)] += 1
    if r % 20 == 0: print("replicate", r, flush=True)
# nearest-neighbour stability on the distance matrices directly (cheaper than tree distances, same question)
for r in range(REP):
    w = np.bincount(rng.integers(0, len(cols), len(cols)), minlength=len(cols)).astype(np.float32)
    Xw = X * w; inter = Xw @ X.T; tot = Xw.sum(1); union = tot[:, None] + tot[None, :] - inter; J = inter / np.maximum(union, 1e-9)
    np.fill_diagonal(J, -1); nn = J.argmax(1)
    nn_same += np.array([chrom[nn[i]] == chrom[i] for i in range(len(names))])
PU = pd.DataFrame(dict(unit=names, centromere=chrom, type=["I" if "C_mac_" + c in T1 else "II" for c in chrom],
                       sister_group_same_centromere=in_own / REP, nearest_unit_same_centromere=nn_same / REP))
PU.to_csv(W / "kmer_tree" / "bootstrap_units.tsv", sep="\t", index=False)
CS = pd.DataFrame([dict(group=g, n_units=len(m), bootstrap_support=clade_hits[g] / REP) for g, m in GROUPS.items()])
CS.to_csv(W / "kmer_tree" / "bootstrap_clades.tsv", sep="\t", index=False)

fig, axs = plt.subplots(1, 2, figsize=(16, 6.5), gridspec_kw=dict(width_ratios=[1.0, 1.6]), facecolor="white")
ax = axs[0]; cs = CS.iloc[::-1]
ax.barh(range(len(cs)), cs.bootstrap_support * 100, color=[LC.TYPE2 if g in ("LG08", "LG25", "LG16", "type II", "LG08+LG25") else LC.TYPE1 for g in cs.group])
ax.set_yticks(range(len(cs))); ax.set_yticklabels([f"{g} ({n})" for g, n in zip(cs.group, cs.n_units)], fontsize=10)
ax.axvline(70, color="#888", ls="--", lw=1); ax.set_xlim(0, 105); ax.set_xlabel("bootstrap support (%, 200 replicates)", fontsize=10)
for i, v in enumerate(cs.bootstrap_support * 100): ax.text(v + 1, i, f"{v:.0f}", va="center", fontsize=9)
ax.set_title("A  Clade support (k-mer bootstrap)", loc="left", fontsize=13, fontweight="bold")
for s in ("top", "right"): axs[0].spines[s].set_visible(False)
ax = axs[1]; order = [lab(c) for c in T1 + T2]
for k, c in enumerate(order):
    v = PU[PU.centromere == c]
    ax.scatter(k + np.random.default_rng(k).uniform(-0.25, 0.25, len(v)), v.nearest_unit_same_centromere * 100, s=22,
               color=LC.TYPE1 if "C_mac_" + c in T1 else LC.TYPE2, edgecolor="white", linewidths=0.3)
ax.set_xticks(range(len(order))); ax.set_xticklabels(order, fontsize=10); ax.set_ylim(-3, 103)
ax.set_ylabel("replicates in which the unit's closest unit\nis from the same centromere (%)", fontsize=10)
ax.set_title("B  Stability of each unit (each dot = one unit)", loc="left", fontsize=13, fontweight="bold")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.suptitle("Bootstrap stability of the unit tree (15-mers, column resampling, NJ)", fontsize=14, fontweight="bold", x=0.02, ha="left")
fig.tight_layout()
for ext in ("png", "pdf"): fig.savefig(W / "figure" / f"bootstrap.{ext}", dpi=450, bbox_inches="tight")
print(CS.to_string(index=False)); print()
print(PU.groupby("centromere").nearest_unit_same_centromere.agg(["min", "median", lambda x: (x < 0.7).sum()]).rename(columns={"<lambda_0>": "units_below_70pct"}).to_string())
