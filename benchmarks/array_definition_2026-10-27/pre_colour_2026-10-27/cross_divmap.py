"""Cross divergence map between two arrays: conserved chain vs rearrangement (2026-10-27).
Every W-bp window of array A (x) against every window of array B (y): MinHash identity of canonical 15-mers (strand-free) plus the strand of the
best match (share of shared k-mers that are identical in forward orientation; '+' blue, '-' red). Both arrays shown in their family orientation
(array_strand.tsv: a '-' array is reversed so both read in the same direction). A conserved chain = a diagonal run of high identity; a shift,
break, duplication or inversion = diagonal pieces offset / repeated / anti-diagonal; no signal = unrelated. Side tracks: runs (satellite) and
RTE hits of each array. Also the best partner of every window of A in B (identity) along A.
usage: python -P cross_divmap.py <array A> <array B> <png> [W=2000]"""
import sys, csv, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
sys.path.insert(0, "F:/data/SatChain/benchmarks/array_definition_2026-10-27"); from kmer_homogeneity import Fasta
AA, AB, PNG = sys.argv[1:4]; W = int(sys.argv[4]) if len(sys.argv) > 4 else 2000; K = 15; SK = 600; R0 = "F:/data/fClaHyb_func_annotation_v1.0"
D0 = "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/figure/diagnostics/centromeres"
STR = {r["array_id"]: r["strand_vs_family"] for r in csv.DictReader(open(f"{D0}/array_strand.tsv"), delimiter="\t")}
RTE = pd.read_csv("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/split/rte_genome_hits_both.tsv", sep="\t"); RTE["a"] = RTE[["sstart", "send"]].min(1); RTE["b"] = RTE[["sstart", "send"]].max(1)
LUT = np.full(256, 255, dtype=np.uint8)
for i, ch in enumerate(b"ACGT"): LUT[ch] = i; LUT[ch + 32] = i
def codes(seq):
    s = LUT[np.frombuffer(seq, dtype=np.uint8)]; n = len(s) - K + 1
    if n <= 0: return np.zeros(0, np.uint64), np.zeros(0, np.uint64)
    bad = np.convolve((s == 255).astype(np.int32), np.ones(K, dtype=np.int32), "valid") > 0; s64 = np.where(s == 255, 0, s).astype(np.uint64)
    f = np.zeros(n, np.uint64); r = np.zeros(n, np.uint64)
    for j in range(K): f = (f << np.uint64(2)) | s64[j:j + n]; r = r | ((np.uint64(3) - s64[j:j + n]) << np.uint64(2 * j))
    return f[~bad], r[~bad]
def load(aid):
    sp = aid.split("_")[0]; C = pd.read_csv(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_satelliteDNA/fClaHyb_{sp}.satDNA_v46.catalogue.tsv", sep="\t")
    x = C[C.array_id == aid].sort_values("start0"); fa = Fasta(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.fa"); lo, hi = int(x.start0.min()), int(x.end.max())
    seq = fa.get(x.chrom.iloc[0], lo, hi).upper(); rev = STR.get(aid) == "-"
    runs = [((hi - b) if rev else (a - lo), (hi - a) if rev else (b - lo)) for a, b in zip(x.start0, x.end)]
    h = RTE[(RTE.chrom == x.chrom.iloc[0]) & (RTE.b > lo) & (RTE.a < hi)]; rte = [((hi - b) if rev else (a - lo), (hi - a) if rev else (b - lo)) for a, b in zip(h.a, h.b)]
    if rev: seq = seq.translate(bytes.maketrans(b"ACGTN", b"TGCAN"))[::-1]
    wins = []
    for i in range(0, len(seq) - W + 1, W):
        f, r = codes(seq[i:i + W]); c = np.minimum(f, r); u, ix = np.unique(c, return_index=True)
        hsh = (u * np.uint64(0x9E3779B97F4A7C15)) ^ (u >> np.uint64(29)); o = np.argsort(hsh)[:SK]
        wins.append((set(hsh[o].tolist()), dict(zip(hsh[o].tolist(), (f[ix][o] == u[o]).tolist()))))
    return x.chrom.iloc[0], lo, hi, wins, runs, rte, rev
ca, la, ha, wa, ra, ta, va = load(AA); cb, lb, hb, wb, rb, tb, vb = load(AB)
M = np.zeros((len(wb), len(wa))); SG = np.zeros_like(M)
for i, (sa, oa) in enumerate(wa):
    for j, (sb, ob) in enumerate(wb):
        sh = sa & sb
        if not sh: continue
        J = len(sh) / SK; M[j, i] = (2 * J / (1 + J)) ** (1 / K)
        SG[j, i] = np.mean([oa[h] == ob[h] for h in sh])                      # same orientation of the k-mer in both windows
fig = plt.figure(figsize=(14, 15)); fig.patch.set_facecolor("white")
ax = fig.add_axes([0.12, 0.25, 0.72, 0.68]); X = len(wa) * W / 1e3; Y = len(wb) * W / 1e3
rgb = np.ones((len(wb), len(wa), 3)); s = np.clip((M - 0.75) / 0.25, 0, 1)
rgb[..., 0] = 1 - s * (SG < 0.5) * 0.15 - s * (SG >= 0.5) * 0.85; rgb[..., 1] = 1 - s * 0.75; rgb[..., 2] = 1 - s * (SG < 0.5) * 0.85 - s * (SG >= 0.5) * 0.15
ax.imshow(rgb, extent=[0, X, Y, 0], aspect="auto", interpolation="nearest")
ax.set_xlabel(f"{AA} ({ca.split('_')[-1]}, kb, family orientation{' - reversed' if va else ''})"); ax.set_ylabel(f"{AB} ({cb.split('_')[-1]}, kb{' - reversed' if vb else ''})")
ax.set_title(f"{AA} ({ca.split('_')[-1]}) x {AB} ({cb.split('_')[-1]}): {W // 1000}-kb window identity (darker = higher, from 0.75); blue = same orientation, red = inverted\n"
             "diagonal = conserved chain; offset / repeated / anti-diagonal pieces = rearrangement", loc="left", fontsize=11, fontweight="bold")
axt = fig.add_axes([0.12, 0.935, 0.72, 0.02], sharex=ax); axt.axis("off")
for a, b in ra: axt.add_patch(Rectangle((a / 1e3, 0), (b - a) / 1e3, 1, color="#2c3e50", lw=0))
for a, b in ta: axt.add_patch(Rectangle((a / 1e3, 1.05), (b - a) / 1e3, 0.6, color="#f2c14e", lw=0))
axt.set_ylim(0, 1.7)
axr = fig.add_axes([0.845, 0.25, 0.02, 0.68], sharey=ax); axr.axis("off")
for a, b in rb: axr.add_patch(Rectangle((0, a / 1e3), 1, (b - a) / 1e3, color="#2c3e50", lw=0))
for a, b in tb: axr.add_patch(Rectangle((1.05, a / 1e3), 0.6, (b - a) / 1e3, color="#f2c14e", lw=0))
axr.set_xlim(0, 1.7)
ax2 = fig.add_axes([0.12, 0.06, 0.72, 0.13], sharex=ax)
best = M.max(0); ax2.plot((np.arange(len(wa)) + 0.5) * W / 1e3, best, color="#444"); ax2.set_ylim(0.6, 1.01); ax2.grid(alpha=0.3)
ax2.set_ylabel("best identity in B"); ax2.set_xlabel(f"{AA} (kb)"); ax2.set_title(f"every window of {AA}: identity of its best partner in {AB} (low = sequence of A missing from B)", loc="left", fontsize=10)
fig.savefig(PNG, dpi=105)
print(f"wrote {PNG}: {len(wa)} x {len(wb)} windows; windows of A with a partner >= 0.9: {np.mean(best >= 0.9):.2f}, >= 0.8: {np.mean(best >= 0.8):.2f}; "
      f"inverted share among matches >= 0.85: {np.mean(SG[M >= 0.85] < 0.5) if (M >= 0.85).any() else float('nan'):.2f}")