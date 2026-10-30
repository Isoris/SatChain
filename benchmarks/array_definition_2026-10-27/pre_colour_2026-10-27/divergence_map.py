"""Pairwise divergence map of a chromosome region: where divergence collapses = an array (2026-10-27).
Windows of w bp (default 2 kb) along the region; each window -> canonical k-mer set (k = 15) -> MinHash bottom-s sketch (s = 500).
Jaccard J between every pair of windows from the sketches; identity estimate (Mash) = (2J / (1 + J))^(1/k); divergence = 1 - identity.
Boundaries: insulation score = mean identity in the d x d square straddling each window boundary (d = 5 windows); a boundary is a local
minimum of the insulation score whose drop below the mean of its two flanking maxima is >= delta (0.05). Arrays = stretches between boundaries
with mean within-block identity >= 0.8 (block of homogenised sequence). Off-diagonal blocks = the same sequence elsewhere in the region.
Panels: identity heatmap (upper triangle), catalogue runs and assembly gaps below, insulation score with boundaries.
usage: python -P divergence_map.py <species> <chrom> <start> <end> <png> [w=2000] [k=15]"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
SP, CH, S0, E0, PNG = sys.argv[1], sys.argv[2], int(float(sys.argv[3])), int(float(sys.argv[4])), sys.argv[5]
W = int(sys.argv[6]) if len(sys.argv) > 6 else 2000; K = int(sys.argv[7]) if len(sys.argv) > 7 else 15; SK = 500; D = 5; DELTA = 0.05
R = f"F:/data/fClaHyb_func_annotation_v1.0/fClaHyb_{SP}/"
LUT = np.full(256, 255, dtype=np.uint8)
for i, ch in enumerate(b"ACGT"): LUT[ch] = i; LUT[ch + 32] = i
idx = {l.split()[0]: tuple(map(int, l.split()[1:5])) for l in open(R + f"fClaHyb_{SP}_genome/fClaHyb_{SP}_LG.fa.fai")}
fh = open(R + f"fClaHyb_{SP}_genome/fClaHyb_{SP}_LG.fa", "rb")
def get(c, a, b):
    L, off, lb, lby = idx[c]; fh.seek(off + (a // lb) * lby + a % lb)
    return fh.read((b - a) + ((b - a) // lb + 2) * (lby - lb)).replace(b"\n", b"").replace(b"\r", b"")[:b - a]
def sketch(seq):
    s = LUT[np.frombuffer(seq, dtype=np.uint8)]; n = len(s) - K + 1
    if n <= 0: return np.zeros(0, dtype=np.uint64)
    bad = np.convolve((s == 255).astype(np.int32), np.ones(K, dtype=np.int32), "valid") > 0
    s64 = np.where(s == 255, 0, s).astype(np.uint64); f = np.zeros(n, dtype=np.uint64); r = np.zeros(n, dtype=np.uint64)
    for j in range(K): f = (f << np.uint64(2)) | s64[j:j + n]; r = r | ((np.uint64(3) - s64[j:j + n]) << np.uint64(2 * j))
    c = np.unique(np.minimum(f, r)[~bad]); h = (c * np.uint64(0x9E3779B97F4A7C15)) ^ (c >> np.uint64(29))   # cheap integer hash
    return np.sort(h)[:SK]
seq = get(CH, S0, E0); nw = len(seq) // W
sk = [sketch(seq[i * W:(i + 1) * W]) for i in range(nw)]
I = np.full((nw, nw), np.nan)
for i in range(nw):
    for j in range(i, nw):
        a, b = sk[i], sk[j]
        if len(a) < 50 or len(b) < 50: continue
        u = np.union1d(a, b)[:SK]; J = len(np.intersect1d(np.intersect1d(a, b), u)) / len(u)
        I[i, j] = I[j, i] = (2 * J / (1 + J)) ** (1 / K) if J > 0 else 0.0
ins = np.full(nw, np.nan)
for i in range(D, nw - D): ins[i] = np.nanmean(I[i - D:i, i:i + D])
bnd = []
for i in range(D + 1, nw - D - 1):
    if np.isnan(ins[i]) or ins[i] > np.nanmin(ins[i - 2:i + 3]): continue
    lm, rm = np.nanmax(ins[max(D, i - 3 * D):i]), np.nanmax(ins[i + 1:min(nw - D, i + 3 * D)])
    if (lm + rm) / 2 - ins[i] >= DELTA: bnd.append(i)
edges = [0] + bnd + [nw]; blocks = []
# v2: an array = contiguous windows whose insulation (mean identity across the boundary) stays above T. A composite array (satellite blocks +
# spacers) has a lattice of off-diagonal stripes and a mean identity of only ~0.4, but a clear insulation plateau; the old ">= 0.8 block" rule
# missed it. Dips shorter than D windows inside a plateau are bridged. T swept 0.10 / 0.15 / 0.20.
ARR = {}
for T in (0.10, 0.15, 0.20):
    on = np.nan_to_num(ins, nan=0) >= T; segs = []; i = 0
    while i < nw:
        if on[i]:
            j = i
            while j + 1 < nw and (on[j + 1] or (np.any(on[j + 1:j + 1 + D]) and not np.all(~on[j + 1:j + 1 + D]))): j += 1
            if j - i + 1 >= 3: segs.append((i, j + 1))
            i = j + 1
        else: i += 1
    ARR[T] = segs
blocks = [(a, b, 1.0) for a, b in ARR[0.15]]                                 # reported default T = 0.15
CAT = pd.read_csv(R + f"fClaHyb_{SP}_satelliteDNA/fClaHyb_{SP}.satDNA_v46.catalogue.tsv", sep="\t")
runs = CAT[(CAT.chrom == CH) & (CAT.end > S0) & (CAT.start0 < E0)]
GP = [int(l.split()[1]) for l in open(R + f"fClaHyb_{SP}_genome/fClaHyb_{SP}_LG.gaps.detgaps.bed") if l.split()[0] == CH and S0 < int(l.split()[1]) < E0]
fig = plt.figure(figsize=(14, 15)); fig.patch.set_facecolor("white")
ax = fig.add_axes([0.08, 0.36, 0.84, 0.58]); x = lambda w_: (S0 + w_ * W) / 1e6
im = ax.imshow(np.where(np.triu(np.ones_like(I), 0) > 0, I, np.nan), cmap="magma_r", vmin=0.7, vmax=1.0, extent=[x(0), x(nw), x(nw), x(0)], interpolation="nearest")
for a, b, m in blocks:
    if not np.isnan(m) and m >= 0.8: ax.add_patch(Rectangle((x(a), x(a)), x(b) - x(a), x(b) - x(a), fill=False, ec="#2a78d6", lw=1.5))
ax.set_ylabel("Mb"); ax.set_title(f"{SP} {CH.split('_')[-1]} {S0 / 1e6:.2f}-{E0 / 1e6:.2f} Mb: pairwise identity of {W // 1000}-kb windows (k = {K} MinHash); "
                                   f"blue = arrays (insulation plateau >= 0.15)", loc="left", fontsize=10.5, fontweight="bold")
fig.colorbar(im, ax=ax, fraction=0.03, label="identity (Mash estimate)")
ax2 = fig.add_axes([0.08, 0.27, 0.735, 0.06], sharex=ax)
for r in runs.itertuples(): ax2.add_patch(Rectangle((r.start0 / 1e6, 0.55), (r.end - r.start0) / 1e6, 0.35, color="#c0392b", lw=0))
for g in GP: ax2.axvline(g / 1e6, color="black", lw=2)
ax2.set_ylim(0, 1); ax2.set_yticks([0.72]); ax2.set_yticklabels(["catalogue runs"], fontsize=8)
for a, b, m in blocks:
    if not np.isnan(m) and m >= 0.8: ax2.add_patch(Rectangle((x(a), 0.08), x(b) - x(a), 0.35, color="#2a78d6", lw=0))
ax2.text(x(0), 0.25, "divergence-collapse blocks ", ha="right", fontsize=8, color="#2a78d6")
ax3 = fig.add_axes([0.08, 0.06, 0.735, 0.18], sharex=ax)
ax3.plot([x(i) for i in range(nw)], ins, color="#444")
for T, ls in ((0.10, ":"), (0.15, "-"), (0.20, "--")): ax3.axhline(T, color="#2a78d6", lw=0.8, ls=ls)
ax3.set_ylabel("insulation (mean identity across)"); ax3.set_xlabel(f"{CH} (Mb)"); ax3.grid(alpha=0.3)
fig.savefig(PNG, dpi=105)
sat_runs = sum(min(r.end, E0) - max(r.start0, S0) for r in runs.itertuples()) / 1e3; blk = sum((b - a) * W for a, b, m in blocks if not np.isnan(m) and m >= 0.8) / 1e3
print(f"wrote {PNG}: {nw} windows; catalogue runs {sat_runs:.0f} kb ({runs.start0.min() / 1e6:.3f}-{runs.end.max() / 1e6:.3f} Mb)" if len(runs) else f"wrote {PNG}")
for T, segs in ARR.items(): print(f"   T = {T:.2f}: arrays " + ", ".join(f"{x(a):.3f}-{x(b):.3f} Mb ({(b - a) * W / 1e3:.0f} kb)" for a, b in segs if (b - a) * W >= 10000))