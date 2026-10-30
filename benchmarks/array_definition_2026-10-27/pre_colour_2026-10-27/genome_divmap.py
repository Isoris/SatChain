"""Genome-wide pairwise identity of the major centromeres (window x window), Mac and Gar together (2026-10-27).
Regions = the 54 major centromeres (centromeres/index.tsv: first run start - 20 kb .. last run end + 20 kb). Windows of W bp (5 kb); canonical
15-mers (strand-free: half the arrays are '-' vs their family); bottom-s MinHash sketch (s = 1000) per window. Shared sketch hashes between two
windows counted through an inverted index (hash -> windows); J ~ shared / s; identity (Mash) = (2J / (1 + J))^(1/k).
Outputs: genome_divmap_windows.npz (identity matrix + window table), genome_divmap_centromeres.tsv (per pair of centromeres: 95th percentile and
max of window-pair identity, share of windows of A with a partner >= 0.9 in B), genome_divmap.png (window map, centromere matrix).
usage: python -P genome_divmap.py [W=5000]"""
import os, sys, collections, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
B = os.path.dirname(os.path.abspath(__file__)); R0 = "F:/data/fClaHyb_func_annotation_v1.0"; D0 = "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/figure/diagnostics/centromeres"
W = int(sys.argv[1]) if len(sys.argv) > 1 else 5000; K = 15; SK = 1000; FL = 20000
sys.path.insert(0, B); from kmer_homogeneity import Fasta, kcodes
I = pd.read_csv(f"{D0}/index.tsv", sep="\t"); win = []; sk = []
for sp in ("Mac", "Gar"):
    C = pd.read_csv(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_satelliteDNA/fClaHyb_{sp}.satDNA_v46.catalogue.tsv", sep="\t"); fa = Fasta(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.fa")
    for r in I[I.species == sp].itertuples():
        x = C[C.array_id == r.array_id]; lo = max(0, int(x.start0.min()) - FL); hi = min(fa.idx[r.chrom][0], int(x.end.max()) + FL)
        seq = fa.get(r.chrom, lo, hi)
        for i in range(0, len(seq) - W + 1, W):
            c = np.unique(kcodes(seq[i:i + W], K)); h = (c * np.uint64(0x9E3779B97F4A7C15)) ^ (c >> np.uint64(29)); h = np.sort(h)[:SK]
            win.append(dict(species=sp, chrom=r.chrom, cen=f"{sp} {r.chrom.split('_')[-1]}", family=r.family, start=lo + i, n=len(h)))
            sk.append(h)
T = pd.DataFrame(win); n = len(T); print(n, "windows", flush=True)
inv = collections.defaultdict(list)
for wi, h in enumerate(sk):
    for v in h.tolist(): inv[v].append(wi)
S = np.zeros((n, n), dtype=np.uint16)
for lst in inv.values():
    if 1 < len(lst) <= 3000:
        a = np.array(lst); S[np.ix_(a, a)] += 1
J = S / SK; np.fill_diagonal(J, 1.0)
ID = np.where(J > 0, (2 * J / (1 + J)) ** (1 / K), 0.0).astype(np.float32)
low = T.n.values < 200; ID[low, :] = np.nan; ID[:, low] = np.nan                    # windows with few k-mers (N / low complexity)
np.savez_compressed(f"{B}/genome_divmap_windows.npz", ID=ID, start=T.start.values, cen=T.cen.values)
cens = list(dict.fromkeys(T.cen)); idx = {c: np.where(T.cen.values == c)[0] for c in cens}; rows = []
for a in cens:
    for b in cens:
        if a == b: continue
        M = ID[np.ix_(idx[a], idx[b])]
        if np.all(np.isnan(M)): continue
        rows.append(dict(A=a, B=b, fam_A=T.family[idx[a][0]], fam_B=T.family[idx[b][0]], p95=round(float(np.nanpercentile(M, 95)), 3), max=round(float(np.nanmax(M)), 3),
                         share_A_windows_with_partner_ge_0_9=round(float(np.mean(np.nanmax(M, 1) >= 0.9)), 3)))
P = pd.DataFrame(rows); P.to_csv(f"{B}/genome_divmap_centromeres.tsv", sep="\t", index=False)
Mx = P.pivot(index="A", columns="B", values="share_A_windows_with_partner_ge_0_9").reindex(index=cens, columns=cens)
fig = plt.figure(figsize=(20, 9.5)); fig.patch.set_facecolor("white")
ax = fig.add_axes([0.03, 0.08, 0.44, 0.84]); ax.imshow(ID, cmap="magma_r", vmin=0.75, vmax=1.0, interpolation="nearest")
bnd = [idx[c][0] for c in cens] + [n]; [ax.axhline(b, color="#2a78d6", lw=0.3) for b in bnd]; [ax.axvline(b, color="#2a78d6", lw=0.3) for b in bnd]
ax.set_xticks([(idx[c][0] + idx[c][-1]) / 2 for c in cens]); ax.set_xticklabels(cens, rotation=90, fontsize=5.5); ax.set_yticks(ax.get_xticks()); ax.set_yticklabels(cens, fontsize=5.5)
ax.set_title(f"A  window x window identity ({W // 1000}-kb windows, 15-mer MinHash, strand-free) across the 54 major centromeres", loc="left", fontsize=11, fontweight="bold")
ax = fig.add_axes([0.53, 0.08, 0.42, 0.84]); im = ax.imshow(Mx.values, cmap="viridis", vmin=0, vmax=1, interpolation="nearest")
ax.set_xticks(range(len(cens))); ax.set_xticklabels(cens, rotation=90, fontsize=6); ax.set_yticks(range(len(cens))); ax.set_yticklabels([f"{c} {T.family[idx[c][0]]}" for c in cens], fontsize=6)
fig.colorbar(im, ax=ax, fraction=0.04, label="share of row-centromere windows with a >= 0.9-identity partner in the column centromere")
ax.set_title("B  which centromeres carry the same sequence", loc="left", fontsize=11, fontweight="bold")
fig.savefig(f"{B}/genome_divmap.png", dpi=110)
pd.set_option("display.width", 200)
print(P.sort_values("share_A_windows_with_partner_ge_0_9", ascending=False).head(25).to_string(index=False))