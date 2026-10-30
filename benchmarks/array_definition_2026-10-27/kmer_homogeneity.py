"""k-mer homogeneity of run-based arrays (2026-10-27).
Arrays = run-based (B, G = 20 kb: neighbouring runs of the same family, never across an assembly gap; compare_array_methods.py), >= 20 kb of
satellite. Sequence = the runs only (interruptions excluded). Canonical k-mers (k = 15 / 21 / 31; strand-free), counted with numpy 2-bit codes.
Per array and k:
  top_unit_cov = share of all k-mer occurrences carried by the L most frequent k-mers, L = unit_bp (catalogue). A perfectly homogeneous
                 array repeats the same L k-mers -> ~1. Variants, divergence and HOR structure lower it.
  n90_per_unit = distinct k-mers needed to cover 90 % of occurrences / L. ~1 when homogeneous; grows with the number of distinct variants.
Outputs: kmer_homogeneity_<sp>.tsv, kmer_homogeneity.png (ranking and scatter), printed top 15 per species.
usage: python -P kmer_homogeneity.py"""
import os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "F:/data/satellite_FISH/runs_of_record_v4.6_2026-10-27/code"); import palette_v46 as PAL
sys.path.insert(0, "F:/data/SatChain/satchain"); from gaps import Gaps
B = os.path.dirname(os.path.abspath(__file__)); R = "F:/data/fClaHyb_func_annotation_v1.0"
LUT = np.full(256, 255, dtype=np.uint8)
for i, ch in enumerate(b"ACGT"): LUT[ch] = i; LUT[ch + 32] = i
class Fasta:
    def __init__(self, p):
        self.fh = open(p, "rb"); self.idx = {l.split()[0]: tuple(map(int, l.split()[1:5])) for l in open(p + ".fai")}
    def get(self, c, a, b):
        L, off, lb, lby = self.idx[c]; self.fh.seek(off + (a // lb) * lby + a % lb)
        return self.fh.read((b - a) + ((b - a) // lb + 2) * (lby - lb)).replace(b"\n", b"").replace(b"\r", b"")[:b - a]
def kcodes(seq, k):
    s = LUT[np.frombuffer(seq, dtype=np.uint8)]; n = len(s) - k + 1
    if n <= 0: return np.zeros(0, dtype=np.uint64)
    bad = np.convolve((s == 255).astype(np.int32), np.ones(k, dtype=np.int32), "valid") > 0
    f = np.zeros(n, dtype=np.uint64); r = np.zeros(n, dtype=np.uint64); s64 = np.where(s == 255, 0, s).astype(np.uint64)
    for j in range(k):
        f = (f << np.uint64(2)) | s64[j:j + n]; r = r | ((np.uint64(3) - s64[j:j + n]) << np.uint64(2 * j))
    return np.minimum(f, r)[~bad]
def plot_scatter(D):
    fig, axs = plt.subplots(1, 2, figsize=(17, 6.5)); fig.patch.set_facecolor("white")
    for ax, sp in zip(axs, ("Mac", "Gar")):
        d = D[D.species == sp]
        for cen, col, lab in ((False, "#d0ccbf", "not centromeric"), (True, PAL.CEN_EDGE, "centromeric (L1)")):
            x = d[d.centromeric == cen]; ax.scatter(x.sat_kb, x.top_unit_cov_k21, s=16, color=col, edgecolor=PAL.dark(PAL.CEN_EDGE, 0.6) if cen else "white", linewidths=0.5, label=lab)
        for r in d.sort_values("top_unit_cov_k21", ascending=False).head(8).itertuples(): ax.annotate(f"{r.family} {r.chrom.split('_')[-1]}", (r.sat_kb, r.top_unit_cov_k21), fontsize=7.5)
        ax.set_xscale("log"); ax.set_xlabel("satellite kb (run-based array)"); ax.set_ylabel("top-unit k-mer coverage (k = 21)"); ax.grid(alpha=0.3); ax.legend(fontsize=9, frameon=False)
        ax.set_title(f"{sp}: k-mer homogeneity of run-based arrays (n = {len(d)}); 1 = the same unit repeated", loc="left", fontsize=11, fontweight="bold")
    fig.tight_layout(); fig.savefig(f"{B}/kmer_homogeneity.png", dpi=450); fig.savefig(f"{B}/kmer_homogeneity.png".replace(".png", ".pdf"))

if __name__ == "__main__":
    if os.environ.get("PLOT_ONLY"):                              # redraw from the saved tables (colour changes only)
        plot_scatter(pd.concat([pd.read_csv(f"{B}/kmer_homogeneity_{sp}.tsv", sep="\t") for sp in ("Mac", "Gar")])); sys.exit()
    rows = []
    for sp in ("Mac", "Gar"):
        C = pd.read_csv(f"{R}/fClaHyb_{sp}/fClaHyb_{sp}_satelliteDNA/fClaHyb_{sp}.satDNA_v46.catalogue.tsv", sep="\t").sort_values(["chrom", "start0"]).reset_index(drop=True)
        GP = Gaps(f"{R}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.gaps.detgaps.bed"); fa = Fasta(f"{R}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.fa")
        bid, cur, prev = [], 0, None
        for r in C.itertuples():
            if not (prev is not None and prev.chrom == r.chrom and r.family_v46 == prev.family_v46 and 0 <= r.start0 - prev.end <= 20000 and not GP.between(r.chrom, prev.end, r.start0)): cur += 1
            bid.append(cur); prev = r
        C["B"] = bid
        for b, g in C.groupby("B"):
            if g.bp.sum() < 20000: continue
            seq = b"N".join(fa.get(r.chrom, int(r.start0), int(r.end)) for r in g.itertuples()); L = max(1, int(g.unit_bp.mode().iloc[0]))
            rec = dict(species=sp, array=f"{sp}_B{b:05d}", step02_ids=",".join(sorted(g.array_id.unique()))[:60], chrom=g.chrom.iloc[0], start=int(g.start0.min()), end=int(g.end.max()),
                       runs=len(g), sat_kb=round(g.bp.sum() / 1e3, 1), unit_bp=L, family=g.family_v46.iloc[0], name=str(g.name_v46.iloc[0]),
                       centromeric=(g.putative_centromeric_array_L1 == "yes").any(), lineage=g.lineage.iloc[0])
            for k in (15, 21, 31):
                _, cnt = np.unique(kcodes(seq, k), return_counts=True); cnt = np.sort(cnt)[::-1]; tot = cnt.sum(); cs = np.cumsum(cnt)
                rec[f"top_unit_cov_k{k}"] = round(float(cs[min(L, len(cs)) - 1] / tot), 3); rec[f"n90_per_unit_k{k}"] = round(float((np.searchsorted(cs, 0.9 * tot) + 1) / L), 2)
            rows.append(rec)
    D = pd.DataFrame(rows)
    for sp in ("Mac", "Gar"): D[D.species == sp].sort_values("top_unit_cov_k21", ascending=False).to_csv(f"{B}/kmer_homogeneity_{sp}.tsv", sep="\t", index=False)
    pd.set_option("display.width", 260); pd.set_option("display.max_colwidth", 40)
    cols = ["array", "chrom", "start", "runs", "sat_kb", "unit_bp", "family", "name", "centromeric", "top_unit_cov_k15", "top_unit_cov_k21", "top_unit_cov_k31", "n90_per_unit_k21"]
    for sp in ("Mac", "Gar"):
        d = D[D.species == sp].sort_values("top_unit_cov_k21", ascending=False); print(f"== {sp}: {len(d)} run-based arrays >= 20 kb; most homogeneous (k = 21):"); print(d[cols].head(15).to_string(index=False)); print()
        rk = d[[f"top_unit_cov_k{k}" for k in (15, 21, 31)]].rank(ascending=False); print("   rank agreement across k (Spearman k15-k21, k21-k31):",
              round(rk.iloc[:, 0].corr(rk.iloc[:, 1]), 2), round(rk.iloc[:, 1].corr(rk.iloc[:, 2]), 2)); print()
    plot_scatter(D)
