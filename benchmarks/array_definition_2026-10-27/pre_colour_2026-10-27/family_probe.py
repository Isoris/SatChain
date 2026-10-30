"""One satellite family for primer / probe design (generalised from garsat04.py): where it is in both genomes, and every figure about it (2026-10-27). Windows Python only.
k-mer set = canonical 21-mers of all catalogue runs of family GarSat04 (Gar), minus k-mers inside microsatellites. Both genomes scanned in 10-kb
windows; signal = share of a window's 21-mers that are in the set. Locus = window with share >= t (t swept 0.1 / 0.3 / 0.5; 0.3 drawn).
Gar windows outside GarSat04 runs and any Mac window above t = sites where a GarSat04-derived probe could also bind (cross-hybridisation in the hybrid).
Outputs (GarSat04/): ideogram_{FAM}.png (both genomes), garsat04_arrays.tsv, garsat04_windows.tsv (windows >= 0.1), {FAM}_figures.tsv (index of
existing figures about GarSat04) + {FAM}_figures_sheet.png (contact sheet), summary printed.
usage: python -P family_probe.py <species Gar|Mac> <family_v46> <out dir>"""
import os, re, glob, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.image as mpimg
from matplotlib.patches import Rectangle
SPF, FAM, B = sys.argv[1], sys.argv[2], sys.argv[3]; os.makedirs(B, exist_ok=True); BM = "F:/data/SatChain/benchmarks/array_definition_2026-10-27"; sys.path.insert(0, BM); R0 = "F:/data/fClaHyb_func_annotation_v1.0"
from kmer_homogeneity import Fasta, kcodes
W = 10000; BS = chr(92)
C = pd.read_csv(f"{R0}/fClaHyb_{SPF}/fClaHyb_{SPF}_satelliteDNA/fClaHyb_{SPF}.satDNA_v46.catalogue.tsv", sep="\t"); X = C[C.family_v46 == FAM].sort_values(["chrom", "start0"])
FA = {sp: Fasta(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.fa") for sp in ("Gar", "Mac")}
def mask_micro(s):
    t = s.decode().upper(); b = bytearray(t.encode())
    for k in range(1, 7):
        for m in re.finditer(("(([ACGT]{%d})" + BS + "2{%d,})") % (k, max(2, 12 // k)), t):
            if m.end() - m.start() >= 16: b[m.start():m.end()] = b"N" * (m.end() - m.start())
    return bytes(b)
ks = np.unique(np.concatenate([np.unique(kcodes(mask_micro(FA[SPF].get(r.chrom, int(r.start0), int(r.end))), 21)) for r in X.itertuples()]))
print(f"{SPF} {FAM}: {len(X)} runs, {X.bp.sum() / 1e6:.2f} Mb, {X.chrom.nunique()} chromosomes, {len(ks):,} distinct 21-mers", flush=True)
rows = []
for sp in ("Gar", "Mac"):
    fa = FA[sp]
    for c in sorted(fa.idx, key=lambda c: int(c.split("LG")[-1])):
        L = fa.idx[c][0]
        for off in range(0, L, 5_000_000):
            seq = fa.get(c, off, min(L, off + 5_000_000 + 20)); codes = kcodes(seq, 21)
            if not len(codes): continue
            hit = np.isin(codes, ks).astype(np.int32); n = len(hit) // W
            if n == 0: continue
            sh = hit[:n * W].reshape(n, W).mean(1)
            for i in np.where(sh >= 0.1)[0]: rows.append(dict(species=sp, chrom=c, start=off + i * W, share=round(float(sh[i]), 3)))
    print("scanned", sp, flush=True)
D = pd.DataFrame(rows)
runs = {c: list(zip(g.start0, g.end)) for c, g in X.groupby("chrom")}
D["in_family_run"] = [sp == SPF and any(a < s + W and b > s for a, b in runs.get(c, [])) for sp, c, s in zip(D.species, D.chrom, D.start)]
D.to_csv(f"{B}/{FAM}_windows.tsv", sep="\t", index=False)
for t in (0.1, 0.3, 0.5):
    d = D[D.share >= t]
    o = "Mac" if SPF == "Gar" else "Gar"
    print(f"t = {t}: {SPF} windows {int((d.species == SPF).sum())} (outside {FAM} runs {int(((d.species == SPF) & ~d.in_family_run).sum())}, chromosomes "
          f"{d[(d.species == SPF) & ~d.in_family_run].chrom.str.split('_').str[-1].unique().tolist()[:12]}); {o} windows {int((d.species == o).sum())} on "
          f"{d[d.species == o].chrom.str.split('_').str[-1].unique().tolist()[:15]}")
# array summary
H = pd.read_csv(f"{BM}/array_levels_{SPF}.tsv", sep="\t").set_index("array_id"); ST = pd.read_csv("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/figure/diagnostics/centromeres/array_strand.tsv", sep="\t").set_index("array_id")
S = X.groupby("array_id").agg(chrom=("chrom", "first"), runs=("run_id", "size"), sat_kb=("bp", lambda v: round(v.sum() / 1e3, 1)), start=("start0", "min"), end=("end", "max"),
                               unit_bp=("unit_bp", lambda v: int(v.mode().iloc[0])), centromeric=("putative_centromeric_array_L1", lambda v: (v == "yes").any())).reset_index()
S["rlvl1"] = S.array_id.map(H.rlvl1); S["rlvl2"] = S.array_id.map(H.rlvl2); S["kmer_recurrence"] = S.array_id.map(H.kmer_recurrence); S["strand_vs_family"] = S.array_id.map(ST.strand_vs_family)
S.sort_values("sat_kb", ascending=False).to_csv(f"{B}/{FAM}_arrays.tsv", sep="\t", index=False)
pd.set_option("display.width", 220); print(S.sort_values("sat_kb", ascending=False).head(20).to_string(index=False))
# ideogram, both genomes side by side
fig, axs = plt.subplots(1, 2, figsize=(22, 13)); fig.patch.set_facecolor("white")
cm = plt.get_cmap("magma_r")
for ax, sp in zip(axs, ("Gar", "Mac")):
    fa = FA[sp]; chroms = sorted(fa.idx, key=lambda c: int(c.split("LG")[-1])); n = len(chroms)
    for i, c in enumerate(chroms):
        y = n - 1 - i; ax.add_patch(Rectangle((0, y - 0.15), fa.idx[c][0] / 1e6, 0.3, color="#eceae4", ec="#9a9890", lw=0.6))
        for r in D[(D.species == sp) & (D.chrom == c)].itertuples():
            ax.add_patch(Rectangle((r.start / 1e6, y - 0.32), W / 1e6, 0.64, color=cm(min(1, r.share / 0.8)), lw=0))
        if sp == SPF:
            for a, b in runs.get(c, []): ax.add_patch(Rectangle((a / 1e6, y + 0.36), max((b - a) / 1e6, 0.05), 0.12, color="#2a78d6", lw=0))
        ax.text(-0.5, y, c.split("_")[-1], ha="right", va="center", fontsize=9)
    ax.set_xlim(-2.5, max(v[0] for v in fa.idx.values()) / 1e6 + 1); ax.set_ylim(-0.8, n - 0.2); ax.set_yticks([]); ax.set_xlabel("Mb")
    for s_ in ("left", "right", "top"): ax.spines[s_].set_visible(False)
    tot = D[(D.species == sp) & (D.share >= 0.3)]
    ax.set_title(f"{sp}: 10-kb windows with {FAM} 21-mers (colour = share, from 0.1; >= 0.3 in {len(tot)} windows)" + (f"; blue = catalogue {FAM} runs" if sp == SPF else
                 " - possible cross-hybridisation in the hybrid"), loc="left", fontsize=11, fontweight="bold")
sm = plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(0, 0.8)); fig.colorbar(sm, ax=axs, fraction=0.015, label=f"share of the window's 21-mers that are {FAM} k-mers")
fig.suptitle(f"{FAM} ({X.name_v46.mode().iloc[0]}) for probe design: {len(X)} runs, {X.bp.sum() / 1e6:.2f} Mb on {X.chrom.nunique()} {SPF} chromosomes; k-mer set = all its runs, microsatellites removed",
             x=0.01, ha="left", fontsize=13, fontweight="bold")
fig.savefig(f"{B}/ideogram_{FAM}.png", dpi=100)
# index of existing figures about GarSat04
ids = set(S.array_id); figs = []
for p in glob.glob("F:/data/**/*.png", recursive=True) if False else []: pass
roots = ["F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/figure/diagnostics/centromeres", "F:/data/satellite_FISH/motif_views_v2", BM]
for root in roots:
    for p in glob.glob(os.path.join(root, "**", "*.png"), recursive=True):
        q = p.replace("\\", "/"); b = os.path.basename(q)
        if any(a in b for a in ids) or FAM in b or b.startswith(("genome_divmap", f"ideogram_{SPF}_rlvl", f"centromeres_rlvl1_overview_{SPF}", "kmer_homogeneity")):
            figs.append(q)
F = pd.DataFrame({"figure": sorted(set(figs))}); F.to_csv(f"{B}/{FAM}_figures.tsv", sep="\t", index=False); print(len(F), f"figures about {FAM} indexed")
n = len(F); cols = 6; rr = (n + cols - 1) // cols; fig, axs = plt.subplots(rr, cols, figsize=(cols * 4.2, rr * 3.4)); fig.patch.set_facecolor("white")
for ax in np.ravel(axs): ax.axis("off")
for i, q in enumerate(F.figure):
    ax = np.ravel(axs)[i]; ax.imshow(mpimg.imread(q)); ax.set_title(f"{i + 1}. {os.path.basename(q)[:46]}", fontsize=7)
fig.tight_layout(); fig.savefig(f"{B}/{FAM}_figures_sheet.png", dpi=75)