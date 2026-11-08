"""Sensitivity of the part-by-part divergence results (2026-10-23; user: stability under sampling regimes and k-mer sizes).
Grid: k in {9,11,13,15,17,21} x window in {250,500,1000} bp x pairs per comparison in {100,300,1000} x regime {all units, 13 units per centromere}
x 5 seeds. Claims, evaluated in every run:
 C1 homogenization: median within-centromere / median between-centromere (same type), pooled over element windows < 0.2
 C2 RT core (5,560-7,400): dispersed copies more diverged than centromeric copies of different centromeres
 C3 RT core = lowest between-centromere divergence among element regions (core vs non-core windows, satellite stretch excluded)
 C4 insert (between type II centromeres) < RT core between centromeres
 C5 built-in satellite stretch (2,929-3,388) = highest between-centromere divergence among element regions
Outputs: kmer_tree/sensitivity_runs.tsv, kmer_tree/sensitivity_summary.tsv, figure/sensitivity.png/.pdf"""
import pathlib, itertools, random, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/code"); import ltr_colours as LC

W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22"); A = pathlib.Path("F:/data/fClaHyb_func_annotation_v1.0")
R = pathlib.Path("F:/data/satellite_FISH/runs_of_record_v4.5.3_2026-10-18/Mac")
T1 = ["C_mac_LG07", "C_mac_LG04", "C_mac_LG01", "C_mac_LG23", "C_mac_LG02"]; T2 = ["C_mac_LG08", "C_mac_LG25", "C_mac_LG16"]
KS, WINS, NPS, REGIMES, SEEDS = [9, 11, 13, 15, 17, 21], [250, 500, 1000], [100, 300, 1000], ["all units", "13 per centromere"], range(5)
CORE, SATS = (5560, 7400), (2929, 3388)
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
def kset(s, K): return {min(s[i:i + K], s[i:i + K].translate(comp)[::-1]) for i in range(len(s) - K + 1) if "N" not in s[i:i + K]}
def mash(a, b, K):
    if len(a) < 20 or len(b) < 20: return np.nan
    j = len(a & b) / len(a | b); return 1.0 if j == 0 else max(0.0, -1 / K * np.log(2 * j / (1 + j)))

U = pd.read_csv(W / "unit" / "units.tsv", sep="\t"); U = U[U.chrom.isin(T1 + T2)].reset_index(drop=True)
eh = pd.read_csv(W / "copies" / "Mac.hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
eh["lo"] = eh[["s0", "s1"]].min(1); eh["hi"] = eh[["s0", "s1"]].max(1)
ih = pd.read_csv(W / "insert" / "Mac.insert_hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
ih["lo"] = ih[["s0", "s1"]].min(1); ih["hi"] = ih[["s0", "s1"]].max(1)
cp = pd.read_csv(W / "core" / "Mac.core_copies.tsv", sep="\t"); arm = cp[~cp.in_cen_1Mb]

def window_seq(h, w0, w1):
    parts = []
    for r in h.itertuples():
        a, b = max(r.q0, w0), min(r.q1, w1)
        if b - a < 40: continue
        f = lambda q: r.s0 + (q - r.q0) * (r.s1 - r.s0) / max(1, r.q1 - r.q0)
        x0, x1 = sorted((int(f(a)), int(f(b)))); parts.append(fa.get(r.chrom, x0, x1))
    return "".join(parts)

# sequences per item, per window size (cached once)
items = []
for r in U.itertuples():
    a, b = int(r.start), int(r.start + r.period)
    hh = eh[(eh.chrom == r.chrom) & (eh.hi > a) & (eh.lo < b)]; hi = ih[(ih.chrom == r.chrom) & (ih.hi > a) & (ih.lo < b) & (ih.len >= 300)]
    items.append(dict(group=r.chrom, type="I" if r.chrom in T1 else "II", kind="cen", h=hh,
                      insert="".join(fa.get(r.chrom, max(s, a), min(e, b)) for s, e in zip(hi.lo, hi.hi))))
for r in arm.itertuples():
    items.append(dict(group=f"arm_{r.chrom}_{r.lo}", type="-", kind="arm", insert="",
                      h=eh[(eh.chrom == r.chrom) & (eh.hi > r.lo - 8000) & (eh.lo < r.hi + 8000)]))
seqs = {ws: [{w0: window_seq(it["h"], w0, w0 + ws) for w0 in range(0, 8164, ws)} for it in items] for ws in WINS}
print("sequences cached", flush=True)

runs = []
for K in KS:
    for ws in WINS:
        wins = list(range(0, 8164, ws))
        S = [{w0: kset(s, K) for w0, s in d.items()} for d in seqs[ws]]
        Si = [kset(it["insert"], K) for it in items]
        for regime in REGIMES:
            for NPAIR in NPS:
                for seed in SEEDS:
                    rnd = random.Random(seed)
                    cen = [i for i, x in enumerate(items) if x["kind"] == "cen"]
                    if regime == "13 per centromere":
                        cen = [i for c in T1 + T2 for i in rnd.sample([j for j in cen if items[j]["group"] == c], 13)]
                    armi = [i for i, x in enumerate(items) if x["kind"] == "arm"]
                    within = [(i, j) for i, j in itertools.combinations(cen, 2) if items[i]["group"] == items[j]["group"]]
                    btw = [(i, j) for i, j in itertools.combinations(cen, 2) if items[i]["group"] != items[j]["group"] and items[i]["type"] == items[j]["type"]]
                    bt2 = [(i, j) for i, j in btw if items[i]["type"] == "II"]
                    disp = list(itertools.combinations(armi, 2))
                    smp = lambda P: rnd.sample(P, min(NPAIR, len(P)))
                    within, btw, bt2, disp = smp(within), smp(btw), smp(bt2), smp(disp)
                    med = lambda P, w0: np.nanmedian([mash(S[i][w0], S[j][w0], K) for i, j in P]) if P else np.nan
                    prof_w = {w0: med(within, w0) for w0 in wins}; prof_b = {w0: med(btw, w0) for w0 in wins}
                    core = [w0 for w0 in wins if w0 >= CORE[0] - ws / 2 and w0 + ws <= CORE[1] + ws / 2]
                    sat = [w0 for w0 in wins if w0 < SATS[1] and w0 + ws > SATS[0]]
                    other = [w0 for w0 in wins if w0 not in core and w0 not in sat]
                    b_core = np.nanmedian([prof_b[w] for w in core]); b_sat = np.nanmax([prof_b[w] for w in sat]); b_other = np.nanmedian([prof_b[w] for w in other])
                    d_core = np.nanmedian([med(disp, w) for w in core])
                    ins = np.nanmedian([mash(Si[i], Si[j], K) for i, j in bt2])
                    ratio = np.nanmedian(list(prof_w.values())) / np.nanmedian(list(prof_b.values()))
                    runs.append(dict(k=K, window=ws, pairs=NPAIR, regime=regime, seed=seed, within_over_between=ratio, between_core=b_core,
                                     between_other=b_other, between_satstretch=b_sat, dispersed_core=d_core, insert_between_typeII=ins,
                                     C1=ratio < 0.2, C2=d_core > b_core, C3=b_core < b_other, C4=ins < b_core, C5=b_sat > b_other))
        print("k", K, "window", ws, "done", flush=True)
RN = pd.DataFrame(runs); RN.to_csv(W / "kmer_tree" / "sensitivity_runs.tsv", sep="\t", index=False)
CL = {"C1": "homogenization: within/between < 0.2", "C2": "RT core: dispersed > centromeric (between centromeres)",
      "C3": "RT core most conserved element region", "C4": "insert least diverged (< RT core)", "C5": "built-in satellite stretch most diverged"}
SUM = pd.DataFrame([dict(claim=c, description=d, runs=len(RN), share_true=RN[c].mean(),
                         by_k="; ".join(f"k{k}:{RN[RN.k == k][c].mean():.2f}" for k in KS),
                         by_window="; ".join(f"{w}:{RN[RN.window == w][c].mean():.2f}" for w in WINS),
                         by_regime="; ".join(f"{g}:{RN[RN.regime == g][c].mean():.2f}" for g in REGIMES)) for c, d in CL.items()])
SUM.to_csv(W / "kmer_tree" / "sensitivity_summary.tsv", sep="\t", index=False)

fig, axs = plt.subplots(2, 3, figsize=(17, 9.5), facecolor="white"); axs = axs.ravel()
metrics = [("within_over_between", "within / between (homogenization)", 0.2), ("between_core", "RT core, between centromeres", None),
           ("dispersed_core", "RT core, dispersed copies", None), ("insert_between_typeII", "insert, between type II centromeres", None),
           ("between_satstretch", "built-in satellite stretch, between centromeres", None)]
for ax, (m, t, thr) in zip(axs, metrics):
    for ws, col in zip(WINS, [LC.REL[3], LC.REL[2], LC.REL[0]]):
        g = RN[RN.window == ws].groupby("k")[m]; ax.errorbar(KS, g.median().values, yerr=[g.median() - g.quantile(0.05), g.quantile(0.95) - g.median()],
                                                              fmt="-o", color=col, label=f"{ws}-bp windows", capsize=3, ms=4)
    if thr: ax.axhline(thr, color=LC.INK, ls="--", lw=1)
    ax.set_xlabel("k-mer size"); ax.set_title(t, fontsize=11, fontweight="bold", loc="left")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axs[0].legend(fontsize=8, frameon=False)
ax = axs[5]; ax.axis("off")
txt = "Claims held (share of all runs)\n\n" + "\n".join(f"{c}  {r.share_true * 100:5.1f}%  {CL[c]}" for c, r in zip(CL, SUM.itertuples()))
txt += f"\n\n{len(RN)} runs: k {KS}, windows {WINS} bp,\npairs {NPS}, regimes: all units / 13 per centromere, 5 seeds.\nBars: 5th–95th percentile over pairs, regimes, seeds."
ax.text(0, 1, txt, va="top", fontsize=10, family="monospace")
fig.suptitle("Sensitivity of the divergence results to k-mer size, window size and sampling", fontsize=14, fontweight="bold", x=0.02, ha="left")
fig.tight_layout()
for ext in ("png", "pdf"): fig.savefig(W / "figure" / f"sensitivity.{ext}", dpi=450, bbox_inches="tight")
print(SUM[["claim", "share_true", "by_k", "by_window", "by_regime"]].to_string(index=False))
