"""Anchored, homologous unit boundaries (2026-10-23; user: fix the boundaries so the same segment is compared; compare only what is comparable).
Anchor = a fixed position of the element consensus inside the RT core (default 6,500; 6,000 and 7,000 tested), kept by type I and type II units.
For every Mac element hit covering the anchor (consensus coordinate q0 <= A <= q1), the genome position of the anchor is interpolated inside the hit,
strand-aware; anchors closer than 1 kb are merged. Centromeric anchors (centromere +-1 Mb) are sorted per chromosome; a unit is the interval between
two consecutive anchors on the same strand, 8-25 kb long, oriented so the element reads forward (minus-strand arrays are read from the higher anchor).
Checks: strand consistency per array, period agreement between anchors 6,000 / 6,500 / 7,000, shift of the old (core-start) boundaries.
Outputs: unit_anchored/units_A6500.tsv (+ A6000, A7000), unit_anchored/checks.txt, figure/anchored_units.png/.pdf"""
import pathlib, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/code"); import ltr_colours as LC

W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22"); R = pathlib.Path("F:/data/satellite_FISH/runs_of_record_v4.5.3_2026-10-18/Mac")
OUT = W / "unit_anchored"; OUT.mkdir(exist_ok=True)
T1 = ["C_mac_LG07", "C_mac_LG04", "C_mac_LG01", "C_mac_LG23", "C_mac_LG02"]; T2 = ["C_mac_LG08", "C_mac_LG25", "C_mac_LG16"]
lab = lambda c: "LG" + c.split("LG")[-1]
eh = pd.read_csv(W / "copies" / "Mac.hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
cen = pd.read_csv(R / "centromere_check.tsv", sep="\t").set_index("chrom")

def anchors(A):
    h = eh[(eh.q0 <= A) & (eh.q1 >= A)].copy()
    h["pos"] = (h.s0 + (A - h.q0) * (h.s1 - h.s0) / (h.q1 - h.q0).clip(lower=1)).round().astype(int)
    h["strand"] = np.where(h.s1 >= h.s0, "+", "-")
    h = h[[c in cen.index and cen.loc[c, "cen_start"] - 1e6 < p < cen.loc[c, "cen_end"] + 1e6 for c, p in zip(h.chrom, h.pos)]]
    out = []
    for (c, st), g in h.sort_values("pos").groupby(["chrom", "strand"]):
        last = -1e9
        for r in g.sort_values("bits", ascending=False).sort_values("pos").itertuples():
            if r.pos - last >= 1000: out.append((c, st, r.pos, r.pid)); last = r.pos
    return pd.DataFrame(out, columns=["chrom", "strand", "pos", "pid"])

def units(A):
    an = anchors(A); rows = []
    for c, g in an.groupby("chrom"):
        g = g.sort_values("pos").reset_index(drop=True)
        for i in range(len(g) - 1):
            a, b = g.loc[i], g.loc[i + 1]
            if a.strand != b.strand: continue
            P = b.pos - a.pos
            if not 8000 <= P <= 25000: continue
            start, end = (a.pos, b.pos) if a.strand == "+" else (b.pos, a.pos)     # oriented: unit starts at the anchor and reads with the element
            rows.append(dict(chrom=c, strand=a.strand, lo=a.pos, hi=b.pos, oriented_start=start, oriented_end=end, period=P, anchor=A))
    return an, pd.DataFrame(rows)

res = {A: units(A) for A in (6000, 6500, 7000)}
for A, (an, u) in res.items(): u.to_csv(OUT / f"units_A{A}.tsv", sep="\t", index=False)
an, u = res[6500]
lines = [f"anchor 6,500: {len(an)} centromeric anchors; {len(u)} units (8-25 kb, same strand) on {u.chrom.nunique()} chromosomes"]
st = an.groupby("chrom").strand.agg(lambda s: f"{(s == '+').sum()}+/{(s == '-').sum()}-")
lines.append("anchor strands per chromosome: " + ", ".join(f"{lab(c)} {v}" for c, v in st.items()))
per = u.groupby("chrom").period.agg(["size", "median", "std"]).round(0)
lines.append("units per chromosome (n, median period, sd):\n" + per.to_string())
# period agreement between anchors: match units by nearest lo
agree = []
for A in (6000, 7000):
    u2 = res[A][1]
    for c in u.chrom.unique():
        a1 = u[u.chrom == c].sort_values("lo"); a2 = u2[u2.chrom == c].sort_values("lo")
        if a2.empty: continue
        for r in a1.itertuples():
            j = (a2.lo - r.lo).abs().idxmin(); d = abs(a2.loc[j, "lo"] - r.lo)
            if d < 3000: agree.append((A, abs(a2.loc[j, "period"] - r.period)))
ag = pd.DataFrame(agree, columns=["anchor", "period_diff"])
lines.append("period difference vs anchor 6,500 (bp, median / 90th pct): " +
             "; ".join(f"A{A}: {g.period_diff.median():.0f} / {g.period_diff.quantile(.9):.0f} (n={len(g)})" for A, g in ag.groupby("anchor")))
old = pd.read_csv(W / "unit" / "units.tsv", sep="\t")
sh = []
for r in old.itertuples():
    n = u[u.chrom == r.chrom]
    if len(n): sh.append(np.min(np.abs(n.lo.values - r.start)))
sh = np.array(sh)
lines.append(f"old unit starts (core-copy start) vs anchored boundaries: median shift {np.median(sh):.0f} bp, 90th pct {np.percentile(sh, 90):.0f} bp, "
             f"> 1 kb in {(sh > 1000).mean() * 100:.0f}% of old units")
(OUT / "checks.txt").write_text("\n".join(lines) + "\n", encoding="utf-8"); print("\n".join(lines))

fig, axs = plt.subplots(1, 3, figsize=(17, 5.2), facecolor="white")
order = [c for c in T1 + T2 if c in set(u.chrom)]
ax = axs[0]
for k, c in enumerate(order):
    v = u[u.chrom == c].period / 1000
    ax.scatter(k + np.random.default_rng(k).uniform(-0.2, 0.2, len(v)), v, s=14, color=LC.TYPE1 if c in T1 else LC.TYPE2, edgecolor="white", linewidths=0.3)
ax.set_xticks(range(len(order))); ax.set_xticklabels([lab(c) for c in order]); ax.set_ylabel("unit length (kb), anchored")
ax.set_title("A  Anchored unit lengths per centromere", loc="left", fontweight="bold")
ax = axs[1]
for A, col in ((6000, LC.REL[3]), (7000, LC.REL[0])):
    g = ag[ag.anchor == A].period_diff; ax.hist(g.clip(upper=2000), bins=40, color=col, histtype="stepfilled", label=f"anchor {A:,} vs 6,500")
ax.set_xlabel("difference in unit length (bp; capped at 2 kb)"); ax.set_ylabel("units"); ax.legend(frameon=False, fontsize=9)
ax.set_title("B  Same units with a different anchor", loc="left", fontweight="bold")
ax = axs[2]; ax.hist(np.clip(sh, 0, 8000), bins=40, color=LC.REL[2])
ax.set_xlabel("shift of the old unit start from the anchored boundary (bp)"); ax.set_ylabel("old units")
ax.set_title("C  How far the old boundaries were off", loc="left", fontweight="bold")
for a in axs:
    for s in ("top", "right"): a.spines[s].set_visible(False)
fig.suptitle("Unit boundaries fixed at a homologous anchor (element position 6,500, RT core)", fontsize=14, fontweight="bold", x=0.02, ha="left")
fig.tight_layout()
for ext in ("png", "pdf"): fig.savefig(W / "figure" / f"anchored_units.{ext}", dpi=450, bbox_inches="tight")
