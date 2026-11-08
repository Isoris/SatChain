"""Step 34a (2026-10-24): compare the step-33 quick split with SatChain v4.6 (runs of record + interruptions)."""
import pathlib, pandas as pd, numpy as np
D = pathlib.Path("F:/data/fClaHyb_func_annotation_v1.0/fClaHyb_Mac/fClaHyb_Mac_satelliteDNA"); O = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/split")
C = pd.read_csv(D / "fClaHyb_Mac.satDNA_v46.catalogue.tsv", sep="\t"); I = pd.read_csv(D / "fClaHyb_Mac.satDNA_v46.interruptions_between_runs.tsv", sep="\t")
REG = pd.read_csv(O / "regions.tsv", sep="\t"); V = pd.read_csv(O / "monomer_validation.tsv", sep="\t"); S = pd.read_csv(O / "spacers.tsv", sep="\t")
V["valid"] = V.score >= 0.30
def inreg(df, a="start0", b="end"):
    m = np.zeros(len(df), bool)
    for r in REG.itertuples(): m |= ((df.chrom == r.chrom) & (df[a] < r.end) & (df[b] > r.start)).values
    return df[m]
c = inreg(C); print("v4.6 runs in the 19 regions by family/name:"); print(c.groupby(["family", "name_v46", "lineage"]).bp.agg(["size", "sum"]).sort_values("sum", ascending=False).head(8).to_string())
print("unit_bp of MacSat03 runs:", c[c.family == "MacSat03"].unit_bp.describe()[["min", "50%", "max"]].to_dict())
def frac_in(lo, hi, runs):
    out = []
    for ch, a, b in zip(lo.index.map(lambda i: V.chrom[i]), lo, hi):
        r = runs[(runs.chrom == ch) & (runs.start0 < b) & (runs.end > a)]
        out.append(sum(min(b, y) - max(a, x) for x, y in zip(r.start0, r.end)) / (b - a))
    return np.array(out)
r03 = c[c.family == "MacSat03"]
for lab, sel in (("valid (my check)", V.valid), ("false (my check)", ~V.valid)):
    v = V[sel]; f = frac_in(v.lo, v.hi, r03); print(f"{lab}: {len(v)} monomer calls, mean fraction inside v4.6 MacSat03 runs {f.mean():.2f}, share >= 0.5 inside {np.mean(f >= .5):.2f}")
i = inreg(I); i = i[i.chrom.isin(c.chrom.unique())]
print("\nv4.6 interruptions in the regions:", len(i), "| class:", i["class"].value_counts().to_dict(), "| bp median", int(i.bp.median()))
# reciprocal match of my TE spacers to v4.6 interruptions
t = S[S.src_bp >= 300]; hit = 0; ovs = []
for r in t.itertuples():
    ch = "C_mac_" + r.chrom; x = i[(i.chrom == ch) & (i.start0 < r.spacer_hi) & (i.end > r.spacer_lo)]
    if len(x):
        o = max(min(r.spacer_hi, b) - max(r.spacer_lo, a) for a, b in zip(x.start0, x.end)); ovs.append(o / max(r.spacer_bp, 1)); hit += 1
print(f"my TE spacers: {len(t)}; overlapping a v4.6 interruption: {hit}; median overlap fraction {np.median(ovs):.2f}")
