"""Truncation / deletion state of the source RTE in every centromeric unit (source coordinates 1-3,117; LINE 3' end = 1, 5' end = 3,117)."""
import pathlib, pandas as pd, numpy as np
W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22")
B = pd.read_csv(W / "homologous" / "src_vs_units.tsv", sep="\t", header=None, names=["unit", "q0", "q1", "s0", "s1", "pid", "len", "bits"])
L = 3117; rows = []; prof = {}
for u, d in B.groupby("unit"):
    cov = np.zeros(L + 1, bool)
    # best non-overlapping hits in query space
    taken = []
    for r in d.sort_values("bits", ascending=False).itertuples():
        if any(r.q0 < b - 10 and r.q1 > a + 10 for a, b in taken): continue
        taken.append((r.q0, r.q1)); cov[r.q0:r.q1 + 1] = True
    prof[u] = cov[1:]
    bins = cov[1:].reshape(-1, 1) if False else None
    present = np.where(cov)[0]
    rows.append(dict(unit=u, chrom=u.split("|")[0], covered=int(cov.sum()), q3end=int(present.min()), q5end=int(present.max()),
                     has_3tail=present.min() <= 30, has_5end=present.max() >= L - 30,
                     gap_1700_2680=cov[1800:2600].mean() < 0.2, mid_present=round(cov[1750:2650].mean(), 2)))
R = pd.DataFrame(rows)
g = R.groupby("chrom").agg(units=("unit", "size"), covered_median=("covered", "median"), q3end=("q3end", "median"), q5end=("q5end", "median"),
                           tail_pct=("has_3tail", "mean"), full5_pct=("has_5end", "mean"), mid_gap_pct=("gap_1700_2680", "mean"))
T2 = {"LG08", "LG25", "LG16"}
g["type"] = ["II" if c in T2 else "I" for c in g.index]
pd.set_option("display.width", 200); print(g.round(2).sort_values(["type", "covered_median"]).to_string())
# 100-bp presence profile per chromosome
P = pd.DataFrame({u: v[:3100].reshape(31, 100).mean(1) > 0.5 for u, v in prof.items()}).T
P["chrom"] = [u.split("|")[0] for u in P.index]
pp = P.groupby("chrom").mean()
print("\n100-bp source bins present (>=50% of units), 3' end left -> 5' end right; # present, . absent")
for c, r in pp.iterrows(): print(f"  {c:5s} {'II' if c in T2 else 'I ':2s} " + "".join("#" if v >= 0.5 else "." for v in r.values))
R.to_csv(W / "homologous" / "source_state_per_unit.tsv", sep="\t", index=False); pp.to_csv(W / "homologous" / "source_bins_per_chrom.tsv", sep="\t")
