import pathlib, pandas as pd, numpy as np
W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22")
for sp in ("Mac", "Gar"):
    cc = pd.read_csv(W / "core" / f"{sp}.core_copies.tsv", sep="\t")
    eh = pd.read_csv(W / "copies" / f"{sp}.hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
    eh["lo"] = eh[["s0", "s1"]].min(1); eh["hi"] = eh[["s0", "s1"]].max(1)
    eh = eh[~((eh[["q0", "q1"]].min(1) >= 2880) & (eh[["q0", "q1"]].max(1) <= 3400))]
    res = []
    for g, d in eh.groupby("chrom"):
        d = d.sort_values("lo"); c = cc[cc.chrom == g]
        for r in c.itertuples():
            # collinear chain around the core copy: hits within 3 kb gaps, same orientation as core
            x = d[(d.hi > r.lo - 6000) & (d.lo < r.hi + 6000)]
            ori = "+" if r.strand == "+" else "-"
            x = x[(x.s1 > x.s0) == (ori == "+")]
            if x.empty: continue
            res.append(dict(sp=sp, chrom=g, lo=r.lo, cen=r.in_cen_1Mb, qmin=x[["q0", "q1"]].min(1).min(), qmax=x[["q0", "q1"]].max(1).max(),
                            covered=sum(min(b, 8164) - a for a, b in zip(x[["q0", "q1"]].min(1), x[["q0", "q1"]].max(1)))))
    R = pd.DataFrame(res)
    for cen, d in R.groupby("cen"):
        print(f"\n{sp} {'centromeric' if cen else 'arm'} n={len(d)}")
        print("  qmin (low-coordinate end) top bins:", (d.qmin // 50 * 50).value_counts().head(6).to_dict())
        print("  qmax (high-coordinate end) top bins:", (d.qmax // 50 * 50).value_counts().head(6).to_dict())
        print("  span qmax-qmin quantiles:", d.eval("qmax-qmin").quantile([.1, .5, .9, 1]).round().to_dict())
    R.to_csv(W / "core" / f"{sp}.copy_ends.tsv", sep="\t", index=False)
