import pathlib, pandas as pd
exec(open("32e_extract.py").read().split("random.seed")[0])
U = pd.read_csv(W / "unit_anchored" / "units_A6500.tsv", sep="\t").sort_values(["chrom", "lo"])
fa = Fa(G["Mac"]); out = open(W / "seq" / "units_sample.fa", "w")
for c, g in U.groupby("chrom"):
    for r in g.iloc[[len(g) // 3, 2 * len(g) // 3]].itertuples():
        s = fa.get(c, r.lo, r.hi); s = s if r.strand == "+" else rc(s)
        out.write(f">{c.replace('C_mac_', '')}|{r.lo}\n{s}\n")
out.close(); print("ok")
