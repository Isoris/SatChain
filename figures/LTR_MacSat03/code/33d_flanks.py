"""Step 33d: extract the non-satellite, non-RTE flanks of each TE spacer (array orientation) for identification."""
import pathlib, pandas as pd
exec(open(pathlib.Path(__file__).parent / "32e_extract.py").read().split("random.seed")[0])
O = W / "split"; fa = Fa(G["Mac"]); S = pd.read_csv(O / "spacers.tsv", sep="\t"); S = S[S.src_bp >= 300]
out = open(O / "flanks.fa", "w"); n = 0
for r in S.itertuples():
    c = "C_mac_" + r.chrom
    L = (r.spacer_lo, r.te_lo); R = (r.te_hi, r.spacer_hi)
    if r.orient == "-": L, R = R, L
    for side, (a, b) in (("left", L), ("right", R)):
        if b - a < 200: continue
        s = fa.get(c, a, b); s = s if r.orient == "+" else rc(s)
        out.write(f">{r.chrom}|{a}|{side}|{b-a}\n{s}\n"); n += 1
print(n, "flanks")
