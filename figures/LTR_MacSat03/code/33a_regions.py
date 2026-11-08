"""Step 33a (2026-10-24): array regions (from anchored units, <= 5 kb gaps, +-15 kb), written as FASTA for TE mapping.
Coordinates in names: chrom:start (0-based)."""
import pathlib, pandas as pd
exec(open(pathlib.Path(__file__).parent / "32e_extract.py").read().split("random.seed")[0])
O = W / "split"; O.mkdir(exist_ok=True)
U = pd.read_csv(W / "unit_anchored" / "units_A6500.tsv", sep="\t").sort_values(["chrom", "lo"])
reg = []
for c, g in U.groupby("chrom"):
    s = e = None
    for r in g.itertuples():
        if s is None: s, e = r.lo, r.hi
        elif r.lo - e <= 5000: e = max(e, r.hi)
        else: reg.append((c, s, e)); s, e = r.lo, r.hi
    reg.append((c, s, e))
fa = Fa(G["Mac"]); out = open(O / "regions.fa", "w"); rows = []
for c, s, e in reg:
    a, b = s - 15000, e + 15000; out.write(f">{c}:{a}\n{fa.get(c, a, b)}\n"); rows.append((c, a, b))
pd.DataFrame(rows, columns=["chrom", "start", "end"]).to_csv(O / "regions.tsv", sep="\t", index=False)
print(len(rows), "regions", sum(b - a for _, a, b in rows) / 1e6, "Mb")
