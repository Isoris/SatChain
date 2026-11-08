"""Composition of the spacer flanks: bp covered by other satellite families, by other TE library entries (by class), by the EDTA unit consensus."""
import pathlib, pandas as pd, numpy as np, re
O = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/split")
cols = ["q", "s", "q0", "q1", "s0", "s1", "pid", "len", "e", "bits"]
fl = {}; k = None
for l in open(O / "flanks.fa"):
    l = l.strip()
    if l.startswith(">"): k = l[1:]; fl[k] = 0
    else: fl[k] += len(l)
def cov(df, key=None):
    out = {}
    for q, d in df.groupby("q"):
        m = np.zeros(fl[q], bool)
        for a, b in zip(d.q0, d.q1): m[a - 1:b] = True
        out[q] = m
    return out
sat = pd.read_csv(O / "fl_vs_sat.tsv", sep="\t", header=None, names=cols); sat["fam"] = sat.s.str.split().str[0]
te = pd.read_csv(O / "fl_vs_te.tsv", sep="\t", header=None, names=cols); te["cls"] = te.s.str.split("#").str[1].fillna("?")
un = pd.read_csv(O / "fl_vs_edta_unit.tsv", sep="\t", header=None, names=cols)
rows = []
for q, L in fl.items():
    ch, a, side, ln = q.split("|")
    r = dict(chrom=ch, side=side, bp=L)
    for name, df in (("MacSat03", sat[sat.fam == "MacSat03"]), ("other_sat", sat[sat.fam != "MacSat03"]), ("TE_lib", te), ("EDTA_unit", un)):
        d = df[df.q == q]; m = np.zeros(L, bool)
        for x, y in zip(d.q0, d.q1): m[x - 1:y] = True
        r[name] = m.mean()
    d = te[te.q == q]
    r["top_TE"] = d.groupby("cls").len.sum().idxmax() if len(d) else ""
    d = sat[(sat.q == q) & (sat.fam != "MacSat03")]; r["top_sat"] = d.groupby("fam").len.sum().idxmax() if len(d) else ""
    rows.append(r)
R = pd.DataFrame(rows); pd.set_option("display.width", 250)
g = R.groupby(["chrom", "side"]).agg(n=("bp", "size"), bp=("bp", "median"), MacSat03=("MacSat03", "median"), other_sat=("other_sat", "median"),
                                      TE_lib=("TE_lib", "median"), EDTA_unit=("EDTA_unit", "median"),
                                      top_TE=("top_TE", lambda x: x.value_counts().index[0] if (x != "").any() else ""), top_sat=("top_sat", lambda x: x.value_counts().index[0]))
print(g.round(2).to_string()); R.to_csv(O / "flank_composition.tsv", sep="\t", index=False)
print("\nTE classes over all flanks (bp):", te.groupby("cls").len.sum().sort_values(ascending=False).head(8).to_dict())
print("other satellite families over all flanks (bp):", sat[sat.fam != "MacSat03"].groupby("fam").len.sum().sort_values(ascending=False).head(8).to_dict())
