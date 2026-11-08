"""Extract dispersed copies with both canonical ends (q ~5,550 and ~7,820) + 150 bp flanks, in consensus orientation."""
import pathlib, pandas as pd, random
W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22")
class Fa:
    def __init__(self, p):
        self.f = open(p, "rb"); self.i = {}
        for l in open(str(p) + ".fai"):
            n, ln, off, lb, lby = l.split("\t")[:5]; self.i[n] = (int(ln), int(off), int(lb), int(lby))
    def get(self, c, a, b):
        ln, off, lb, lby = self.i[c]; a = max(0, int(a)); b = min(ln, int(b)); s0 = off + (a // lb) * lby + a % lb; s1 = off + (b // lb) * lby + b % lb
        self.f.seek(s0); return self.f.read(s1 - s0).replace(b"\n", b"").decode().upper()
rc = lambda s: s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
G = {"Gar": "F:/data/fClaHyb_func_annotation_v1.0/fClaHyb_Gar/fClaHyb_Gar_genome/fClaHyb_Gar_LG.fa", "Mac": "F:/data/fClaHyb_func_annotation_v1.0/fClaHyb_Mac/fClaHyb_Mac_genome/fClaHyb_Mac_LG.fa"}
random.seed(1); out = open(W / "seq" / "dispersed_fulllength.fa", "w"); n = {}
for sp in ("Gar", "Mac"):
    fa = Fa(G[sp]); E = pd.read_csv(W / "core" / f"{sp}.copy_ends.tsv", sep="\t")
    eh = pd.read_csv(W / "copies" / f"{sp}.hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
    full = E[(E.qmin.between(5540, 5600)) & (E.qmax.between(7790, 7850)) & (~E.cen)]
    rows = full.sample(min(40, len(full)), random_state=1) if sp == "Gar" else full
    k = 0
    for r in rows.itertuples():
        # genomic coordinates of q 5550 and q 7820 from the hits of this copy
        h = eh[(eh.chrom == r.chrom) & (eh[["s0", "s1"]].max(1) > r.lo - 6000) & (eh[["s0", "s1"]].min(1) < r.lo + 8000)]
        a = h[h[["q0", "q1"]].min(1) == r.qmin].iloc[0]; b = h[h[["q0", "q1"]].max(1) == r.qmax].iloc[0]
        fwd = a.s1 > a.s0
        ga = a.s0 if fwd else a.s0; gb = b.s1
        lo, hi = sorted((ga, gb))
        if not 1800 < hi - lo < 3500: continue
        s = fa.get(r.chrom, lo - 151, hi + 150); s = s if fwd else rc(s)
        out.write(f">{sp}|{r.chrom}|{lo}|{hi}|q{r.qmin}-{r.qmax}\n{s}\n"); k += 1
    n[sp] = (len(full), k)
print("full-length dispersed copies (available, written):", n)
