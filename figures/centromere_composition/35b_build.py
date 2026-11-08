"""Step 35b (2026-10-24): data for the enlarged-centromere composition figure (both haplotypes).
Region per chromosome = largest cluster of SatChain v4.6 putative centromeric runs (runs <= 1 Mb apart) +-150 kb.
Tracks: v4.6 satellite runs (family, strand by reverse_fraction), EDTA TE features with < 50% overlap to satellite runs (insertions),
assembly gaps, homogeneity profile h(x) = Jaccard of canonical k-mer sets of the windows [x-W, x) and [x, x+W) (step 5 kb;
W = 10/20/50 kb, k = 15/21/31). Homogenized core = maximal stretch with h >= tau (tau = 0.3/0.5/0.7), gaps <= 10 kb merged, >= 50 kb.
Outputs in F:/data/satellite_FISH/centromere_composition_2026-10-24/."""
import pathlib, gzip, collections, numpy as np, pandas as pd, sys
exec(open(pathlib.Path(__file__).parent / "32e_extract.py").read().split("random.seed")[0])
OUT = pathlib.Path("F:/data/satellite_FISH/centromere_composition_2026-10-24")
STEP = 5000; WS = (10000, 20000, 50000); KS = (15, 21, 31); TAUS = (0.3, 0.5, 0.7)
RC = str.maketrans("ACGT", "TGCA")
def kmer_hashes(s, k):
    r = s.translate(RC)[::-1]; n = len(s); out = np.empty(max(0, n - k + 1), np.int64); bad = np.zeros(len(out), bool)
    for i in range(n - k + 1):
        a = s[i:i + k]
        if "N" in a: bad[i] = True; continue
        b = r[n - i - k:n - i]; out[i] = hash(a if a < b else b)
    return out, bad
regs, sats, tes, homs, cores, gaps = [], [], [], [], [], []
for sp in ("Mac", "Gar"):
    base = pathlib.Path(f"F:/data/fClaHyb_func_annotation_v1.0/fClaHyb_{sp}")
    C = pd.read_csv(base / "_satelliteDNA" / f"fClaHyb_{sp}.satDNA_v46.catalogue.tsv", sep="\t"); c = C[C.putative_centromeric_array_L1 == "yes"]
    fa = Fa(G[sp]); sizes = dict(l.split()[:2] for l in open(base / "_genome" / f"fClaHyb_{sp}_LG.chrom.sizes"))
    te = pd.read_csv(base / "_repeats_EDTA" / f"fClaHyb_{sp}.TEanno.gff3", sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 6, 8], names=["chrom", "type", "s", "e", "strand", "attr"])
    te = te[~te.type.isin(["repeat_region", "target_site_duplication", "long_terminal_repeat"])]
    gp = pd.read_csv(base / "_genome" / f"fClaHyb_{sp}_LG.gaps.detgaps.bed", sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "s", "e"])
    for ch, d in c.groupby("chrom"):
        d = d.sort_values("start0"); cl = []
        for r in d.itertuples():
            if cl and r.start0 - cl[-1][1] <= 1_000_000: cl[-1][1] = max(cl[-1][1], r.end); cl[-1][2] += r.bp
            else: cl.append([r.start0, r.end, r.bp])
        b = max(cl, key=lambda x: x[2]); lo = max(0, b[0] - 150000); hi = min(int(sizes[ch]), b[1] + 150000)
        regs.append(dict(sp=sp, chrom=ch, lo=lo, hi=hi, cen_lo=b[0], cen_hi=b[1]))
        sr = C[(C.chrom == ch) & (C.end > lo) & (C.start0 < hi)]
        for r in sr.itertuples(): sats.append(dict(sp=sp, chrom=ch, s=r.start0, e=r.end, strand="-" if r.reverse_fraction > 0.5 else "+", family=r.name_v46, lineage=r.lineage, element=r.element_class, gene=r.gene_flag, identity=r.mean_identity_conf, unit=r.unit_bp))
        satmask = np.zeros(hi - lo, bool)
        for r in sr.itertuples(): satmask[max(0, r.start0 - lo):max(0, min(hi, r.end) - lo)] = True
        t = te[(te.chrom == ch) & (te.e > lo) & (te.s < hi)]
        for r in t.itertuples():
            a, z = max(lo, r.s - 1), min(hi, r.e)
            if z <= a or satmask[a - lo:z - lo].mean() >= 0.5: continue
            cls = r.type; nm = [x.split("=")[1] for x in r.attr.split(";") if x.startswith("Name=")]
            tes.append(dict(sp=sp, chrom=ch, s=a, e=z, strand=r.strand, type=cls, name=nm[0] if nm else ""))
        for r in gp[(gp.chrom == ch) & (gp.e > lo) & (gp.s < hi)].itertuples(): gaps.append(dict(sp=sp, chrom=ch, s=r.s, e=r.e))
        seq = fa.get(ch, lo, hi)
        for k in KS:
            hsh, bad = kmer_hashes(seq, k)
            for W in WS:
                for x in range(lo + W, hi - W + 1, STEP):
                    i = x - lo; A = hsh[i - W:i - k + 1][~bad[i - W:i - k + 1]]; B = hsh[i:i + W - k + 1][~bad[i:i + W - k + 1]]
                    if len(A) < W / 2 or len(B) < W / 2: h = np.nan
                    else:
                        A = np.unique(A); B = np.unique(B); inter = len(np.intersect1d(A, B, assume_unique=True)); h = inter / (len(A) + len(B) - inter)
                    homs.append((sp, ch, x, W, k, h))
        print(sp, ch, "done", flush=True)
H = pd.DataFrame(homs, columns=["sp", "chrom", "pos", "W", "k", "h"])
for (sp, ch, W, k), d in H.groupby(["sp", "chrom", "W", "k"]):
    d = d.sort_values("pos")
    for tau in TAUS:
        seg = []
        for p, h in zip(d.pos, d.h):
            if h >= tau:
                if seg and p - seg[-1][1] <= 10000 + STEP: seg[-1][1] = p
                else: seg.append([p, p])
        for a, z in seg:
            if z - a >= 50000: cores.append(dict(sp=sp, chrom=ch, W=W, k=k, tau=tau, s=a - STEP // 2, e=z + STEP // 2))
for nm, rows in (("regions", regs), ("satellite_runs", sats), ("te_insertions", tes), ("gaps", gaps), ("cores", cores)): pd.DataFrame(rows).to_csv(OUT / f"{nm}.tsv", sep="\t", index=False)
H.to_csv(OUT / "homogeneity.tsv.gz", sep="\t", index=False)
print("ok")
