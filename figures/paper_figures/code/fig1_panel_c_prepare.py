"""Fig. 1c data cache (run once; ~2 min, < 400 MB): strand blocks of the monomer copies inside v4.6 satellite runs, per centromere region,
exactly the rule of reference_figures/35c_figure.py (copies with identity >= 60, copies fully inside a run, merged when same strand and < 2 kb apart,
majority strand = strand with more bp). Output panels/_data/fig1c_strand_blocks.tsv (sp, chrom, s, e, major)."""
import pandas as pd, numpy as np, pathlib
O = pathlib.Path("F:/data/satellite_FISH/centromere_composition_2026-10-24"); OUT = pathlib.Path(__file__).resolve().parent.parent / "panels" / "_data" / "fig1c_strand_blocks.tsv"
R = pd.read_csv(O / "regions.tsv", sep="\t"); S = pd.read_csv(O / "satellite_runs.tsv", sep="\t"); rows = []
for sp in ("Mac", "Gar"):
    keep = []
    rr = R[R.sp == sp].set_index("chrom")
    for c in pd.read_csv(f"F:/data/satellite_FISH/pipeline_v4_5_fixed/out/fClaHyb_{sp}_satellites/03_monomer_copies.tsv", sep="\t",
                         usecols=["chrom", "genomic_start_0based", "genomic_end", "strand", "best_alignment_identity"], chunksize=500_000):
        c = c[(c.best_alignment_identity >= 60) & c.chrom.isin(rr.index)]
        for ch, d in c.groupby("chrom"):
            d = d[(d.genomic_end > rr.lo[ch]) & (d.genomic_start_0based < rr.hi[ch])]
            if len(d): keep.append(d)
    c = pd.concat(keep)
    for ch, lo, hi in zip(rr.index, rr.lo, rr.hi):
        s = S[(S.sp == sp) & (S.chrom == ch)]; d = c[(c.chrom == ch)].sort_values("genomic_start_0based")
        ins = np.zeros(len(d), bool)
        for r in s.itertuples(): ins |= ((d.genomic_start_0based >= r.s) & (d.genomic_end <= r.e)).values
        d = d[ins]; out = []
        for r in d.itertuples():
            if out and out[-1][2] == r.strand and r.genomic_start_0based - out[-1][1] < 2000: out[-1][1] = r.genomic_end
            else: out.append([r.genomic_start_0based, r.genomic_end, r.strand])
        plus = sum(b - a for a, b, t in out if t == "+"); minus = sum(b - a for a, b, t in out if t == "-"); maj = "+" if plus >= minus else "-"
        rows += [(sp, ch, a, b, t == maj) for a, b, t in out]
    print(sp, "done", flush=True)
pd.DataFrame(rows, columns=["sp", "chrom", "s", "e", "major"]).to_csv(OUT, sep="\t", index=False)
