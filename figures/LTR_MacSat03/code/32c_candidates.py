import pathlib, pandas as pd, numpy as np
W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22")
S = pd.read_csv(W / "tree" / "spine_arm_copies.tsv", sep="\t")
cc = pd.read_csv(W / "core" / "Mac.core_copies.tsv", sep="\t")
full = pd.read_csv(W / "copies" / "Mac.copies.tsv", sep="\t")
mono = pd.read_csv("F:/data/satellite_FISH/pipeline_v4_5_fixed/out/fClaHyb_Mac_satellites/36_monomers_final.tsv", sep="\t", usecols=["chrom", "genomic_start_0based", "genomic_end", "family_id"])
eh = pd.read_csv(W / "copies" / "Mac.hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"]); eh["lo"] = eh[["s0", "s1"]].min(1); eh["hi"] = eh[["s0", "s1"]].max(1)
rows = []
for r in S.itertuples():
    sp, lg, p, w, s = r.copy.split("|"); ch = f"C_mac_{lg}"; p = int(p)
    c = cc[(cc.chrom == ch) & (cc.lo == p)].iloc[0]
    f = full[(full.chrom == ch) & (full.lo <= c.hi) & (full.hi >= c.lo)]
    fl = f.iloc[0] if len(f) else None
    # element-consensus coverage within +-10 kb, excluding built-in satellite-only hits
    h = eh[(eh.chrom == ch) & (eh.hi > c.lo - 10000) & (eh.lo < c.hi + 10000)]
    h = h[~((h[["q0", "q1"]].min(1) >= 2880) & (h[["q0", "q1"]].max(1) <= 3400))]
    cov = np.zeros(8164, bool)
    for x in h.itertuples(): cov[min(x.q0, x.q1) - 1:max(x.q0, x.q1)] = True
    m = mono[(mono.chrom == ch) & (mono.genomic_end > c.lo - 20000) & (mono.genomic_start_0based < c.hi + 20000)]
    fams = m.family_id.value_counts()
    rows.append(dict(level=r.spine_level, copy=r.copy, core_pid=round(c.pid, 1), full_copy_bp=(fl.hi - fl.lo) if fl is not None else np.nan,
                     full_qrange=f"{fl.qlo}-{fl.qhi}" if fl is not None else "", element_cov_10kb=round(cov.mean(), 2),
                     monomers_20kb=len(m), sat_families=";".join(f"{k}:{v}" for k, v in fams.head(3).items())))
D = pd.DataFrame(rows); pd.set_option("display.width", 250); print(D.to_string(index=False))
D.to_csv(W / "tree" / "spine_arm_context.tsv", sep="\t", index=False)
