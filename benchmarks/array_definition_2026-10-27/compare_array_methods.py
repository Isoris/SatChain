"""Compare two definitions of a satellite array on the v4.6 catalogue (2026-10-27).
A (current, inherited): array = v4.5 step-02 candidate interval (array_id). Step 02: 10-kb windows, 5-kb step, base identity at lag 40-1000 bp
  >= 0.34, positive windows merged when <= 20 kb apart, edges refined to the TRF same-period union. Every run inside the interval belongs to it.
B (SatChain, from the runs): array = chain of neighbouring catalogue runs (56 + val56, sequence-checked) of the SAME family_v46, interruption
  <= G bp, never across an assembly gap (detgaps BED, zero-length joins included). G swept 10 / 20 / 50 kb.
Per G: correspondence of A and B arrays (1:1, A split by B, A merged by B, complex), why A is split (mixed family / assembly gap / interruption
  > G), arrays and bp, and the major centromere per chromosome (largest centromeric array by run bp; gene arrays allowed) under A and under B.
Outputs: compare_<G>.tsv per array, summary.tsv, centromeres.tsv.
usage: python -P compare_array_methods.py"""
import os, sys, collections, pandas as pd, numpy as np
sys.path.insert(0, "F:/data/SatChain/satchain"); from gaps import Gaps
R = "F:/data/fClaHyb_func_annotation_v1.0"; O = os.path.dirname(os.path.abspath(__file__))
SUM, CEN = [], []
for sp in ("Mac", "Gar"):
    C = pd.read_csv(f"{R}/fClaHyb_{sp}/fClaHyb_{sp}_satelliteDNA/fClaHyb_{sp}.satDNA_v46.catalogue.tsv", sep="\t").sort_values(["chrom", "start0"]).reset_index(drop=True)
    GP = Gaps(f"{R}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.gaps.detgaps.bed")
    A = C.groupby("array_id").agg(chrom=("chrom", "first"), runs=("run_id", "size"), bp=("bp", "sum"), lo=("start0", "min"), hi=("end", "max"),
                                  n_fam=("family_v46", "nunique"), cen=("putative_centromeric_array_L1", lambda v: (v == "yes").any()))
    A["gap_inside"] = [any(GP.between(c, int(a), int(b)) for _ in [0]) for c, a, b in zip(A.chrom, A.lo, A.hi)]
    for G in (10000, 20000, 50000):
        bid, cur, why = [], 0, collections.Counter()
        prev = None
        for r in C.itertuples():
            new = True
            if prev is not None and prev.chrom == r.chrom:
                d = r.start0 - prev.end; same = r.family_v46 == prev.family_v46; gap = GP.between(r.chrom, prev.end, r.start0)
                new = not (same and 0 <= d <= G and not gap)
                if new and prev.array_id == r.array_id:                     # a break inside an A array: record why (first reason that applies)
                    why["family change" if not same else "assembly gap" if gap else f"interruption > {G // 1000} kb"] += 1
            if new: cur += 1
            bid.append(f"{sp}_B{G // 1000}k_{cur:05d}"); prev = r
        C[f"B{G}"] = bid
        AB = C.groupby("array_id")[f"B{G}"].nunique(); BA = C.groupby(f"B{G}").array_id.nunique()
        kind = {}
        for a, nb in AB.items():
            bs = C.loc[C.array_id == a, f"B{G}"].unique(); merged = any(BA[b] > 1 for b in bs)
            kind[a] = "1:1" if nb == 1 and not merged else ("A split by B" if nb > 1 and not merged else ("A merged by B" if nb == 1 else "complex"))
        Ak = pd.Series(kind); B = C.groupby(f"B{G}").agg(chrom=("chrom", "first"), bp=("bp", "sum"), runs=("run_id", "size"), fam=("family_v46", "first"),
                                                         cen=("putative_centromeric_array_L1", lambda v: (v == "yes").any()))
        A_ = A.assign(kind=Ak, n_B=AB); A_.to_csv(os.path.join(O, f"compare_{sp}_{G // 1000}k.tsv"), sep="\t")
        multi = A_[A_.runs >= 2]
        SUM.append(dict(species=sp, G_kb=G // 1000, A_arrays=len(A), B_arrays=len(B), A_multirun=len(multi),
                        **{f"A_{k}": int((multi.kind == k).sum()) for k in ("1:1", "A split by B", "A merged by B", "complex")},
                        A_mixed_family=int((A.n_fam > 1).sum()), A_with_assembly_gap=int(A.gap_inside.sum()),
                        **{f"breaks_inside_A: {k}": v for k, v in why.items()},
                        median_A_bp=int(A.bp.median()), median_B_bp=int(B.bp.median())))
        if G == 20000:
            for chrom, g in C[C.putative_centromeric_array_L1 == "yes"].groupby("chrom"):
                a_best = g.groupby("array_id").bp.sum().idxmax(); b_best = g.groupby(f"B{G}").bp.sum().idxmax()
                ra, rb = g[g.array_id == a_best], g[g[f"B{G}"] == b_best]
                CEN.append(dict(species=sp, chrom=chrom, A_array=a_best, A_runs=len(ra), A_kb=round(ra.bp.sum() / 1e3), A_families=",".join(sorted(ra.family_v46.unique())),
                                B_array=b_best, B_runs=len(rb), B_kb=round(rb.bp.sum() / 1e3), B_family=rb.family_v46.iloc[0],
                                same_runs=set(ra.run_id) == set(rb.run_id), B_share_of_A=round(len(set(ra.run_id) & set(rb.run_id)) / len(ra), 2)))
S = pd.DataFrame(SUM).fillna(0); S.to_csv(os.path.join(O, "summary.tsv"), sep="\t", index=False)
Z = pd.DataFrame(CEN); Z.to_csv(os.path.join(O, "centromeres.tsv"), sep="\t", index=False)
pd.set_option("display.width", 260); pd.set_option("display.max_columns", 30)
print(S.T.to_string()); print()
print(Z[~Z.same_runs].to_string(index=False)); print("\nmajor centromere identical under A and B (G = 20 kb):", Z.groupby("species").same_runs.sum().to_dict(), "of", Z.groupby("species").size().to_dict())