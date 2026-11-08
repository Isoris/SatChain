"""Step 36a (2026-10-24): units defined from the satellite side. Boundary = start of a v4.6 MacSat03 run (copy-exact, HMM-trimmed);
unit = run_i + the following interruption, up to the start of run_{i+1} (only when that interruption is the v4.6 interruption between them).
Compared with the element-anchored units of step 20 (anchor = EDTA consensus position 6,500)."""
import pathlib, pandas as pd, numpy as np
D = pathlib.Path("F:/data/fClaHyb_func_annotation_v1.0/fClaHyb_Mac/fClaHyb_Mac_satelliteDNA"); W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22")
C = pd.read_csv(D / "fClaHyb_Mac.satDNA_v46.catalogue.tsv", sep="\t"); I = pd.read_csv(D / "fClaHyb_Mac.satDNA_v46.interruptions_between_runs.tsv", sep="\t")
C = C.set_index("run_id"); I = I[I.run_left.map(C.family).eq("MacSat03") & I.run_right.map(C.family).eq("MacSat03")]
F = pd.read_csv(W / "split" / "v46_fragments.tsv", sep="\t")
rows = []
for r in I.itertuples():
    L, R = C.loc[r.run_left], C.loc[r.run_right]
    rows.append(dict(chrom=r.chrom.replace("C_mac_", ""), unit_lo=L.start0, unit_hi=R.start0, sat_bp=L.end - L.start0, sat_copies=L.n_conf, spacer_bp=r.bp))
U = pd.DataFrame(rows); U["period"] = U.unit_hi - U.unit_lo
A = pd.read_csv(W / "unit_anchored" / "units_A6500.tsv", sep="\t"); A["chrom"] = A.chrom.str.replace("C_mac_", ""); A["period"] = A.hi - A.lo
def stats(x): return f"n={len(x)}, median {x.median():.0f}, IQR {x.quantile(.25):.0f}-{x.quantile(.75):.0f}, CV {x.std() / x.mean():.2f}"
for ch in ("LG04", "LG01", "LG07", "LG19", "LG08", "LG25"):
    print(f"{ch}  satellite-side units: {stats(U[U.chrom == ch].period)} | spacer {U[U.chrom == ch].spacer_bp.median():.0f} bp, satellite {U[U.chrom == ch].sat_bp.median():.0f} bp ({U[U.chrom == ch].sat_copies.median():.0f} copies)")
    print(f"      element-anchored units: {stats(A[A.chrom == ch].period)}")
print("\nall MacSat03 centromeres, satellite-side unit period CV:", U.groupby("chrom").period.apply(lambda x: round(x.std() / x.mean(), 2)).to_dict())
U.to_csv(W / "split" / "units_satellite_side.tsv", sep="\t", index=False)
