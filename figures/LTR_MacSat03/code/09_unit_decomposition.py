"""Decompose the Mac centromeric composite unit (2026-10-23; user: each chromosome has its own period, so something in the interval differs).
Unit = interval from one centromeric element core copy to the next on the same chromosome (8-25 kb apart).
Per unit: satellite bp (validated runs, satellite class, S1/S2; MacSat03 bp and monomer count), element bp (hits of the full element consensus excluding
its built-in satellite stretch 2,929-3,388), other bp. Per chromosome: medians and regression of period on parts ('subtraction').
Calibrated nodes: the element consensus cut into 100-bp bins (shared coordinates); for every chromosome, the share of its units containing each bin.
Outputs in unit/: units.tsv, per_chromosome.tsv, node_presence.tsv, summary.txt."""
import pathlib, numpy as np, pandas as pd
W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22"); R = pathlib.Path("F:/data/satellite_FISH/runs_of_record_v4.5.3_2026-10-18/Mac")
(W / "unit").mkdir(exist_ok=True)
cp = pd.read_csv(W / "core" / "Mac.core_copies.tsv", sep="\t"); cp = cp[cp.in_cen_1Mb].sort_values(["chrom", "lo"])
runs = pd.read_csv(R / "runs_of_record.tsv", sep="\t"); runs = runs[(runs.unit_class == "satellite") & runs.tier.isin(["S1", "S2"])]
hits = pd.read_csv(W / "copies" / "Mac.hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
hits = hits[~((hits.q0 >= 2900) & (hits.q1 <= 3420))]                     # drop hits that are only the satellite stretch inside the consensus
hits["lo"] = hits[["s0", "s1"]].min(1); hits["hi"] = hits[["s0", "s1"]].max(1)
NB = 8164 // 100 + 1


def union_bp(iv, a, b):
    iv = sorted((max(x, a), min(y, b)) for x, y in iv if y > a and x < b); t = 0; cs = ce = None
    for x, y in iv:
        if ce is None or x > ce:
            if ce is not None: t += ce - cs
            cs, ce = x, y
        else: ce = max(ce, y)
    return t + (ce - cs if ce is not None else 0)


units = []; nodes = {}
for ch, g in cp.groupby("chrom"):
    lo = g.lo.values
    rc = runs[runs.chrom == ch]; hc = hits[hits.chrom == ch]
    for a, b in zip(lo[:-1], lo[1:]):
        P = b - a
        if not 8000 <= P <= 25000: continue
        r = rc[(rc.end > a) & (rc.start0 < b)]
        sat = union_bp(zip(r.start0, r.end), a, b)
        m3 = union_bp(zip(r[r.family == "MacSat03"].start0, r[r.family == "MacSat03"].end), a, b)
        h = hc[(hc.hi > a) & (hc.lo < b)]
        el = union_bp(zip(h.lo, h.hi), a, b)
        bins = np.zeros(NB, bool)
        for q0, q1 in zip(h.q0, h.q1): bins[q0 // 100: q1 // 100 + 1] = True
        units.append(dict(chrom=ch, start=a, period=P, satellite_bp=sat, MacSat03_bp=m3, MacSat03_monomers=round(m3 / 224, 1), element_bp=el,
                          other_bp=max(0, P - sat - el), element_bins=int(bins.sum())))
        nodes.setdefault(ch, []).append(bins)
U = pd.DataFrame(units); U.to_csv(W / "unit" / "units.tsv", sep="\t", index=False)
per = U.groupby("chrom").agg(units=("period", "size"), period=("period", "median"), satellite_bp=("satellite_bp", "median"),
                             MacSat03_monomers=("MacSat03_monomers", "median"), element_bp=("element_bp", "median"), other_bp=("other_bp", "median"),
                             element_bins=("element_bins", "median")).reset_index()
per = per[per.units >= 5].sort_values("period"); per.to_csv(W / "unit" / "per_chromosome.tsv", sep="\t", index=False)
# subtraction: how much of the between-chromosome period variation each part explains (per-chromosome medians)
var = per.period.var(); expl = {k: np.corrcoef(per.period, per[k])[0, 1] ** 2 if per[k].std() > 0 else 0 for k in ("satellite_bp", "element_bp", "other_bp")}
slope = np.polyfit(per.satellite_bp, per.period, 1)[0] if per.satellite_bp.std() > 0 else np.nan
# calibrated nodes: share of units per chromosome carrying each 100-bp element bin
NP = pd.DataFrame({ch: np.mean(np.vstack(v), 0) for ch, v in nodes.items() if len(v) >= 5}, index=[f"{i * 100}-{i * 100 + 99}" for i in range(NB)]).T
NP.round(2).to_csv(W / "unit" / "node_presence.tsv", sep="\t")
core = [c for c in NP.columns if (NP[c] >= 0.8).all()]
variable = [c for c in NP.columns if (NP[c] >= 0.8).any() and (NP[c] <= 0.2).any()]
txt = [f"units (8-25 kb, centromeric, consecutive): {len(U)} on {U.chrom.nunique()} chromosomes", "",
       "per chromosome (medians):", per.to_string(index=False), "",
       f"share of between-chromosome period variance tracked by each part (r^2 of medians): " + ", ".join(f"{k} {v:.2f}" for k, v in expl.items()),
       f"slope period ~ satellite_bp: {slope:.2f} (1.0 = every extra satellite bp adds one bp of period)", "",
       f"calibrated element nodes (100-bp bins): present in >= 80% of units on every chromosome: {len(core)}; "
       f"chromosome-variable (>= 80% on some, <= 20% on others): {len(variable)}",
       "variable nodes: " + ", ".join(variable[:40])]
(W / "unit" / "summary.txt").write_text("\n".join(txt) + "\n", encoding="utf-8"); print("\n".join(txt))
