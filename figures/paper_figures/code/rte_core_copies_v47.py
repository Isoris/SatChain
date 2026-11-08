"""RTE core copies inside / outside the centromere, by the family group of the centromere (data behind Figure 2c), on the SatChain v4.7 runs of record (2026-11-01).
Same copy chaining and the same sweep as panels_rte/panel_c.py and thesis_ready/code/rte_inside_outside_v46.py (code of 06b_core_copies.py), but:
  centromere intervals and dominant family  = satchain_v4.7_2026-11-01/runs_of_record/{Mac,Gar}/centromere_check.tsv (cen_start0, cen_end, dominant_family; v4.7 = same intervals as v4.6);
  group of a centromere                     = v4.7 seq_cluster of its dominant family: C2.1 -> "cen" (the GarSat03 / MacSat03 group, the RTE-carrying centromeric families),
                                              C2.2 -> "5S" (GarSat04, the 5S-derived satellite), anything else (gene-array-dominated, 'none', unsupported) -> "other";
  runs near a copy                          = v4.7 runs of record, unit_class satellite, tier S1 (v4.7 has no S2 tier); "near cen" = a C2.1 run within 5 kb.
The chromosomes without v4.7 rows (Mac LG14 / 17 / 18, Gar LG23 / 27 / 28) have dominant_family 'none' and fall in "other".
Reads: LTR_MacSat03_2026-10-22/core/{Mac,Gar}.core_hits.tsv (BLAST hits of the RT-bearing core 5,560-7,400, 1,841 bp).  No work at import."""
import sys, pathlib, numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); import satellite_style as ST

W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22")
CL = 1841
COVS = (0.5, 0.7, 0.9)          # copy definition: share of the 1,841-bp core covered by the chained hits
PADS = (0.5e6, 1e6, 2e6)        # window around the centromere interval
REF_COV, REF_PAD = 0.7, 1e6
GROUPS = ("cen", "5S", "other")


def chain(h):
    t = pd.read_csv(W / "core" / f"{h}.core_hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
    t["strand"] = np.where(t.s1 >= t.s0, "+", "-"); t["lo"] = t[["s0", "s1"]].min(axis=1); t["hi"] = t[["s0", "s1"]].max(axis=1)
    cps = []
    for (c, st), g in t.sort_values("lo").groupby(["chrom", "strand"]):
        cur = None
        for r in g.itertuples():
            colin = cur is not None and ((st == "+" and -100 <= r.q0 - cur["qhi"] <= 500) or (st == "-" and -100 <= cur["qlo"] - r.q1 <= 500))
            if cur is not None and r.lo - cur["hi"] <= 500 and colin:
                cur.update(hi=max(cur["hi"], r.hi), qlo=min(cur["qlo"], r.q0), qhi=max(cur["qhi"], r.q1))
            else:
                if cur: cps.append(cur)
                cur = dict(chrom=c, strand=st, lo=r.lo, hi=r.hi, qlo=r.q0, qhi=r.q1)
        if cur: cps.append(cur)
    cp = pd.DataFrame(cps); cp["cov"] = (cp.qhi - cp.qlo + 1) / CL; cp["bp"] = cp.hi - cp.lo + 1
    return cp


def group_of(fam):
    c = ST.family_cluster(str(fam))
    return "cen" if c == "C2.1" else "5S" if c == "C2.2" else "other"


def load(h):
    cp = chain(h)
    rr = ST.runs(h); rr = rr[(rr.unit_class == "satellite") & rr.tier.isin(["S1", "S2"])].copy(); rr["grp"] = rr.seq_cluster.map({"C2.1": "cen", "C2.2": "5S"}).fillna("other")
    cen = pd.read_csv(ST.RUNS / h / "centromere_check.tsv", sep="\t").set_index("chrom")
    cen["group"] = [group_of(f) for f in cen.dominant_family]
    f = cp[(cp["cov"] >= min(COVS)) & (cp.bp <= 1.5 * CL)].copy()          # superset of every definition
    near = []
    for r in f.itertuples():
        d = rr[(rr.chrom == r.chrom) & (rr.end > r.lo - 5000) & (rr.start0 < r.hi + 5000)]; near.append((d.grp == "cen").any())
    f["cen_5kb"] = near
    f["cen_group"] = [cen.loc[c, "group"] if c in cen.index else "NA" for c in f.chrom]
    f["cs"] = [cen.loc[c, "cen_start0"] if c in cen.index else np.nan for c in f.chrom]; f["ce"] = [cen.loc[c, "cen_end"] if c in cen.index else np.nan for c in f.chrom]
    f["hap"] = h
    return f, cen


def rows(f, cen, h, cov, pad):
    out = []
    g = f[f["cov"] >= cov]
    for grp in GROUPS:
        chs = cen.index[cen.group == grp]; gg = g[g.chrom.isin(chs)]
        inside = (gg.lo < gg.ce + pad) & (gg.hi > gg.cs - pad)
        out.append(dict(haplotype=h, centromere_group=grp, cov_min=cov, window_Mb=pad / 1e6, n_chromosomes=len(chs), chromosomes_with_copies=gg.chrom.nunique(),
                        chromosomes_with_copy_inside=gg[inside].chrom.nunique(), core_copies=len(gg), inside=int(inside.sum()), outside=int((~inside).sum()),
                        pct_inside=100 * inside.sum() / len(gg) if len(gg) else np.nan, next_to_cen=int(gg.cen_5kb.sum()),
                        pct_next_to_cen=100 * gg.cen_5kb.sum() / len(gg) if len(gg) else np.nan, chromosomes=";".join(c.split("_")[-1] for c in chs)))
    return out


def table():
    allrows, data = [], {}
    for h in ("Mac", "Gar"):
        f, cen = load(h); data[h] = (f, cen)
        for cov in COVS:
            for pad in PADS: allrows += rows(f, cen, h, cov, pad)
    return pd.DataFrame(allrows), data
