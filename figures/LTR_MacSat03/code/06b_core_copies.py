"""Core copies (2026-10-23): collinear chains of core hits (core = consensus 5,560-7,400, 1,841 bp) covering >= 70% of the core and spanning
<= 1.5x core. Per haplotype: core/<hap>.core_copies.tsv, core/<hap>.core_regions.tsv (for extraction), core/summary.txt with counts, identity,
lineage A satellite neighbourhood and the centromere-lineage table (copies inside centromere +-1 Mb by satellite lineage of that centromere)."""
import pandas as pd, numpy as np, pathlib
W = pathlib.Path("/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22"); R = pathlib.Path("/mnt/f/data/satellite_FISH/runs_of_record_v4.5.3_2026-10-18")
CL = 1841
LIN = {"GarSat03": "A", "GarSat09": "A", "GarSat10": "A", "MacSat03": "A", "MacSat51": "A", "MacSat162": "A", "GarSat04": "B", "GarSat49": "B", "MacSat06": "B"}
out = []; tab = []
for h in ("Mac", "Gar"):
    t = pd.read_csv(W / "core" / f"{h}.core_hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
    t["strand"] = np.where(t.s1 >= t.s0, "+", "-"); t["lo"] = t[["s0", "s1"]].min(1); t["hi"] = t[["s0", "s1"]].max(1)
    cps = []
    for (c, st), g in t.sort_values("lo").groupby(["chrom", "strand"]):
        cur = None
        for r in g.itertuples():
            colin = cur is not None and ((st == "+" and -100 <= r.q0 - cur["qhi"] <= 500) or (st == "-" and -100 <= cur["qlo"] - r.q1 <= 500))
            if cur is not None and r.lo - cur["hi"] <= 500 and colin:
                cur.update(hi=max(cur["hi"], r.hi), qlo=min(cur["qlo"], r.q0), qhi=max(cur["qhi"], r.q1), aln=cur["aln"] + r.len, idw=cur["idw"] + r.pid * r.len)
            else:
                if cur: cps.append(cur)
                cur = dict(chrom=c, strand=st, lo=r.lo, hi=r.hi, qlo=r.q0, qhi=r.q1, aln=r.len, idw=r.pid * r.len)
        if cur: cps.append(cur)
    cp = pd.DataFrame(cps); cp["cov"] = (cp.qhi - cp.qlo + 1) / CL; cp["pid"] = cp.idw / cp.aln; cp["bp"] = cp.hi - cp.lo + 1
    cp["core_copy"] = (cp["cov"] >= 0.7) & (cp.bp <= 1.5 * CL)
    rr = pd.read_csv(R / h / "runs_of_record.tsv", sep="\t"); rr = rr[(rr.unit_class == "satellite") & rr.tier.isin(["S1", "S2"])]
    cen = pd.read_csv(R / h / "centromere_check.tsv", sep="\t").set_index("chrom")
    f = cp[cp.core_copy].copy()
    def near(r, pad=5000):
        d = rr[(rr.chrom == r.chrom) & (rr.end > r.lo - pad) & (rr.start0 < r.hi + pad)]; return ",".join(sorted(set(d.family)))
    f["satellite_within_5kb"] = [near(r) for r in f.itertuples()]
    f["lineageA_within_5kb"] = f.satellite_within_5kb.apply(lambda s: any(LIN.get(x) == "A" for x in s.split(",") if x))
    f["in_cen_1Mb"] = [r.chrom in cen.index and r.lo < cen.loc[r.chrom, "cen_end"] + 1e6 and r.hi > cen.loc[r.chrom, "cen_start"] - 1e6 for r in f.itertuples()]
    f["cen_lineage"] = [LIN.get(cen.loc[r.chrom, "top_family_in_cen"], "other") if r.chrom in cen.index else "NA" for r in f.itertuples()]
    f.drop(columns=["idw"]).to_csv(W / "core" / f"{h}.core_copies.tsv", sep="\t", index=False)
    with open(W / "core" / f"{h}.core_regions.tsv", "w") as o:
        for r in f.itertuples():
            fam = cen.loc[r.chrom, "top_family_in_cen"] if r.chrom in cen.index else "NA"
            o.write(f"{r.chrom}:{r.lo}-{r.hi}\t{r.strand}\t{h}|{r.chrom.split('_')[-1]}|{r.lo}|{'cen' if r.in_cen_1Mb else 'arm'}|{fam}\n")
    out += [f"== {h}: core hits chained into {len(cp)} segments; core copies (>= 70% of core, span <= 1.5x): {len(f)}",
            f"identity to the (Mac-built) core: median {f.pid.median():.1f}% (quartiles {f.pid.quantile(.25):.1f}-{f.pid.quantile(.75):.1f})",
            f"next to lineage A satellite (<= 5 kb): {int(f.lineageA_within_5kb.sum())} of {len(f)}; inside a centromere (+-1 Mb): {int(f.in_cen_1Mb.sum())}",
            "per chromosome (top 8): " + ", ".join(f"{c.split('_')[-1]}:{n}" for c, n in f.chrom.value_counts().head(8).items()), ""]
    for lin in ("A", "B", "other"):
        ch = [c for c in cen.index if LIN.get(cen.loc[c, "top_family_in_cen"], "other") == lin]; ff = f[f.chrom.isin(ch)]
        tab.append((h, lin, len(ch), len(ff), int(ff.in_cen_1Mb.sum()), int(ff.lineageA_within_5kb.sum())))
tb = pd.DataFrame(tab, columns=["haplotype", "centromere_lineage", "n_chromosomes", "core_copies", "in_centromere_pm1Mb", "next_to_lineageA"])
tb.to_csv(W / "core" / "by_centromere_lineage.tsv", sep="\t", index=False)
txt = "\n".join(out) + "\n" + tb.to_string(index=False) + "\n"; (W / "core" / "summary.txt").write_text(txt); print(txt)
