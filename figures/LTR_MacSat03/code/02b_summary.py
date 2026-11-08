"""Chain blastn HSPs into element copies (same chromosome and strand, consecutive on both query and subject within 2 kb), then summarise per haplotype:
copies by consensus coverage, identity, chromosome, and overlap with satellite runs of record v4.5.3 (lineage A families: MacSat03/51/162, GarSat03/09/10).
Writes copies/<hap>.copies.tsv and copies/summary.txt (2026-10-22)."""
import pandas as pd, numpy as np, pathlib
W = pathlib.Path("/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22"); R = pathlib.Path("/mnt/f/data/satellite_FISH/runs_of_record_v4.5.3_2026-10-18")
QL = 8164; LINA = {"Mac": ["MacSat03", "MacSat51", "MacSat162"], "Gar": ["GarSat03", "GarSat09", "GarSat10"]}
out = []
for h in ("Mac", "Gar"):
    t = pd.read_csv(W / "copies" / f"{h}.hits.tsv", sep="\t", header=None, names=["chrom", "s0", "s1", "q0", "q1", "pid", "len", "e", "bits"])
    t["strand"] = np.where(t.s1 >= t.s0, "+", "-"); t["lo"] = t[["s0", "s1"]].min(1); t["hi"] = t[["s0", "s1"]].max(1)
    copies = []
    for (c, st), g in t.sort_values("lo").groupby(["chrom", "strand"]):
        cur = None
        for r in g.itertuples():
            # collinear chaining (fixed 2026-10-23): along the subject, query coordinates must advance (+ strand) or recede (- strand) by <= 2 kb steps
            colin = cur and ((st == "+" and -200 <= r.q0 - cur["qhi"] <= 2000) or (st == "-" and -200 <= cur["qlo"] - r.q1 <= 2000))
            if cur and r.lo - cur["hi"] <= 2000 and colin:
                cur["hi"] = max(cur["hi"], r.hi); cur["qlo"] = min(cur["qlo"], r.q0); cur["qhi"] = max(cur["qhi"], r.q1); cur["aln"] += r.len; cur["idw"] += r.pid * r.len
            else:
                if cur: copies.append(cur)
                cur = dict(chrom=c, strand=st, lo=r.lo, hi=r.hi, qlo=r.q0, qhi=r.q1, aln=r.len, idw=r.pid * r.len)
        if cur: copies.append(cur)
    cp = pd.DataFrame(copies); cp["cov"] = (cp.qhi - cp.qlo + 1) / QL; cp["pid"] = cp.idw / cp.aln; cp["bp"] = cp.hi - cp.lo + 1
    rr = pd.read_csv(R / h / "runs_of_record.tsv", sep="\t"); rr = rr[(rr.unit_class == "satellite") & rr.tier.isin(["S1", "S2"])]
    la = rr[rr.family.isin(LINA[h])]
    def near(row, df, pad=5000):
        d = df[(df.chrom == row.chrom) & (df.end > row.lo - pad) & (df.start0 < row.hi + pad)]
        return ",".join(sorted(set(d.family)))[:60]
    cp["satellite_within_5kb"] = [near(r, rr) for r in cp.itertuples()]
    cp["lineageA_within_5kb"] = [bool(near(r, la)) for r in cp.itertuples()]
    cp.drop(columns=["idw"]).to_csv(W / "copies" / f"{h}.copies.tsv", sep="\t", index=False)
    cp["near_full"] = (cp["cov"] >= 0.7) & (cp.bp >= 0.6 * QL) & (cp.bp <= 1.4 * QL)     # length bound added 2026-10-23
    full = cp[cp.near_full]; part = cp[(cp["cov"] >= 0.1) & ~cp.near_full]
    out += [f"== {h} (fClaHyb_{h}_LG)",
            f"copies (chained HSPs): {len(cp)}; total bp {cp.bp.sum() / 1e6:.2f} Mb",
            f"near-full (>= 70% of consensus): {len(full)}; median identity {full.pid.median():.1f}%; partial (10-70%): {len(part)}",
            f"near-full copies next to lineage A satellite (<= 5 kb): {int(full.lineageA_within_5kb.sum())} of {len(full)}",
            "near-full copies per chromosome (top 8): " + ", ".join(f"{c.split('_')[-1]}:{n}" for c, n in full.chrom.value_counts().head(8).items()),
            "identity of near-full copies, quartiles: " + ", ".join(f"{q:.1f}" for q in full.pid.quantile([0.25, 0.5, 0.75])) if len(full) else "no near-full copies", ""]
(W / "copies" / "summary.txt").write_text("\n".join(out) + "\n"); print("\n".join(out))
