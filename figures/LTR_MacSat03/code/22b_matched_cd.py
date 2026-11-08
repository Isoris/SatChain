"""Step 22b (2026-10-24): fair L-c vs L-d. For unit i, another-chromosome unit j and a Gar copy g, divergence i-j and i-g on the SAME
element windows (present in all three). Also the nearest Gar copy (min over 20) vs nearest other-chromosome unit (min over 20).
Output: homologous/matched_cd.tsv + summary per window set."""
import pathlib
src = open(pathlib.Path(__file__).parent / "22_per_unit_levels.py", encoding="utf-8").read()
exec(src[:src.index("rows = []")])
def div(fa_, fb_, ws):
    P = Q = n = 0
    for w in ws:
        p, q, m = pq(fa_["win"][w]["seq"], fb_["win"][w]["seq"]); P += p; Q += q; n += m
    return k2p(P, Q, n), n
rows = []
for i, r in U.iterrows():
    others = random.sample(list(U.index[U.chrom != r["chrom"]]), 20); gs = random.sample(range(len(G)), 20)
    for j, g in zip(others, gs):
        ws = sorted(set(F[i]["win"]) & set(F[j]["win"]) & set(G[g]["win"]))
        if not ws: continue
        dij, nij = div(F[i], F[j], ws); dig, nig = div(F[i], G[g], ws)
        rows.append(dict(unit=i, chrom=lab(r["chrom"]), type=r["type"], other_type=U.loc[j, "type"], gar_copy=g, windows=len(ws),
                         first_window=ws[0], last_window=ws[-1], mac_other_chrom=dij, gar=dig, n_cols=min(nij, nig)))
    if i % 25 == 0: print("unit", i, flush=True)
M = pd.DataFrame(rows); M.to_csv(OUT / "matched_cd.tsv", sep="\t", index=False)
M["gar_minus_mac"] = M.gar - M.mac_other_chrom
# per unit medians first (units are the measurement), then across units; also per centromere (n = 12)
pu = M.groupby(["chrom", "unit"])[["mac_other_chrom", "gar", "gar_minus_mac"]].median()
pc = pu.groupby("chrom").median()
print("per unit median, then median over units:\n", pu.median().round(4).to_string())
print("per centromere:\n", pc.round(4).to_string())
print("centromeres where Gar > other Mac chromosome:", int((pc.gar_minus_mac > 0).sum()), "of", len(pc))
print("windows used (median):", M.windows.median(), "| window start range:", M.first_window.min(), "-", M.last_window.max())
nr = M.groupby("unit").apply(lambda d: pd.Series(dict(nearest_mac=d.mac_other_chrom.min(), nearest_gar=d.gar.min())))
print("nearest (min over 20) per unit, median:\n", nr.median().round(4).to_string(), "\nunits whose nearest Gar copy beats nearest other-chromosome unit:", int((nr.nearest_gar < nr.nearest_mac).sum()), "of", len(nr))
