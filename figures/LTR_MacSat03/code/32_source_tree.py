"""Step 32a (2026-10-24): source lineage of the Mac centromeric element in the core tree (716 copies, core 5,560-7,400).
Rooting: midpoint and, as a check, on the largest Gar clade. Reports the clade that holds the Mac centromeric copies, the copies nested
inside it that are NOT centromeric, and its sister lineage."""
import collections, pathlib
from Bio import Phylo
W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22")
T = Phylo.read(W / "tree" / "copies.treefile", "newick")
tips = T.get_terminals()
pos = collections.Counter(t.name.split("|")[3] for t in tips); print("position labels:", dict(pos))
def info(t):
    sp, lg, p, where, sat = t.name.split("|"); return sp, lg, where, sat
macc = [t for t in tips if info(t)[0] == "Mac" and info(t)[2] != "arm"]
print("Mac centromeric tips:", len(macc))
def report(T, label):
    print(f"\n=== rooting: {label}")
    macc = [t for t in T.get_terminals() if info(t)[0] == "Mac" and info(t)[2] != "arm"]; M = set(macc)
    mrca = T.common_ancestor(macc); inside = mrca.get_terminals()
    c = collections.Counter((info(t)[0], "cen" if info(t)[2] != "arm" else "arm") for t in inside)
    print("smallest clade with all Mac centromeric tips:", len(inside), "tips", dict(c))
    # largest clade that is >= 95% Mac-centromeric
    best = None
    for cl in T.find_clades():
        ts = cl.get_terminals(); n = sum(1 for t in ts if t in M)
        if len(ts) >= 10 and n / len(ts) >= 0.95 and (best is None or n > best[1]): best = (cl, n, len(ts))
    cl, n, N = best; print(f"largest >=95% Mac-centromeric clade: {n} of {N} tips, support {cl.confidence}")
    path = T.get_path(cl); 
    for depth, anc in enumerate(reversed(path[:-1][-4:]), 1):
        ts = anc.get_terminals(); C = set(cl.get_terminals()); sis = [t for t in ts if t not in C]
        cc = collections.Counter((info(t)[0], "cen" if info(t)[2] != "arm" else "arm", info(t)[1]) for t in sis)
        print(f"  ancestor -{depth} (support {anc.confidence}): +{len(sis)} tips:", cc.most_common(8))
    return cl
def midpoint(T):
    T2 = Phylo.read(W / "tree" / "copies.treefile", "newick"); T2.root_at_midpoint(); return T2
cl = report(midpoint(T), "midpoint")
T3 = Phylo.read(W / "tree" / "copies.treefile", "newick")
gar = [t for t in T3.get_terminals() if info(t)[0] == "Gar"]
# outgroup = the Gar 266-copy clade: largest clade that is 100% Gar
bestg = max((c for c in T3.find_clades() if all(info(t)[0] == "Gar" for t in c.get_terminals()) and c != T3.root), key=lambda c: c.count_terminals())
T3.root_with_outgroup(bestg); report(T3, f"Gar clade ({bestg.count_terminals()} tips)")
# arm copies nearest to the centromeric clade MRCA (patristic), both species
Tm = midpoint(T); clm = Tm.common_ancestor([t for t in Tm.get_terminals() if t.name in {x.name for x in cl.get_terminals()}])
inner = set(t.name for t in clm.get_terminals())
d = []
for t in Tm.get_terminals():
    if t.name in inner and info(t)[2] != "arm": continue
    if info(t)[2] == "arm" or info(t)[0] == "Gar": d.append((Tm.distance(clm, t), t.name, t.name in inner))
d.sort()
print("\nnon-centromeric / Gar copies closest to the Mac centromeric clade ancestor (patristic, inside-clade flag):")
for x in d[:15]: print(f"  {x[0]:.4f}  {x[1]}  {'INSIDE' if x[2] else ''}")
open(W / "tree" / "source_candidates_tree.tsv", "w").write("distance_to_cen_mrca\tcopy\tinside_cen_clade\n" + "".join(f"{a:.5f}\t{b}\t{c}\n" for a, b, c in d))
