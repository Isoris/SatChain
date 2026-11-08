import collections, pathlib, pandas as pd
from Bio import Phylo
W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22")
T = Phylo.read(W / "tree" / "copies.treefile", "newick")
bestg = max((c for c in T.find_clades() if all(t.name.startswith("Gar") for t in c.get_terminals()) and c != T.root), key=lambda c: c.count_terminals())
T.root_with_outgroup(bestg)
macc = [t for t in T.get_terminals() if t.name.startswith("Mac") and "|cen|" in t.name]
T2 = [t.name for t in macc if t.name.split("|")[1] in ("LG08", "LG25", "LG16")]
def comp(cl):
    c = collections.Counter()
    for t in cl.get_terminals():
        sp, lg, p, w, s = t.name.split("|")
        c["Gar" if sp == "Gar" else ("Mac-cen-II" if t.name in T2 else "Mac-cen-I" if w == "cen" else "Mac-arm")] += 1
    return dict(c)
# spine: from the smallest clade holding all Mac centromeric copies, walk into the child that holds most of them; print the other child
node = T.common_ancestor(macc); lvl = 0
arm_rows = []
while not node.is_terminal() and lvl < 40:
    kids = sorted(node.clades, key=lambda k: sum(1 for t in k.get_terminals() if t in set(macc)), reverse=True)
    main, others = kids[0], kids[1:]
    for o in others:
        print(f"lvl {lvl:2d} support {node.confidence}: off-spine clade {comp(o)}  | main continues with {comp(main)}")
        for t in o.get_terminals():
            if "|arm|" in t.name and t.name.startswith("Mac"): arm_rows.append((lvl, t.name, T.distance(node, t)))
    if sum(1 for t in main.get_terminals() if t in set(macc)) < 0.5 * len(macc): break
    node = main; lvl += 1
pd.DataFrame(arm_rows, columns=["spine_level", "copy", "dist_from_spine_node"]).to_csv(W / "tree" / "spine_arm_copies.tsv", sep="\t", index=False)
