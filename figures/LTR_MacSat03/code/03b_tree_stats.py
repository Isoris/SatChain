"""Read the copy tree (2026-10-22): nearest neighbour of every copy (patristic), species mixing, clade structure, within/between-species distances
(also split by centromeric vs arm copies). Writes tree/tree_stats.txt and tree/nearest_neighbours.tsv."""
import pathlib, itertools, collections
import numpy as np
from Bio import Phylo

W = pathlib.Path("/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22/tree")
t = Phylo.read(str(W / "copies.treefile"), "newick"); t.root_at_midpoint()
tips = t.get_terminals(); names = [x.name for x in tips]; n = len(tips)
sp = [nm.split("|")[0] for nm in names]; pos = [nm.split("|")[3] for nm in names]
# patristic distances via depths and LCA (n ~ 400: fine)
D = np.zeros((n, n))
for i, j in itertools.combinations(range(n), 2):
    D[i, j] = D[j, i] = t.distance(tips[i], tips[j])
np.fill_diagonal(D, np.inf)
nn = D.argmin(1)
mix = np.mean([sp[i] != sp[nn[i]] for i in range(n)])
rows = ["copy\tspecies\tposition\tnearest\tnearest_species\tdistance"] + [f"{names[i]}\t{sp[i]}\t{pos[i]}\t{names[nn[i]]}\t{sp[nn[i]]}\t{D[i, nn[i]]:.4f}" for i in range(n)]
(W / "nearest_neighbours.tsv").write_text("\n".join(rows) + "\n")
np.fill_diagonal(D, np.nan)
def med(a, b, pa=None, pb=None):
    I = [i for i in range(n) if sp[i] == a and (pa is None or pos[i] == pa)]; J = [j for j in range(n) if sp[j] == b and (pb is None or pos[j] == pb)]
    v = D[np.ix_(I, J)]; return np.nanmedian(v), len(I), len(J)
# largest single-species clades (support >= 95)
mono = []
for cl in t.get_nonterminals():
    lv = cl.get_terminals(); s = {x.name.split("|")[0] for x in lv}
    if len(s) == 1 and len(lv) >= 10 and (cl.confidence or 0) >= 90: mono.append((len(lv), s.pop(), cl))
mono.sort(key=lambda x: -x[0]); top = []
for k, s_, cl in mono:                                   # keep maximal clades only
    if not any(cl in [c for c in big.find_clades()] for _, _, big in top): top.append((k, s_, cl))
null = []
rng = np.random.default_rng(1)
for _ in range(200):                                     # mixing expected if species labels were random
    p = rng.permutation(sp); null.append(np.mean([p[i] != p[nn[i]] for i in range(n)]))
out = [f"copies in tree: {n} (Gar {sp.count('Gar')}, Mac {sp.count('Mac')})",
       f"copies whose nearest neighbour is from the other species: {mix * 100:.1f}% (random labels: {np.mean(null) * 100:.1f}%, 95% range {np.percentile(null, 2.5) * 100:.1f}-{np.percentile(null, 97.5) * 100:.1f}%)",
       "median patristic distance:",
       "  Gar-Gar %.3f | Mac-Mac %.3f | Gar-Mac %.3f" % (med('Gar', 'Gar')[0], med('Mac', 'Mac')[0], med('Gar', 'Mac')[0]),
       "  centromeric copies only: Gar-Gar %.3f | Mac-Mac %.3f | Gar-Mac %.3f" % (med('Gar', 'Gar', 'cen', 'cen')[0], med('Mac', 'Mac', 'cen', 'cen')[0], med('Gar', 'Mac', 'cen', 'cen')[0]),
       "  Gar arm copies vs Gar centromeric: %.3f" % med('Gar', 'Gar', 'arm', 'cen')[0],
       f"single-species clades with >= 10 copies and SH-aLRT >= 90 (maximal): {len(top)}",
       *[f"  {s_}: {k} copies" for k, s_, _ in top[:10]],
       f"copies inside such clades: Gar {sum(k for k, s_, _ in top if s_ == 'Gar')} of {sp.count('Gar')}, Mac {sum(k for k, s_, _ in top if s_ == 'Mac')} of {sp.count('Mac')}"]
(W / "tree_stats.txt").write_text("\n".join(out) + "\n"); print("\n".join(out))
