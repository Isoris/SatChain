#!/usr/bin/env python3
"""Figure 1 panel a (92 x 56 mm): monomer trees of the two satellite groups, drawn from the v4.5 IQ-TREE treefiles in the house style, with FOUR resolution levels.
5S-related tree: GarSat04 (5S-derived satellite), the 5S rDNA gene arrays GarSat49 and MacSat06, + 5S rRNA gene references.  Centromeric tree (the GarSat03 / MacSat03 group,
the RTE-carrying centromeric families; SatChain v4.7 sequence cluster C2.1): GarSat03, GarSat09, GarSat10, MacSat03, MacSat162.  The word "lineage" is dropped (SatChain v4.7, 2026-11-01);
the treefile stems (lineageA, lineageB_5S_derived) are the fixed v4.5 input names and are not renamed.
Layout = equal-angle + equal-daylight (copy of thesis_defense_slides_tpl/code_R34/daylight_tree.py, algorithm unchanged; layout cached in panels/_cache).
Level 0  tips: fill = family (palette_v46.family, own colour per family), marker = species (circle Gar, square Mac, diamond 5S gene).
Level 1  sub-clades R1.. (centromeric tree) / S1.. (5S-related tree): maximal subtrees whose patristic diameter <= CUT (substitutions/site) AND whose UFBoot >= 95 AND >= 5 copies.  CUT is swept
         (23 values, 0.10-1.20) and taken as the middle of the longest plateau of identical clade sets (rule in build_model); label = clade name + copies (multi-chromosome clades) or + chromosome.
Level 2  chromosome origin: one ring segment per chromosome array inside each clade (width = copies); hue = species (purple Gar, green Mac); dark = array flagged centromeric in the v4.7
         runs of record (centromeric_L1), pale = other, open = 5S rRNA genes (no array).  Level 3: UFBoot dots on the backbone nodes (nodes that join clades), the node where the 5S gene clade
         attaches to the satellite clades is ringed in gold with its support.  Species mixing = the segment hues of each ring (both hues = intermixed clade).  Evidence level L1.
usage: draw(fig, left_mm, top_mm, w_mm, h_mm, letter) -> dict ; python fig1_panel_a.py renders it alone (panels/_standalone) and writes numbers / sweep / caption."""
import sys, pickle, pathlib, csv, collections
import numpy as np
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "code"))
import satellite_style as ST
from matplotlib.collections import LineCollection
from matplotlib.patches import Wedge

PHY = pathlib.Path(r"F:\data\satellite_FISH\share_LITE_v4.5_2026-10-15\phylogeny")
RUNS = ST.RUNS                                    # SatChain v4.7 runs of record (centromeric_L1 flag of the array covering each copy)
CACHE = HERE / "_cache"
# tree keys: "B" = the 5S-related tree (treefile lineageB_5S_derived), "A" = the centromeric tree (treefile lineageA); labels and clade prefixes carry no lineage wording
LIN = {"B": ("lineageB_5S_derived", "5S-related families"), "A": ("lineageA", "Centromeric families")}
PREFIX = {"B": "S", "A": "R"}                     # clade names: S1.. (5S-related tree), R1.. (centromeric, RTE-carrying tree)
GROUPCOL = {"B": ST.FIVE_S, "A": ST.CEN}
TREE_NAME = {"B": "5S-related", "A": "centromeric"}
GOLD5S = "#e6a817"            # palette_v46.RDNA['5S rDNA clusters']
GOLD_DARK = "#9a6b00"
SPECIES_MARK = {"Gar": "o", "Mac": "s", "5S": "D"}
MODEL_VER = 11                # 11 = v4.7 runs of record for the centromeric flag ('na' = unknown) + R/S clade prefixes (2026-11-01)
CUT_GRID = [round(x, 2) for x in np.arange(0.10, 1.2001, 0.05)]      # patristic-diameter cut sweep (substitutions / site)
SUP_MIN, MIN_SIZE, JACC = 95, 5, 0.8
SUP_GRID = (70, 80, 90, 95, 99)
BURIED = 3.5                  # mm: a clade whose ring would float further than this beyond its own copies is ringed locally


# ---------------------------------------------------------------- layout (copy of daylight_tree.layout, algorithm unchanged)
def _layout(tree, iters=30, min_bl=0.0005):
    clades = list(tree.find_clades(order="preorder")); idx = {id(c): i for i, c in enumerate(clades)}; n = len(clades)
    parent = [-1] * n; kids = [[] for _ in range(n)]
    for c in clades:
        for ch in c.clades: parent[idx[id(ch)]] = idx[id(c)]; kids[idx[id(c)]].append(idx[id(ch)])
    bl = np.array([max(c.branch_length or 0.0, min_bl) if parent[i] >= 0 else 0.0 for i, c in enumerate(clades)])
    desc = [None] * n; nleaf = np.zeros(n)
    for i in range(n - 1, -1, -1):
        d = [i]
        for k in kids[i]: d.extend(desc[k])
        desc[i] = np.array(d); nleaf[i] = 1 if not kids[i] else sum(nleaf[k] for k in kids[i])
    pos = np.zeros((n, 2)); wedge = {0: (0.0, 2 * np.pi)}
    for i in range(n):
        a0, a1 = wedge[i]; cur = a0
        for k in kids[i]:
            span = (a1 - a0) * nleaf[k] / nleaf[i]; mid = cur + span / 2
            pos[k] = pos[i] + bl[k] * np.array([np.cos(mid), np.sin(mid)]); wedge[k] = (cur, cur + span); cur += span
    inner = [i for i in range(n) if len(kids[i]) >= 2]; allnodes = np.arange(n)
    for _ in range(iters):
        for v in inner:
            groups = []; outside = np.setdiff1d(allnodes, desc[v])
            if len(outside): groups.append((outside, True))
            for k in kids[v]: groups.append((desc[k], False))
            if len(groups) < 3: continue
            info = []
            for mem, isp in groups:
                d = pos[mem] - pos[v]; ang = np.arctan2(d[:, 1], d[:, 0]); ref = np.arctan2(*(d.mean(0)[::-1]))
                rel = (ang - ref + np.pi) % (2 * np.pi) - np.pi; info.append([ref + rel.min(), ref + rel.max(), mem, isp])
            info = [info[j] for j in sorted(range(len(info)), key=lambda j: info[j][0] % (2 * np.pi))]
            tot = 2 * np.pi - sum((g[1] - g[0]) for g in info)
            if tot <= 0: continue
            gap = tot / len(info); pi = next((j for j, g in enumerate(info) if g[3]), 0); cur = info[pi][1]
            for off in range(1, len(info)):
                g = info[(pi + off) % len(info)]; new0 = cur + gap; delta = (new0 - g[0] + np.pi) % (2 * np.pi) - np.pi
                c, s = np.cos(delta), np.sin(delta); R = np.array([[c, -s], [s, c]]); mem = g[2]
                pos[mem] = pos[v] + (pos[mem] - pos[v]) @ R.T; cur = new0 + (g[1] - g[0])
    return clades, pos, parent


def load_tree(key):
    """-> dict(names, xy (tip coordinates), segs (edge segments)); layout cached because the daylight passes take ~1 min"""
    stem = LIN[key][0]; cf = CACHE / f"layout_{stem}.pkl"
    if cf.exists(): return pickle.load(open(cf, "rb"))
    from Bio import Phylo
    tree = Phylo.read(str(PHY / stem / f"{stem}.treefile"), "newick")
    clades, pos, parent = _layout(tree, 30)
    tip = [i for i, c in enumerate(clades) if c.is_terminal()]
    d = dict(stem=stem, names=[clades[i].name for i in tip], xy=pos[tip],
             segs=np.array([[pos[parent[i]], pos[i]] for i in range(len(clades)) if parent[i] >= 0]), n_nodes=len(clades),
             tree_len=float(sum((c.branch_length or 0) for c in clades)))
    CACHE.mkdir(exist_ok=True); pickle.dump(d, open(cf, "wb")); return d


def parse(name):
    sp, fam, chrom, pos = name.split("|"); return sp, fam, chrom


def families(d):
    info = [parse(n) for n in d["names"]]
    fams = sorted({i[1] for i in info}, key=lambda f: (f == "5S_rRNA", f[:3], -sum(i[1] == f for i in info)))
    return info, fams


def fcol(f): return GOLD5S if f == "5S_rRNA" else ST.family_colour(f)



# ---------------------------------------------------------------- tree data (layout cache)
def load_tree(key):
    """-> dict(names, xy (tip coordinates), segs (edge segments, row j = edge into preorder node j+1: [parent xy, child xy])); layout cached because the daylight passes take ~1 min"""
    stem = LIN[key][0]; cf = CACHE / f"layout_{stem}.pkl"
    if cf.exists(): return pickle.load(open(cf, "rb"))
    from Bio import Phylo
    tree = Phylo.read(str(PHY / stem / f"{stem}.treefile"), "newick")
    clades, pos, parent = _layout(tree, 30)
    tip = [i for i, c in enumerate(clades) if c.is_terminal()]
    d = dict(stem=stem, names=[clades[i].name for i in tip], xy=pos[tip],
             segs=np.array([[pos[parent[i]], pos[i]] for i in range(len(clades)) if parent[i] >= 0]), n_nodes=len(clades),
             tree_len=float(sum((c.branch_length or 0) for c in clades)))
    CACHE.mkdir(exist_ok=True); pickle.dump(d, open(cf, "wb")); return d


def node_xy(d):
    """layout coordinates of ALL nodes in preorder (node 0 = root), rebuilt from the cached edge list (row j = edge into node j+1)"""
    segs = d["segs"]; P = np.zeros((len(segs) + 1, 2)); P[0] = segs[0][0]; P[1:] = segs[:, 1]; return P


def parse(name):
    sp, fam, chrom, pos = name.split("|"); return sp, fam, chrom


def families(d):
    info = [parse(n) for n in d["names"]]
    fams = sorted({i[1] for i in info}, key=lambda f: (f == "5S_rRNA", f[:3], -sum(i[1] == f for i in info)))
    return info, fams


def fcol(f): return GOLD5S if f == "5S_rRNA" else ST.family_colour(f)


# ---------------------------------------------------------------- sub-clades (level 1), chromosome origin (level 2), backbone support (level 3)
def _tip_centromeric(names):
    """True / False / None per tip: centromeric_L1 of the v4.7 run of record covering the copy (same species and chromosome; the run of the tip's family preferred, any run otherwise because
    v4.6 / v4.7 re-assigned the family of a few arrays); None for 5S rRNA gene references and copies outside every run (incl. the chromosomes absent from the v4.7 tables)"""
    import pandas as pd
    R = {sp: pd.read_csv(RUNS / sp / "runs_of_record.tsv", sep="\t", usecols=["chrom", "start0", "end", "family", "centromeric_L1"]) for sp in ("Gar", "Mac")}
    out = []
    for n in names:
        sp, f, c, p = n.split("|"); p = int(p)
        if f == "5S_rRNA": out.append(None); continue
        r = R[sp]; x = r[(r.chrom == f"C_{sp.lower()}_{c}") & (r.start0 - 1 <= p) & (r.end + 1 >= p)]
        if not len(x): out.append(None); continue
        y = x[x.family == f]; y = y if len(y) else x; flag = str(y.iloc[0].centromeric_L1)
        out.append(None if flag == "na" else flag == "yes")               # v4.7: gene arrays (GarSat49, MacSat06) carry 'na' = not assessed -> unknown (pale ring)
    return out


def _cut(m, t, smin=SUP_MIN, minsize=MIN_SIZE):
    """maximal subtrees (preorder node ids) with patristic diameter <= t, UFBoot >= smin and >= minsize copies; the root (arbitrary trifurcation of the unrooted ML tree) is never a clade"""
    out = []; st = list(m["kids"][0])
    while st:
        i = st.pop(); d = m["desc"][i]
        if m["all5s"][i]: continue                                   # the 5S rRNA gene references are a reference group, never a satellite sub-clade
        if len(d) >= 2 and m["diam"][i] <= t and m["sup"][i] is not None and m["sup"][i] >= smin: out.append(i); continue
        st.extend(m["kids"][i])
    return sorted([i for i in out if len(m["desc"][i]) >= minsize])


def _keep(P1, P2):
    """fraction of the tips assigned in partition P1 (list of tip sets) that keep their clade partners in P2 (best Jaccard of their clade with a clade of P2 >= JACC)"""
    tot = sum(len(c) for c in P1)
    if not tot: return float("nan")
    k = 0
    for c in P1:
        if P2 and max(len(c & d) / len(c | d) for d in P2) >= JACC: k += len(c)
    return k / tot


def build_model(key):
    stem = LIN[key][0]; cf = CACHE / f"model_{stem}.pkl"
    if cf.exists():
        mo = pickle.load(open(cf, "rb"))
        if mo.get("ver") == MODEL_VER: return mo
    from Bio import Phylo
    tree = Phylo.read(str(PHY / stem / f"{stem}.treefile"), "newick")
    cl = list(tree.find_clades(order="preorder")); n = len(cl); idx = {id(c): i for i, c in enumerate(cl)}
    kids = [[idx[id(k)] for k in c.clades] for c in cl]; par = [-1] * n
    for i in range(n):
        for k in kids[i]: par[k] = i
    bl = np.array([(c.branch_length or 0.0) if par[i] >= 0 else 0.0 for i, c in enumerate(cl)]); dep = np.zeros(n)
    for i in range(1, n): dep[i] = dep[par[i]] + bl[i]
    tipi = [i for i in range(n) if not kids[i]]; tpos = {i: k for k, i in enumerate(tipi)}; mt = len(tipi); names = [cl[i].name for i in tipi]
    desc = [None] * n; diam = np.zeros(n); D = np.zeros((mt, mt)); dt = np.array([dep[i] for i in tipi])
    for v in range(n - 1, -1, -1):
        if not kids[v]: desc[v] = [tpos[v]]; continue
        desc[v] = [x for k in kids[v] for x in desc[k]]; diam[v] = max(diam[k] for k in kids[v])
        for a in range(len(kids[v])):
            for b in range(a + 1, len(kids[v])):
                da, db = desc[kids[v][a]], desc[kids[v][b]]
                blk = (dt[da] - dep[v])[:, None] + (dt[db] - dep[v])[None, :]; D[np.ix_(da, db)] = blk; D[np.ix_(db, da)] = blk.T; diam[v] = max(diam[v], blk.max())
    sup = [c.confidence if not c.is_terminal() else None for c in cl]
    is5 = np.array(["|5S_rRNA|" in x for x in names]); m = dict(kids=kids, desc=desc, diam=diam, sup=sup, all5s=[bool(is5[desc[i]].all()) for i in range(n)])
    # ---- sweep of the patristic-diameter cut
    parts = {}; nodes_at = {}
    for t in CUT_GRID: nodes_at[t] = _cut(m, t); parts[t] = [set(desc[i]) for i in nodes_at[t]]
    sig = {t: frozenset(frozenset(c) for c in parts[t]) for t in CUT_GRID}
    runs = []; cur = [CUT_GRID[0]]
    for t in CUT_GRID[1:]:
        if sig[t] == sig[cur[-1]]: cur.append(t)
        else: runs.append(cur); cur = [t]
    runs.append(cur)
    ok = lambda P: len(P) >= 4 and max(len(c) for c in P) <= 0.6 * sum(len(c) for c in P)             # a non-trivial partition: >= 4 clades, none holding > 60 % of the assigned copies
    best = max([r for r in runs if ok(parts[r[0]])], key=lambda r: (len(r), sum(len(c) for c in parts[r[0]]), r[0]))          # longest plateau of IDENTICAL clade sets; ties: more copies assigned, then coarser
    cut = float(best[(len(best) - 1) // 2]); Pstar = parts[cut]; nodes = nodes_at[cut]
    sweep = []
    for k, t in enumerate(CUT_GRID):
        P = parts[t]
        sweep.append(dict(cut=t, n_clades=len(P), tips_assigned=sum(len(c) for c in P), frac_assigned=round(sum(len(c) for c in P) / mt, 3), largest=max((len(c) for c in P), default=0),
                          keep_vs_chosen=round(_keep(Pstar, P), 3), keep_vs_previous=round(_keep(parts[CUT_GRID[k - 1]], P), 3) if k else "", keep_vs_next=round(_keep(P, parts[CUT_GRID[k + 1]]), 3) if k < len(CUT_GRID) - 1 else "",
                          same_set_as_chosen=int(sig[t] == sig[cut]), plateau=int(t in best)))
    supsweep = []
    for s_ in SUP_GRID:
        P = [set(desc[i]) for i in _cut(m, cut, smin=s_)]
        supsweep.append(dict(ufboot_min=s_, cut=cut, n_clades=len(P), tips_assigned=sum(len(c) for c in P), keep_vs_chosen=round(_keep(Pstar, P), 3)))
    # ---- clade table
    tip_cen = _tip_centromeric(names); order = sorted(nodes, key=lambda i: (-len(desc[i]), np.mean(desc[i])))
    cla = []
    for k, i in enumerate(order):
        tl = sorted(desc[i]); nn = [names[x] for x in tl]; sp = collections.Counter(x.split("|")[0] for x in nn); fam = collections.Counter(x.split("|")[1] for x in nn)
        grp = collections.OrderedDict()
        for x in tl:
            s_, f_, c_, _ = names[x].split("|"); g = grp.setdefault((s_, c_), dict(sp=s_, chrom=c_, n=0, cen=[]))
            g["n"] += 1
            if tip_cen[x] is not None: g["cen"].append(tip_cen[x])
        groups = [dict(sp=g["sp"], chrom=g["chrom"], n=g["n"], cen=(float(np.mean(g["cen"])) if g["cen"] else None)) for g in grp.values()]
        cross, exp_ = _nn_cross([names[x] for x in tl], D[np.ix_(tl, tl)])
        cla.append(dict(name=f"{PREFIX[key]}{k + 1}", node=i, tips=tl, n=len(tl), species=dict(sp), families=dict(fam), dom_family=fam.most_common(1)[0][0], groups=groups, n_chrom=len(groups),
                        n_cen_arrays=sum(1 for g in groups if g["cen"] is not None and g["cen"] >= 0.5), n_other_arrays=sum(1 for g in groups if g["cen"] is not None and g["cen"] < 0.5),
                        n_unknown_arrays=sum(1 for g in groups if g["cen"] is None), nn_cross_species=cross, nn_cross_expected_random=exp_,
                        root_ufboot=sup[i], diameter=round(float(diam[i]), 4), mean_pairwise=round(float(D[np.ix_(tl, tl)].sum() / (len(tl) * (len(tl) - 1))), 4)))
    pc = [t for t in CUT_GRID if 0.30 <= t <= 1.0]
    for c in cla:                                                       # persistence: share of cuts 0.30-1.00 at which a clade with Jaccard >= 0.8 to this one exists
        cs_ = set(c["tips"]); c["persistence"] = round(float(np.mean([max((len(cs_ & d_) / len(cs_ | d_) for d_ in parts[t]), default=0) >= JACC for t in pc])), 3)
    # ---- backbone nodes = nodes outside every clade that join >= 2 sub-trees containing a clade
    inside = set(); has = [False] * n
    for c in cla:
        st = [c["node"]]
        while st: i = st.pop(); inside.add(i); st.extend(kids[i])
    cn = {c["node"] for c in cla}
    for i in range(n - 1, -1, -1): has[i] = (i in cn) or any(has[k] for k in kids[i])
    backbone = [dict(node=i, ufboot=sup[i], n_clades_below=sum(1 for c in cla if set(desc[c["node"]]) <= set(desc[i])))
                for i in range(1, n) if i not in inside and len(kids[i]) >= 2 and sum(has[k] for k in kids[i]) >= 2]
    attach = None; five = None; t5 = [x for x in range(mt) if is5[x]]
    if t5:
        v = max((i for i in range(n) if set(t5) <= set(desc[i])), key=lambda i: len(_anc(par, i)))          # MRCA of all 5S rRNA genes
        grp = collections.OrderedDict()
        for x in t5:
            s_, f_, c_, _ = names[x].split("|"); g = grp.setdefault((s_, c_), dict(sp=s_, chrom=c_, n=0)); g["n"] += 1
        a = par[v]
        while a > 0 and not any(c["node"] in set(desc_nodes(kids, a)) for c in cla): a = par[a]          # lowest ancestor that contains a satellite clade
        below = [c["name"] for c in cla if c["node"] in set(desc_nodes(kids, a))]; chain = []; q = v
        while q != a: q = par[q]; chain.append(sup[q])
        five = dict(name="5S", node=v, tips=sorted(t5), n=len(t5), monophyletic=bool(len(desc[v]) == len(t5)), ufboot=sup[v], diameter=round(float(diam[v]), 4),
                    groups=[dict(sp=g["sp"], chrom=g["chrom"], n=g["n"], cen=None) for g in grp.values()], n_chrom=len(grp), families={"5S_rRNA": len(t5)}, dom_family="5S_rRNA",
                    species=dict(collections.Counter(names[x].split("|")[0] for x in t5)), long_branch=round(float(bl[v]), 4))
        attach = dict(node=a, ufboot=sup[a], clades_below=below, ufboot_chain=chain, n_tips=len(desc[a]), families=dict(collections.Counter(names[x].split("|")[1] for x in desc[a])))
    ix = [x for x in range(mt) if not is5[x]]; nn_global, nn_exp = _nn_cross([names[x] for x in ix], D[np.ix_(ix, ix)])
    inclade = {x for c in cla for x in c["tips"]}
    mo = dict(ver=MODEL_VER, key=key, stem=stem, names=names, n_tips=mt, clades=cla, five=five, cut=cut, plateau=(float(best[0]), float(best[-1])), sweep=sweep, supsweep=supsweep, backbone=backbone, attach=attach,
              unassigned=[names[x] for x in range(mt) if x not in inclade], nn_cross_global=nn_global, nn_cross_global_expected=nn_exp, tree_len=float(bl.sum()), median_pairwise=float(np.median(D[np.triu_indices(mt, 1)])))
    CACHE.mkdir(exist_ok=True); pickle.dump(mo, open(cf, "wb")); return mo


def _nn_cross(nm, Dm):
    """species-mixing test: for every copy, the nearest copy (patristic) from a DIFFERENT chromosome array (same-array neighbours are trivially conspecific); returns (fraction whose nearest
    other-array copy is of the other species, the value expected if species were unrelated to tree position = mean share of other-species copies among the other-array copies); (None, None) if undefined"""
    sp = np.array([x.split("|")[0] for x in nm]); ar = np.array([x.split("|")[0] + x.split("|")[2] for x in nm]); hit = []; ex = []
    for i in range(len(nm)):
        o = ar != ar[i]
        if not o.any(): continue
        j = np.where(o)[0][np.argmin(Dm[i, o])]; hit.append(sp[j] != sp[i]); ex.append(np.mean(sp[o] != sp[i]))
    return (round(float(np.mean(hit)), 3), round(float(np.mean(ex)), 3)) if hit else (None, None)


def desc_nodes(kids, i):
    out = [i]
    for k in kids[i]: out.extend(desc_nodes(kids, k))
    return out


def _anc(par, i):
    out = []
    while i >= 0: out.append(i); i = par[i]
    return out


# ---------------------------------------------------------------- drawing
def mmx(x): return x * 25.4 / 72      # points -> mm


def text_w(ax, s, fs, weight="normal", style="normal"):
    """rendered width of a string in mm (measured with the Agg renderer)"""
    fig = ax.figure; t = ax.text(0, 0, s, fontsize=fs, fontweight=weight, fontstyle=style)
    bb = t.get_window_extent(fig.canvas.get_renderer()); t.remove(); return bb.width / fig.dpi * 25.4


def unpile(P, mask, dmin=1.05, iters=300):
    """beeswarm-like un-piling of the 5S rRNA gene markers: tips closer than dmin mm push each other apart (weak spring back to the true tip position, so the branch end stays readable)"""
    P = P.copy(); idx = np.where(mask)[0]; P0 = P[idx].copy(); Q = P0.copy()
    for _ in range(iters):
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                v = Q[a] - Q[b]; dd = np.hypot(*v)
                if dd < dmin:
                    v = v / dd * 0.5 if dd > 1e-9 else np.array([np.cos(a * 2.4), np.sin(a * 2.4)]) * 0.5; push = (dmin - dd) / 2; Q[a] += v * push / 0.5 * 0.5; Q[b] -= v * push / 0.5 * 0.5
        Q += 0.02 * (P0 - Q)
    P[idx] = Q; return P


def draw_tree(ax, d, box, scale=None):
    """draw one tree into the mm-axes `ax` inside box=(x0, y0, w, h) (mm, y down).  scale = mm per substitution/site (None = fit the box)."""
    x0, y0, w, h = box; xy = d["xy"]; segs = d["segs"]
    allp = segs.reshape(-1, 2); lo, hi = allp.min(0), allp.max(0)
    fit = min(w / (hi[0] - lo[0]), h / (hi[1] - lo[1]))
    sc = fit * 0.98 if scale is None else scale
    cx, cy = x0 + w / 2 - sc * (lo[0] + hi[0]) / 2, y0 + h / 2 + sc * (lo[1] + hi[1]) / 2
    T = lambda p: np.column_stack([cx + sc * p[..., 0], cy - sc * p[..., 1]])
    ax.add_collection(LineCollection([T(s) for s in segs], colors="#B9BEC2", linewidths=0.3, zorder=1, capstyle="round"))
    info, fams = families(d); P = T(xy); P = unpile(P, [t[1] == "5S_rRNA" for t in info])
    for f in fams:
        for sp in ("Gar", "Mac"):
            sel = [i for i, t in enumerate(info) if t[1] == f and t[0] == sp]
            if not sel: continue
            m = SPECIES_MARK["5S" if f == "5S_rRNA" else sp]
            ax.scatter(P[sel, 0], P[sel, 1], s=5 if f == "5S_rRNA" else 3.4, marker=m, color=fcol(f), edgecolor="white", linewidth=0.12,
                       zorder=3 if f == "5S_rRNA" else 2, rasterized=True)
    return sc, T, P, info, fams, [T(s) for s in segs]


def _tint(c, f):
    """blend colour c with white (f = share of white)"""
    import matplotlib.colors as mc
    r = np.array(mc.to_rgb(c)); return tuple(r + (1 - r) * f)


def ring_geometry(P, cla, centre, margin=2.6, thick=0.6):
    """per clade: angular window seen from the tree centre and ring radius = beyond every tip inside the window (own or foreign); overlapping windows are stacked in steps of `thick`+0.3 mm"""
    ang_all = np.arctan2(P[:, 1] - centre[1], P[:, 0] - centre[0]); rad_all = np.hypot(*(P - centre).T); geo = []
    for c in cla:
        a = np.sort(ang_all[c["tips"]]); gaps = np.diff(np.r_[a, a[0] + 2 * np.pi]); j = int(np.argmax(gaps)); start = a[(j + 1) % len(a)]; span = 2 * np.pi - gaps[j]
        r_own = rad_all[c["tips"]].max(); pad = 0.45 / max(r_own, 1.0); start -= pad; span += 2 * pad
        rel = (ang_all - start) % (2 * np.pi); inwin = rel <= span
        geo.append(dict(start=start, span=span, r=max(r_own, rad_all[inwin].max()) + margin + (0.2 if c["dom_family"] == "5S_rRNA" else 0.0), r_own=r_own))
    for c, g in zip(cla, geo):                         # buried clade (ring would float > BURIED mm beyond its own copies): ring it locally, around its own copies, when no foreign copy falls inside
        if g["r"] - g["r_own"] > BURIED:
            q = P[c["tips"]]; c0 = q.mean(0); r0 = np.hypot(*(q - c0).T).max() + 0.8
            foreign = np.delete(P, c["tips"], axis=0)
            if not np.any(np.hypot(*(foreign - c0).T) < r0 + 0.4): g.update(local=True, c=c0, r=r0, start=0.0, span=2 * np.pi, out=(c0 - centre) / max(np.hypot(*(c0 - centre)), 1e-6))
    for _ in range(40):
        moved = False
        for i in range(len(geo)):
            for j in range(i + 1, len(geo)):
                gi, gj = geo[i], geo[j]
                if gi.get("local") or gj.get("local"): continue
                ov = ((gj["start"] - gi["start"]) % (2 * np.pi) <= gi["span"]) or ((gi["start"] - gj["start"]) % (2 * np.pi) <= gj["span"])
                if ov and abs(gi["r"] - gj["r"]) < thick + 0.2:
                    (gj if gj["r"] >= gi["r"] else gi)["r"] += (thick + 0.2) - abs(gi["r"] - gj["r"]); moved = True
        if not moved: break
    return geo, thick


def draw_rings(ax, P, cla, centre, geo, thick):
    """level 2: ring segments, one per chromosome array of the clade (width = copies), ordered by the mean angle of their copies; hue = species, dark = centromeric array, pale = other, open = 5S genes"""
    samples = []
    for c, g in zip(cla, geo):
        c0 = g.get("c", centre); ang_all = np.arctan2(P[:, 1] - c0[1], P[:, 0] - c0[0]); sel = {}
        for x in c["tips"]:
            s_, f_, ch, _ = c["_names"][x].split("|"); sel.setdefault((s_, ch), []).append((ang_all[x] - g["start"]) % (2 * np.pi))
        grp = sorted(c["groups"], key=lambda q: np.mean(sel[(q["sp"], q["chrom"])])); tot = sum(q["n"] for q in grp); a0 = 0.0
        for q in grp:
            w = g["span"] * q["n"] / tot; th1 = np.degrees(g["start"] + a0); th2 = np.degrees(g["start"] + a0 + w); a0 += w
            col = ST.SPECIES[q["sp"]]; five = c["dom_family"] == "5S_rRNA"
            if five: fc, ec, lw = "none", col, 0.35
            elif q["cen"] is None: fc, ec, lw = _tint(col, 0.8), "white", 0.15
            elif q["cen"] >= 0.5: fc, ec, lw = col, "white", 0.15
            else: fc, ec, lw = _tint(col, 0.62), "white", 0.15
            ax.add_patch(Wedge(c0, g["r"], th1, th2, width=thick, facecolor=fc, edgecolor=ec, linewidth=lw, zorder=4))
        t = np.linspace(g["start"], g["start"] + g["span"], max(3, int(g["span"] * g["r"] / 0.4)))
        samples.append(np.column_stack([c0[0] + (g["r"] - thick / 2) * np.cos(t), c0[1] + (g["r"] - thick / 2) * np.sin(t)]))
    return samples


class Placer:
    """free-space label placement on a 0.25 mm occupancy grid inside `cell` (x0, y0, w, h mm): a label box goes to the free position closest to its anchor"""
    RES = 0.25

    def __init__(self, cell):
        self.x0, self.y0, w, h = cell; self.nx, self.ny = int(w / self.RES), int(h / self.RES); self.O = np.zeros((self.ny, self.nx), bool)

    def mark_points(self, pts, rad):
        r = max(0, int(round(rad / self.RES)))
        for x, y in pts:
            i = int((y - self.y0) / self.RES); j = int((x - self.x0) / self.RES)
            if -r <= i < self.ny + r and -r <= j < self.nx + r: self.O[max(0, i - r):i + r + 1, max(0, j - r):j + r + 1] = True

    def mark_box(self, b, pad=0.6):
        i0 = int((b[1] - pad - self.y0) / self.RES); i1 = int(np.ceil((b[3] + pad - self.y0) / self.RES)); j0 = int((b[0] - pad - self.x0) / self.RES); j1 = int(np.ceil((b[2] + pad - self.x0) / self.RES))
        self.O[max(0, i0):max(0, i1), max(0, j0):max(0, j1)] = True

    def place_ray(self, w, h, c0, th0, r0, pad=0.8, dmax=55, dstep=1.5, kmax=10, kstep=0.6):
        """angular de-collision, leader-free: candidate boxes lie on rays from the ring centre at angles th0 +- delta and radii r0 + k*kstep; the box touches the ray point with its nearest side
        (right of it on the right half, left of it on the left half, above / below it at the poles); returns the free candidate with the smallest angular shift (then radial lift), or None"""
        best = None
        for k in range(kmax + 1):
            for dd in np.arange(0, dmax + 1e-9, dstep):
                for sg in ((1,) if dd == 0 else (1, -1)):
                    th = th0 + np.radians(sg * dd); R = r0 + k * kstep; ax_, ay_ = c0[0] + R * np.cos(th), c0[1] + R * np.sin(th); cx, cy = np.cos(th), np.sin(th)
                    x = ax_ if cx > 0.35 else ax_ - w if cx < -0.35 else ax_ - w / 2; y = ay_ if cy > 0.35 else ay_ - h if cy < -0.35 else ay_ - h / 2
                    j0, j1 = int((x - pad - self.x0) / self.RES), int(np.ceil((x + w + pad - self.x0) / self.RES)); i0, i1 = int((y - pad - self.y0) / self.RES), int(np.ceil((y + h + pad - self.y0) / self.RES))
                    if j0 < 0 or i0 < 0 or j1 > self.nx or i1 > self.ny or self.O[i0:i1, j0:j1].any(): continue
                    cost = dd + 2.0 * k
                    if best is None or cost < best[0]: best = (cost, x, y)
            if best is not None and best[0] <= 2.0 * k: break
        return None if best is None else (best[1], best[2])

    def place(self, w, h, anchor, own, others=(), pad=0.3):
        """returns (x, y) lower-left of the box (mm; y is the TOP of the box in y-down mm) closest to `anchor`; own / others = arc sample points of this clade / of the neighbouring clades"""
        wc, hc = int(np.ceil((w + 2 * pad) / self.RES)), int(np.ceil((h + 2 * pad) / self.RES))
        S = np.zeros((self.ny + 1, self.nx + 1)); S[1:, 1:] = np.cumsum(np.cumsum(self.O, 0), 1)
        ni, nj = self.ny - hc + 1, self.nx - wc + 1
        if ni <= 0 or nj <= 0: return None
        I, J = np.meshgrid(np.arange(ni), np.arange(nj), indexing="ij")
        free = (S[I + hc, J + wc] - S[I, J + wc] - S[I + hc, J] + S[I, J]) == 0
        if not free.any(): return None
        ii, jj = I[free], J[free]; cx = self.x0 + (jj + wc / 2) * self.RES; cy = self.y0 + (ii + hc / 2) * self.RES
        c = np.column_stack([cx, cy]); cost = np.hypot(*(c - np.asarray(anchor)).T)
        do = np.min(np.hypot(c[:, None, 0] - own[None, :, 0], c[:, None, 1] - own[None, :, 1]), axis=1)
        cost = 0.6 * cost + 1.6 * do
        if len(others):
            dn = np.min(np.hypot(c[:, None, 0] - others[None, :, 0], c[:, None, 1] - others[None, :, 1]), axis=1); cost = cost + 8.0 * (dn < do)
        k = int(np.argmin(cost)); x = self.x0 + jj[k] * self.RES + pad; y = self.y0 + ii[k] * self.RES + pad
        return x, y


def label_colour(c):
    """name colour = the clade's dominant family colour, darkened for legibility at 6 pt (ties the label to the tips inside the ring)"""
    if c["dom_family"] == "5S_rRNA": return ST.NAVY
    import matplotlib.colors as mc
    return tuple(np.array(mc.to_rgb(ST.family_colour(c["dom_family"]))) * 0.72)


def clade_label(c):
    if c["name"] == "5S": return "5S rRNA genes", str(c["n"])
    """name + copies (multi-chromosome clades) or name + chromosome (single-chromosome clades: all 12 sampled copies of that array)"""
    return c["name"], (c["groups"][0]["chrom"] if c["n_chrom"] == 1 else str(c["n"]))


def draw_lineage(ax, key, cell, tree_box, scale=None, fs=6):
    """one lineage: tree (level 0), clade rings + labels (levels 1-2), backbone support dots + 5S attachment (level 3).  cell = (x0, y0, w, h) mm available to tree + rings + labels"""
    d = load_tree(key); mo = build_model(key); cla = mo["clades"] + ([mo["five"]] if mo["five"] else [])          # + the 5S rRNA gene reference group
    assert mo["names"] == d["names"]
    for c in cla: c["_names"] = mo["names"]
    sc, T, P, info, fams, sg = draw_tree(ax, d, tree_box, scale)
    lo, hi = P.min(0), P.max(0); centre = (lo + hi) / 2
    geo, thick = ring_geometry(P, cla, centre)
    samples = draw_rings(ax, P, cla, centre, geo, thick)
    N = T(node_xy(d))
    # level 3: UFBoot dots on the backbone
    bb = mo["backbone"]; hi_ = [N[b["node"]] for b in bb if b["ufboot"] is not None and b["ufboot"] >= 95]; mid_ = [N[b["node"]] for b in bb if b["ufboot"] is not None and 70 <= b["ufboot"] < 95]
    if hi_: ax.scatter(*np.array(hi_).T, s=4.2, marker="o", color=ST.NAVY, edgecolor="white", linewidth=0.15, zorder=5)
    if mid_: ax.scatter(*np.array(mid_).T, s=4.2, marker="o", color="white", edgecolor=ST.NAVY, linewidth=0.3, zorder=5)
    pl = Placer(cell); allpts = [P]
    for a, b in sg:
        n_ = max(2, int(np.hypot(*(b - a)) / 0.35)); t = np.linspace(0, 1, n_)[:, None]; allpts.append(a + (b - a) * t)
    pl.mark_points(np.vstack(allpts), 0.35)
    for s in samples: pl.mark_points(s, 0.7)
    pl.mark_points(N[[b["node"] for b in bb]], 0.5)
    texts = {}; items = []
    for k, (c, g, s) in enumerate(zip(cla, geo, samples)):
        mid = g["start"] + g["span"] / 2
        anc = g["c"] + (g["r"] + 1.4) * g["out"] if g.get("local") else centre + (g["r"] + 1.4) * np.array([np.cos(mid), np.sin(mid)])
        ray = (g["c"], float(np.arctan2(g["out"][1], g["out"][0])), g["r"] + 0.9) if g.get("local") else (centre, mid, g["r"] + 0.9)
        nm, tail = clade_label(c); wn, wt = text_w(ax, nm, fs, "normal" if nm.startswith("5S") else "bold"), text_w(ax, " " + tail, fs); items.append((c, nm, tail, wn, wt, anc, s, g["span"], ray))
    at = mo["attach"]
    for c, nm, tail, wn, wt, anc, s, span, ray in sorted(items, key=lambda q: -q[7]):
        others = np.vstack([q[6] for q in items if q[0] is not c]); h = fs * 0.3528 * 1.0
        pos = None
        for padq in (0.9, 0.7, 0.5, 0.3):
            pos = pl.place_ray(wn + wt, h, ray[0], ray[1], ray[2], pad=padq, dmax=100, kmax=14)
            if pos is not None: break
        for padq in (0.8, 0.6, 0.45, 0.3):
            if pos is not None: break
            pos = pl.place(wn + wt, h, anc, s, others, pad=padq)
        if pos is None: raise RuntimeError(f"no free label position for {key} {nm}")
        x, y = pos; pl.mark_box((x, y, x + wn + wt, y + h)); texts[nm] = (round(x, 2), round(y, 2))
        ax.text(x, y + h / 2, nm, fontsize=fs, fontweight="normal" if nm.startswith("5S") else "bold", color=label_colour(c), ha="left", va="center", zorder=7)
        ax.text(x + wn, y + h / 2, " " + tail, fontsize=fs, color=ST.MUTED if True else ST.NAVY, ha="left", va="center", zorder=7)
    if at is not None:
        p = N[at["node"]]; ax.scatter([p[0]], [p[1]], s=15, marker="o", facecolor="none", edgecolor=GOLD5S, linewidth=0.6, zorder=6)
        txt = f"{at['ufboot']:g}"; wtx = text_w(ax, txt, fs); h = fs * 0.3528
        pos = pl.place(wtx, h, p + np.array([0.0, 0.0]), np.array([p]), pad=0.2)
        if pos is not None:
            x, y = pos; pl.mark_box((x, y, x + wtx, y + h)); ax.text(x, y + h / 2, txt, fontsize=fs, color=GOLD_DARK, ha="left", va="center", zorder=7); texts["5S_attach"] = (round(x, 2), round(y, 2))
    return dict(sc=sc, P=P, info=info, fams=fams, model=mo, labels=texts, centre=centre, geo=geo)


def family_key(ax, x, y, fams, font=6, w=None):
    """one coloured name per family (bold, family colour), wrapped into rows of at most `w` mm; 5S rDNA gene arrays (GarSat49, MacSat06) and the 5S-derived satellite are tagged"""
    xx = x; yy = y; names = [f for f in fams if f != "5S_rRNA"]
    for f in names:
        tag = " 5S rDNA" if ST.is_gene_array(f) else ""            # GarSat49 / MacSat06 = 5S rDNA gene arrays; GarSat04 (5S-derived satellite) is named in the caption
        wt = text_w(ax, f, font, "bold"); wtag = text_w(ax, tag, font) if tag else 0.0
        if w is not None and xx + wt + wtag > x + w: xx = x; yy += 2.7
        ax.text(xx, yy, f, fontsize=font, fontweight="bold", color=ST.family_colour(f), ha="left", va="center"); xx += wt
        if tag: ax.text(xx, yy, tag, fontsize=font, color=ST.MUTED, ha="left", va="center"); xx += wtag
        xx += 1.6
    return yy


def key_row(ax, x, y, font=6):
    """marker key: Gar circle, Mac square, 5S diamond (grey markers: the colour carries the family)"""
    for lab, m, c in (("C. gariepinus", "o", ST.MUTED), ("C. macrocephalus", "s", ST.MUTED), ("5S rRNA gene", "D", GOLD5S)):
        ax.scatter([x + 0.9], [y], s=7 if m != "D" else 8, marker=m, color=c, edgecolor="white", linewidth=0.15, zorder=3)
        st = "italic" if lab[0] == "C" else "normal"; ax.text(x + 2.2, y, lab, fontsize=font, ha="left", va="center", style=st)
        x += 2.2 + text_w(ax, lab, font, style=st) + 1.6
    return x


def support_key(ax, x, y, font=6):
    ax.scatter([x + 0.7], [y], s=4.2, color=ST.NAVY, edgecolor="white", linewidth=0.15, zorder=3); ax.text(x + 1.7, y, "UFBoot ≥ 95", fontsize=font, ha="left", va="center")
    x2 = x + 1.7 + text_w(ax, "UFBoot ≥ 95", font) + 1.4
    ax.scatter([x2 + 0.7], [y], s=4.2, color="white", edgecolor=ST.NAVY, linewidth=0.3, zorder=3); ax.text(x2 + 1.7, y, "70–94", fontsize=font, ha="left", va="center")
    return x2 + 1.7 + text_w(ax, "70–94", font)


def ring_key(ax, x, y, font=6):
    """key of the clade ring: one segment per chromosome array; hue = species, dark = centromeric array, pale = other"""
    from matplotlib.patches import Rectangle
    ax.text(x, y, "Ring: array per chromosome", fontsize=font, ha="left", va="center"); x += text_w(ax, "Ring: array per chromosome", font) + 1.4
    for lab, col, f in (("Gar", ST.SPECIES["Gar"], 0), ("Mac", ST.SPECIES["Mac"], 0)):
        ax.add_patch(Rectangle((x, y - 0.6), 1.9, 1.2, facecolor=col, edgecolor="none")); ax.text(x + 2.4, y, lab, fontsize=font, ha="left", va="center", style="normal")
        x += 2.4 + text_w(ax, lab, font) + 1.2
    ax.add_patch(Rectangle((x, y - 0.6), 1.9, 1.2, facecolor="#555B61", edgecolor="none")); ax.add_patch(Rectangle((x + 2.0, y - 0.6), 1.9, 1.2, facecolor=_tint("#555B61", 0.62), edgecolor="none"))
    ax.text(x + 4.4, y, "centromeric / other", fontsize=font, ha="left", va="center")
    return x + 4.4 + text_w(ax, "centromeric / other", font)


def scalebar(ax, x, y, mm_per_unit, val=0.05, font=6, label=True):
    L = val * mm_per_unit
    ax.plot([x, x + L], [y, y], color=ST.NAVY, lw=0.8, solid_capstyle="butt", zorder=5)
    if label: ax.text(x + L + 0.8, y, f"{val:g} substitutions/site", fontsize=font, ha="left", va="center")
    return L


def numbers_rows(panel, key, d):
    """family rows (tips per family and species) + clade rows (one per sub-clade) + backbone / 5S rows"""
    info, fams = families(d); mo = build_model(key); rows = []
    for f in fams:
        for sp in ("Gar", "Mac"):
            n = sum(1 for t in info if t[1] == f and t[0] == sp)
            if n: rows.append(dict(panel=panel, row_type="family", tree=TREE_NAME[key], name=f, species=sp, tips=n, tips_in_tree=len(info), tree_file=f"{d['stem']}.treefile",
                                   tree_total_branch_length=round(d["tree_len"], 3), element_class_v47=("gene_array" if ST.is_gene_array(f) else "satellite"), seq_cluster_v47=ST.family_cluster(f)))
    for c in mo["clades"]:
        rows.append(dict(panel=panel, row_type="clade", tree=TREE_NAME[key], name=c["name"], species="+".join(f"{k}:{v}" for k, v in sorted(c["species"].items())), tips=c["n"], tips_in_tree=len(info),
                         tree_file=f"{d['stem']}.treefile", tree_total_branch_length=round(d["tree_len"], 3),
                         families="+".join(f"{k}:{v}" for k, v in sorted(c["families"].items())), n_chromosome_arrays=c["n_chrom"],
                         chromosomes=";".join(f"{g['sp']}_{g['chrom']}:{g['n']}" for g in c["groups"]), centromeric_arrays=c["n_cen_arrays"], other_arrays=c["n_other_arrays"], unknown_arrays=c["n_unknown_arrays"],
                         nn_cross_species=c["nn_cross_species"], nn_cross_expected_random=c["nn_cross_expected_random"], root_ufboot=c["root_ufboot"], diameter=c["diameter"], mean_pairwise=c["mean_pairwise"], persistence_cuts_030_100=c["persistence"],
                         cut=mo["cut"], ufboot_min=SUP_MIN, min_copies=MIN_SIZE))
    for b in mo["backbone"]:
        rows.append(dict(panel=panel, row_type="backbone_node", tree=TREE_NAME[key], name=f"node{b['node']}", ufboot_node=b["ufboot"], clades_below=b["n_clades_below"], tree_file=f"{d['stem']}.treefile"))
    if mo["attach"]:
        a = mo["attach"]; f5 = mo["five"]
        rows.append(dict(panel=panel, row_type="reference_group", tree=TREE_NAME[key], name="5S rRNA genes", species="+".join(f"{k}:{v}" for k, v in sorted(f5["species"].items())), tips=f5["n"], tips_in_tree=len(info),
                         chromosomes=";".join(f"{g['sp']}_{g['chrom']}:{g['n']}" for g in f5["groups"]), n_chromosome_arrays=f5["n_chrom"], root_ufboot=f5["ufboot"], diameter=f5["diameter"], tree_file=f"{d['stem']}.treefile",
                         monophyletic=int(f5["monophyletic"]), long_branch=f5["long_branch"]))
        rows.append(dict(panel=panel, row_type="5S_attachment_node", tree=TREE_NAME[key], name="5S group + " + "+".join(a["clades_below"]), ufboot_node=a["ufboot"], n_tips_below=a["n_tips"],
                         families="+".join(f"{k}:{v}" for k, v in sorted(a["families"].items())), ufboot_chain="/".join(f"{x:g}" for x in a["ufboot_chain"]), clades_below=len(a["clades_below"]), tree_file=f"{d['stem']}.treefile"))
    rows.append(dict(panel=panel, row_type="summary", tree=TREE_NAME[key], name="clade_set", tips=mo["n_tips"], tips_assigned=sum(c["n"] for c in mo["clades"]), n_clades=len(mo["clades"]), cut=mo["cut"],
                     plateau=f"{mo['plateau'][0]}-{mo['plateau'][1]}", nn_cross_species=mo["nn_cross_global"], nn_cross_expected_random=mo["nn_cross_global_expected"], tree_file=f"{d['stem']}.treefile", tree_total_branch_length=round(d["tree_len"], 3)))
    return rows


# layout of the panel (mm)
TREE_MX, TREE_MY = {"B": 6.5, "A": 8.5}, {"B": 8.8, "A": 7.5}
COLS = {"B": (0.5, 40.0), "A": (40.5, 51.0)}       # (x0, width) mm of the two cells: the centromeric tree carries 19 clade labels and gets the wider cell
CELL_Y0, CELL_Y1, FAM_Y, SCALE_Y, K1_Y, K2_Y = 7.4, 41.6, 43.9, 49.2, 52.0, 54.7


def draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None):
    """two trees, each fitted to its own half of the panel, each with its own scale bar; levels 1-3 per tree"""
    ax = ST.ax_mm(fig, left_mm, top_mm, w_mm, h_mm); ax.set_xlim(0, w_mm); ax.set_ylim(h_mm, 0); ax.axis("off")
    if letter: ST.letter(fig, left_mm, top_mm, letter)
    rows = []; scales = {}; labs = {}; models = {}
    for k, key in enumerate(("B", "A")):
        x0, colw = COLS[key]
        ax.text(x0 + 0.3, 5.6, LIN[key][1], fontsize=7, fontweight="bold", color=GROUPCOL[key], ha="left", va="center")
        r = draw_lineage(ax, key, (x0, CELL_Y0, colw, CELL_Y1 - CELL_Y0), (x0 + TREE_MX[key], CELL_Y0 + TREE_MY[key], colw - 2 * TREE_MX[key], CELL_Y1 - CELL_Y0 - 2 * TREE_MY[key]))
        scales[key] = r["sc"]; labs[key] = r["labels"]; models[key] = r["model"]
        family_key(ax, x0 + 0.3, FAM_Y, r["fams"], w=colw - 0.6)
        scalebar(ax, x0 + 0.8, SCALE_Y, r["sc"], label=False); ax.text(x0 + 0.8 + 0.05 * r["sc"] + 0.8, SCALE_Y, "0.05 substitutions/site", fontsize=6, ha="left", va="center")
        rows += numbers_rows("a", key, load_tree(key))
    xe = key_row(ax, 0.8, K1_Y); support_key(ax, xe + 1.0, K1_Y); ring_key(ax, 0.8, K2_Y)
    return dict(rows=rows, mm_per_unit=scales, label_positions_mm=labs, models=models)


def write_tables(num, stem="fig1_panel_a"):
    keys = []
    for r in num["rows"]:
        for k in r:
            if k not in keys: keys.append(k)
    with open(HERE / f"{stem}_numbers.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", restval=""); w.writeheader(); w.writerows(num["rows"])
    sw = []
    for key, mo in num["models"].items():
        for r in mo["sweep"]: sw.append(dict(tree=TREE_NAME[key], sweep="patristic_cut", **r))
        for r in mo["supsweep"]: sw.append(dict(tree=TREE_NAME[key], sweep="ufboot_min_at_chosen_cut", **r))
    keys = []
    for r in sw:
        for k in r:
            if k not in keys: keys.append(k)
    with open(HERE / f"{stem}_sweep.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", restval=""); w.writeheader(); w.writerows(sw)


def write_caption(num, stem="fig1_panel_a"):
    M = num["models"]; B, A_ = M["B"], M["A"]
    def sw(mo): return ", ".join(f"{r['cut']:g}:{r['n_clades']}" for r in mo["sweep"])
    def keep(mo, lo, hi): return [r["keep_vs_chosen"] for r in mo["sweep"] if lo - 1e-9 <= r["cut"] <= hi + 1e-9]
    def pers(mo): return sum(1 for c in mo["clades"] if c["persistence"] >= 0.5), len(mo["clades"])
    def sup(mo): return ", ".join(f">={r['ufboot_min']}: {r['n_clades']}" for r in mo["supsweep"])
    def npl(mo): return len([r for r in mo["sweep"] if r["plateau"]])
    def nat(mo, t): return next(r["n_clades"] for r in mo["sweep"] if abs(r["cut"] - t) < 1e-9)
    mixed = [c["name"] for mo in (B, A_) for c in mo["clades"] if len(c["species"]) > 1]
    a5 = B["attach"]; f5 = B["five"]; cen_arr = sum(c["n_cen_arrays"] for c in A_["clades"]); oth_arr = sum(c["n_other_arrays"] for c in A_["clades"])
    only_other = [c["name"] for c in A_["clades"] if c["n_cen_arrays"] == 0 and c["n_other_arrays"] > 0]
    unk = sum(c["n_unknown_arrays"] for c in A_["clades"])
    miss = "; ".join(f"{sp} {', '.join(c.split('_')[-1] for c in v)}" for sp, v in ST.MISSING_CHROMS.items())
    txt = f"""# Figure 1a caption stub (evidence level L1: sequence) - SatChain v4.7 wording, 2026-11-01

**Trees of the two satellite groups, with four levels of resolution.** Maximum-likelihood trees of the v4.5 monomer copies (IQ-TREE 2.3.6, GTR+G, 1,000 ultrafast bootstraps; treefiles
`share_LITE_v4.5_2026-10-15/phylogeny/{{lineageB_5S_derived,lineageA}}/*.treefile`; the file stems are the v4.5 input names, the word "lineage" is no longer used):
the 5S-related tree, {B['n_tips']} copies (GarSat04 = the 5S-derived satellite, SatChain v4.7 sequence cluster C2.2; the 5S rDNA gene arrays GarSat49 and MacSat06, which are gene arrays
and not satellite families; + 16 5S rRNA gene references); the centromeric tree, {A_['n_tips']} copies (GarSat03, GarSat09, GarSat10, MacSat03, MacSat162 = the GarSat03 / MacSat03 group,
the RTE-carrying centromeric families, sequence cluster C2.1); up to 12 sampled copies per chromosome array. Unrooted trees, equal-daylight layout, one scale bar each (substitutions/site);
the two groups are not joined (no joint alignment; between-group consensus identity 53-58 %). Tips: fill = family, circle = *C. gariepinus*, square = *C. macrocephalus*, diamond = 5S rRNA gene.
(1) **Sub-clades** R1-R{len(A_['clades'])} (centromeric tree) and S1-S{len(B['clades'])} (5S-related tree), numbered by size: maximal sub-trees with patristic diameter <= {B['cut']:g} (S) / {A_['cut']:g} (R)
substitutions/site, UFBoot >= {SUP_MIN} and >= {MIN_SIZE} copies (the 5S reference group is never a clade); {sum(c['n'] for c in B['clades'])}/{B['n_tips']} (S) and {sum(c['n'] for c in A_['clades'])}/{A_['n_tips']} (R) copies fall in a clade; label = name + copies,
or name + chromosome for single-array clades (all sampled copies of that array).
(2) **Chromosome origin**: ring around each clade, one segment per chromosome array (segment width = copies); hue = species (purple *C. gariepinus*, green *C. macrocephalus*); dark = array flagged
centromeric in the v4.7 runs of record (`centromeric_L1`), pale = other arrays, open = 5S rRNA gene loci (not arrays); two hues in one ring = clade with both species. (3) **Support**: UFBoot on the
backbone nodes (nodes joining >= 2 clades): filled >= 95, open 70-94, none < 70; the gold ring marks where the 5S rRNA gene group joins the 5S-related clades (UFBoot {a5['ufboot']:g}).

## Clade rule and sweep (parameter stability; table: `fig1_panel_a_sweep.tsv`)
Cut swept 0.10-1.20 in steps of 0.05 (23 values), clades per cut, S: {sw(B)}. R: {sw(A_)}. The drawn cut is the middle of the longest run of IDENTICAL clade sets (non-trivial partitions only: >= 4 clades,
none holding > 60 % of assigned copies; ties: more copies assigned, then coarser): S {B['cut']:g} (plateau {B['plateau'][0]:g}-{B['plateau'][1]:g}, {npl(B)} consecutive cuts, one clade set), R {A_['cut']:g}
(plateau {A_['plateau'][0]:g}-{A_['plateau'][1]:g}, {npl(A_)} cuts). Stability = share of assigned copies that keep their clade partners (best Jaccard >= {JACC} with a clade of the other cut):
S {min(keep(B, 0.45, 1.1)):.2f}-{max(keep(B, 0.45, 1.1)):.2f} over cuts 0.45-1.10; R {min(keep(A_, 0.65, 0.75)):.2f}-{max(keep(A_, 0.65, 0.75)):.2f} at the plateau, {keep(A_, 0.6, 0.6)[0]:.2f} at cut 0.60 and {keep(A_, 0.8, 0.8)[0]:.2f} at 0.80: the 5S-related clade set is
robust, whereas the {len(A_['clades'])}-clade level of the centromeric tree is a narrow plateau (cut 0.60 gives {nat(A_, 0.6)} clades, 0.80 gives {nat(A_, 0.8)}).
Clades present (Jaccard >= {JACC}) at >= half of the cuts 0.30-1.00: S {pers(B)[0]}/{pers(B)[1]}, R {pers(A_)[0]}/{pers(A_)[1]} (column `persistence_cuts_030_100` of the numbers table); the cut-dependent clades of the centromeric tree are {', '.join(c['name'] for c in A_['clades'] if c['persistence'] < 0.5)} (large multi-chromosome clades that merge / split with the cut).
UFBoot minimum swept at the drawn cut (clades): S {sup(B)}; R {sup(A_)}. The clade set of the centromeric tree depends on the support threshold (95 = conventional); the trees carry UFBoot only (IQ-TREE was run with -B 1000, no SH-aLRT).

## What the levels show
* Species mixing (share of copies whose nearest patristic neighbour on ANOTHER chromosome array belongs to the other species; `nn_cross_species`): over all satellite copies {B['nn_cross_global']:.3f} (S) and {A_['nn_cross_global']:.3f} (R), against {B['nn_cross_global_expected']:.3f} / {A_['nn_cross_global_expected']:.3f} expected if species were unrelated to tree position.
  {len(mixed)} clade(s) contain both species ({', '.join(mixed)}; share inside the clade {', '.join(str(c['nn_cross_species']) for mo in (B, A_) for c in mo['clades'] if len(c['species']) > 1)}, expected {', '.join(str(c['nn_cross_expected_random']) for mo in (B, A_) for c in mo['clades'] if len(c['species']) > 1)}): the two species form separate sub-groups, they do not interleave at copy level (species-specific expansion rather than a shared pool of monomers).
* Chromosome origin: the centromeric-tree clades hold {cen_arr} centromeric and {oth_arr} other arrays; {unk} array(s) without a covering v4.7 run (the v4.7 tables have no rows for {miss}, so copies on those chromosomes are drawn pale = unknown); clades made only of non-centromeric arrays: {', '.join(only_other)}.
* 5S: the {f5['n']} 5S rRNA genes ({f5['species'].get('Gar', 0)} Gar on LG04, {f5['species'].get('Mac', 0)} Mac on LG01 + LG20) form one clade (UFBoot {f5['ufboot']:g}, monophyletic = {bool(f5['monophyletic'])}) on a long branch ({f5['long_branch']:.2f} substitutions/site); the lowest node that joins them to a satellite clade
  (UFBoot {a5['ufboot']:g}, {a5['n_tips']} copies: the 5S genes + MacSat06 copies incl. {', '.join(a5['clades_below'])}) lies inside the MacSat06 part of the tree (MacSat06 = the *C. macrocephalus* 5S rDNA array, 72 % of its units carry an annotated 5S gene in v4.7);
  the nodes on the way there have UFBoot {'/'.join(f'{x:g}' for x in a5['ufboot_chain'])}, so the 5S-satellite relationship is supported at the last node only; the direction of the derivation is not resolved because the trees are unrooted.
* Copies of the centromeric tree that are not in a v4.7 run of their own family take the centromeric flag of the run covering the copy, where one exists.

Sources: treefiles above; `{ST.RUNS_VERSION}`: `runs_of_record/{{Gar,Mac}}/runs_of_record.tsv` (centromeric_L1), `satellite_family_nomenclature.tsv` (element_class, seq_cluster);
`fig1_panel_a_numbers.tsv` (family / clade / backbone / 5S rows), `fig1_panel_a_sweep.tsv`; code `fig1_panel_a.py` (`fig1_panel_a_combined.py` = same content on one common scale).
"""
    (HERE / f"{stem}_caption.md").write_text(txt, encoding="utf-8")


if __name__ == "__main__":
    ST.house_style(); W, H = 92, 56
    fig = ST.new_figure(W, H); num = draw(fig, 0, 0, W, H, "a")
    ST.save(fig, str(HERE / "_standalone" / "fig1_panel_a")); write_tables(num); write_caption(num)
    for k, v in num["models"].items(): print(k, "cut", v["cut"], "plateau", v["plateau"], "clades", len(v["clades"]), "assigned", sum(c["n"] for c in v["clades"]), "/", v["n_tips"])
    print(num["label_positions_mm"])
