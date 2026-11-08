"""Figure 2, panel f (house style module): one interrupted MacSat03 array (C. macrocephalus LG04, array Mac_SAT00318) to scale (100 x 46 mm).
   draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None) -> dict of the numbers shown.
   Standalone:  python fig2_panel_f.py  -> panels/_standalone/fig2_panel_f.{png,pdf}, fig2_panel_f_numbers.tsv
DATA
  runs        satchain_v4.7_2026-11-01/runs_of_record/Mac/runs_of_record.tsv, array_id Mac_SAT00318 (28 v4.7 runs of record, same id and runs as v4.6; NOT array extents). Interruption k = gap between run k and run k+1 (27).
              Unit k = run k + interruption k (= run start -> next run start; split/units_satellite_side.tsv gives the same coordinates, asserted).
  copies      pipeline_v4_5_fixed/out/fClaHyb_Mac_satellites/03_monomer_copies.tsv (SatChain step 03): every monomer copy fully inside a run = one cell, coloured by its class (C01..C09);
              'partial copy' = truncated, or identity < 80 % to its class (sweep 70 / 80 / 90 %).
  RTE         LTR_MacSat03_2026-10-22/split/rte_genome_hits_both.tsv (dc-megablast of the RTE consensus, HSP >= 150 bp), hits inside interruptions; arrow = RTE 5'->3' (consensus q = 1 is the 3' end,
              so a '-' hit points along the chromosome, a '+' hit against it); sweep of the minimum HSP length 150 / 500 / 1000 / 2000 bp.
  ins         split/spacer_content/nonref_inserts.tsv, the 245-bp insertion, placed with the interruption offsets of split/spacer_decomp/LG04_interruptions.fa.
  spacer types  A = interruption carries the 245-bp ins; C = no ins, 12 variants of its own (exact-support split); B = remainder (no ins, no variant of its own). Source tables: split/spacer_types_msa/
              (type_support*.tsv; MAFFT --auto and --localpair alignments). Re-tested here (a) A: share of the ins 15-mers found (both strands) in each interruption, threshold 0.1-0.9;
              (b) B vs C: unsupervised average-linkage cut (2 clusters) of the 15 no-ins interruptions on alignment substitution columns, minor allele count >= 2 / 3 / 4, both alignments.
  centromere  satchain_v4.7_2026-11-01/runs_of_record/Mac/centromere_check.tsv (LG04 interval cen_start0-cen_end, pale grey band; same interval as v4.6)."""
import sys, csv, pathlib, functools, collections
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "code"))
import satellite_style as ST
import fig1_efg_common as C
from matplotlib.patches import Rectangle, Polygon

AID = "Mac_SAT00318"
W = pathlib.Path(r"F:\data\satellite_FISH\LTR_MacSat03_2026-10-22"); RR = ST.RUNS
COPIES = r"F:\data\satellite_FISH\pipeline_v4_5_fixed\out\fClaHyb_Mac_satellites\03_monomer_copies.tsv"
MSA = W / "split" / "spacer_types_msa"
CHROM = "C_mac_LG04"
# monomer-class colours = class_colours.py (same classes, same colours as the stack / linear figures of this array), muted a touch toward white for the paper style
CLASS_ORDER = ["C01", "C02", "C03", "C04", "C06", "C07", "C08", "C09"]
CLASS_COL = {"C01": "#1f4e9c", "C02": "#1b9e8a", "C04": "#7b4fa3", "C03": "#c0392b", "C07": "#7f9f2f", "C09": "#4a90d9", "C08": "#b5468a", "C06": "#2e8b3d"}
PARTIAL = "#b8b6ae"
TYPE_COL = {"A": "#D08C6E", "B": "#B8A9CC", "C": "#71A29C"}                    # spacer types (house terracotta / lavender / sage); a separate lane, never next to the class cells
INS_COL = "#83251D"
BAND = "#e6e2d6"                                                               # pale centromere context
FS = 6.5


def _mute(h, f=0.12):
    from matplotlib.colors import to_rgb
    c = np.array(to_rgb(h)); return tuple(c + (1 - c) * f)


@functools.lru_cache(maxsize=1)
def data():
    R = pd.read_csv(RR / "Mac" / "runs_of_record.tsv", sep="\t"); R = R[R.array_id == AID].sort_values("start0").reset_index(drop=True)
    assert len(R) == 28 and (R.family == "MacSat03").all()
    gaps = list(zip(R.end.values[:-1], R.start0.values[1:]))
    U = pd.read_csv(W / "split" / "units_satellite_side.tsv", sep="\t"); U = U[U.chrom == "LG04"].reset_index(drop=True)
    assert len(U) == 27 and (U.unit_lo.values == R.start0.values[:-1]).all() and (U.unit_hi.values == R.start0.values[1:]).all()
    cen = pd.read_csv(RR / "Mac" / "centromere_check.tsv", sep="\t").set_index("chrom").loc[CHROM]; cen = cen.rename({"cen_start0": "cen_start"})   # v4.7 column name
    # monomer copies inside the runs
    iv = list(zip(R.start0, R.end)); rows = []
    with open(COPIES, newline="") as fh:
        hdr = next(fh).rstrip("\n").split("\t"); ix = {h: i for i, h in enumerate(hdr)}
        for ln in fh:
            if AID not in ln: continue
            p = ln.rstrip("\n").split("\t")
            if p[ix["array_id"]] != AID: continue
            s, e = int(p[ix["genomic_start_0based"]]), int(p[ix["genomic_end"]])
            if any(s >= a and e <= b for a, b in iv):
                rows.append((s, e, p[ix["monomer_class"]], p[ix["complete_or_truncated"]] == "complete", float(p[ix["best_alignment_identity"]] or 0)))
    cop = pd.DataFrame(rows, columns=["s", "e", "cls", "complete", "ident"]).sort_values("s").reset_index(drop=True)
    # RTE hits
    H = pd.read_csv(W / "split" / "rte_genome_hits_both.tsv", sep="\t"); H = H[H.chrom == CHROM].copy()
    H["lo"] = H[["sstart", "send"]].min(axis=1) - 1; H["hi"] = H[["sstart", "send"]].max(axis=1); H["plus"] = H.send > H.sstart
    ia, ib = R.start0.min(), R.end.max(); H = H[(H.lo >= ia) & (H.hi <= ib)]
    k = []
    for r in H.itertuples():                                              # which interruption
        g = [i for i, (a, b) in enumerate(gaps) if r.lo >= a - 1 and r.hi <= b + 1]
        k.append(g[0] + 1 if g else 0)
    H["int"] = k
    # interruption sequences + ins positions
    seqs, off = {}, {}
    for blk in open(W / "split" / "spacer_decomp" / "LG04_interruptions.fa").read().split(">")[1:]:
        h, s = blk.split("\n", 1); nm = h.split("|")[0]; seqs[nm] = s.replace("\n", "").upper(); off[nm] = int(h.split("|")[1].split(":")[1].split("-")[0])
    NI = pd.read_csv(W / "split" / "spacer_content" / "nonref_inserts.tsv", sep="\t"); NI = NI[NI.bp == 245]
    INS = pd.DataFrame([(off[r.interruption] + r.q_start - 1, off[r.interruption] + r.q_end, int(r.interruption[3:])) for r in NI.itertuples()], columns=["lo", "hi", "int"])
    return R, gaps, U, cen, cop, H, seqs, INS


RC = str.maketrans("ACGT", "TGCA")
def _kmers(s, k): return {s[i:i + k] for i in range(len(s) - k + 1)}


@functools.lru_cache(maxsize=1)
def spacer_types():
    """A / B / C for interruptions 2..27 (1 = truncated end piece, not typed) + the re-tests."""
    R, gaps, U, cen, cop, H, seqs, INS = data()
    n = len(gaps); info = {}
    ref = seqs["int03"][3991:4236]; assert len(ref) == 245
    K = 15; kr = _kmers(ref, K) | _kmers(ref.translate(RC)[::-1], K)
    share = {}
    for i in range(2, n + 1):                                             # share of the ins 15-mers present in the interruption, both strands
        s = seqs[f"int{i:02d}"]; ks = _kmers(s, K) | _kmers(s.translate(RC)[::-1], K); share[i] = len(kr & ks) / len(kr)
    A_set = set(int(x) for x in INS["int"]); A_sweep = {}
    for th in (0.1, 0.3, 0.5, 0.7, 0.9): A_sweep[th] = sorted(i for i, v in share.items() if v >= th)
    Bn = [2, 10, 11, 12, 14, 16, 17, 18, 22, 23]; Cn = [8, 9, 13, 15, 21]                   # no-ins interruptions (the two groups to be separated)
    # (b) unsupervised 2-cluster cut of the no-ins interruptions on alignment substitution columns
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    names = [f"int{i:02d}" for i in Bn + Cn]; cut = {}
    for tag in ("auto", "linsi"):
        A = {}; nm = None
        for l in open(MSA / f"aln_{tag}.fa"):
            l = l.strip()
            if l.startswith(">"): nm = l[1:]; A[nm] = []
            else: A[nm].append(l.upper())
        M = np.array([list("".join(A[x])) for x in names]); data_ = np.zeros(M.shape, bool)
        for i in range(len(names)):
            nz = np.where(M[i] != "-")[0]; data_[i, nz.min():nz.max() + 1] = True
        cols = []
        for j in range(M.shape[1]):
            h = [i for i in range(len(names)) if data_[i, j] and M[i, j] != "-"]
            if len(h) < len(names) - 3: continue
            cn = collections.Counter(M[i, j] for i in h)
            if len(cn) == 2: cols.append((j, min(cn.values())))
        for mm in (2, 3, 4):
            cc = [j for j, m in cols if m >= mm]
            D = np.zeros((len(names), len(names)))
            for a in range(len(names)):
                for b in range(len(names)):
                    u = [j for j in cc if M[a, j] != "-" and M[b, j] != "-"]
                    D[a, b] = np.mean([M[a, j] != M[b, j] for j in u]) if u else 0.5
            lab = fcluster(linkage(squareform(D, checks=False), "average"), 2, "maxclust")
            g = {k: frozenset(int(names[i][3:]) for i in range(len(names)) if lab[i] == k) for k in set(lab)}
            cut[(tag, mm)] = (len(cc), sorted(len(v) for v in g.values()), frozenset(Cn) in g.values())
    types = {}
    for i in range(2, n + 1): types[i] = "A" if i in A_set else ("C" if i in Cn else "B")
    assert sorted(i for i in types if i in Bn and types[i] != "B") == [] and sorted(Cn) == sorted(i for i, t in types.items() if t == "C")
    return types, A_set, share, A_sweep, cut


@functools.lru_cache(maxsize=1)
def sweeps():
    R, gaps, U, cen, cop, H, seqs, INS = data(); out = {}
    base = np.where(cop.complete & (cop.ident >= 80), cop.cls, "U")
    for th in (70, 80, 90):
        t = np.where(cop.complete & (cop.ident >= th), cop.cls, "U"); out[f"ident{th}"] = (int((t == "U").sum()), int((t != base).sum()))
    Hf = pd.read_csv(W / "split" / "rte_genome_hits_both.tsv", sep="\t"); Hf = Hf[Hf.chrom == CHROM].copy()
    Hf["lo"] = Hf[["sstart", "send"]].min(axis=1) - 1; Hf["hi"] = Hf[["sstart", "send"]].max(axis=1)
    Hf = Hf[(Hf.lo >= R.start0.min()) & (Hf.hi <= R.end.max())]
    for ml in (150, 500, 1000, 2000):
        g = Hf[Hf.len >= ml]; out[f"rte{ml}"] = (len(g), len(set(H_int(g, gaps))))
    return out


def H_int(g, gaps):
    r = []
    for x in g.itertuples():
        for i, (a, b) in enumerate(gaps):
            if x.lo >= a - 1 and x.hi <= b + 1: r.append(i + 1)
    return r


def _rect(ax, x0, x1, y0, y1, fc, z=2, ec="none", lw=0):
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z))


def _rte(ax, x0, x1, yc, h, forward, head):
    C.arrow(ax, x0, x1, yc, h, ST.RTE, C.GOLD_EDGE, forward=forward, head=head, lw=0.25, z=4)


def draw(fig, left_mm, top_mm, w_mm=100.0, h_mm=46.0, letter=None):
    R, gaps, U, cen, cop, H, seqs, INS = data(); types, A_set, share, A_sweep, cut = spacer_types(); sw = sweeps()
    if letter: ST.letter(fig, left_mm, top_mm, letter)
    ax = C.mm_axes(fig, left_mm, top_mm, w_mm, h_mm)
    T = lambda x, y, s, **k: C.txt(ax, x, y, s, **k)
    # ---- title: bold + italic species
    t1 = ax.text(4.2, 0.2, "One interrupted MacSat03 array, ", fontsize=7, fontweight="bold", ha="left", va="top")
    fig.canvas.draw(); bb = t1.get_window_extent(); inv = ax.transData.inverted(); x_end = inv.transform((bb.x1, bb.y0))[0]
    ax.text(x_end, 0.2, "C. macrocephalus", fontsize=7, fontweight="bold", fontstyle="italic", ha="left", va="top")
    t2 = ax.text(x_end, 0.2, "", fontsize=7)
    # ---- geometry: overview = whole centromere interval (pale band), array to scale
    c0, c1 = float(cen.cen_start), float(cen.cen_end); XL, XR = 2.0, w_mm - 2.0; sc = (XR - XL) / (c1 - c0); X = lambda p: XL + (p - c0) * sc     # mm per bp
    a0, a1 = int(R.start0.min()), int(R.end.max())
    # zoom window: units 7..10
    zu0, zu1 = 7, 10; z0 = int(R.start0[zu0 - 1]) - 300; z1 = int(R.start0[zu1]); zsc = (XR - XL) / (z1 - z0); Z = lambda p: XL + (p - z0) * zsc
    cls_col = {k: _mute(v) for k, v in CLASS_COL.items()}; cls_col["U"] = PARTIAL
    base = np.where(cop.complete & (cop.ident >= 80), cop.cls, "U")
    # ---------------- overview strip ----------------
    yb = 4.9                                                                # zoom bracket (above the strip)
    ZX0, ZX1 = X(z0), X(z1); ax.plot([ZX0, ZX0, ZX1, ZX1], [yb + 1.0, yb, yb, yb + 1.0], color=ST.NAVY, lw=0.5, solid_capstyle="butt", zorder=6)
    T((ZX0 + ZX1) / 2 + 0.0, yb - 0.1, "", size=FS)
    yo, ho = 7.2, 3.2                                                       # overview track: y centre / run height
    _rect(ax, X(c0), X(c1), yo - 2.1, yo + 2.1, BAND, z=1)
    T(XL + 0.4, yo + 3.45, "", size=FS)
    for r in R.itertuples():                                                # runs: cells
        pass
    for s, e, c in zip(cop.s, cop.e, base):
        ax.add_patch(Rectangle((X(s), yo - ho / 2), max(X(e) - X(s), 0.02), ho, facecolor=cls_col[c], edgecolor=cls_col[c], linewidth=0.12, zorder=3))
    for r in H.itertuples():                                                # RTE arrows, orientation = strand
        _rte(ax, X(r.lo), X(r.hi), yo, 1.9, forward=not r.plus, head=0.35)
    ys = yo + 2.35; hb = 1.1
    for i, (g0, g1) in enumerate(gaps, 1):                                  # spacer type bar + letter
        if i in types:
            _rect(ax, X(g0), X(g1), ys, ys + hb, TYPE_COL[types[i]], z=3)
            T((X(g0) + X(g1)) / 2, ys + hb + 1.15, types[i], size=FS, ha="center")
    # scale bar + extent
    yscale = ys + hb + 3.6
    C.hline(ax, X(c0), X(c0) + 50000 * sc, yscale); C.vline(ax, X(c0), yscale - 0.5, yscale + 0.5); C.vline(ax, X(c0) + 50000 * sc, yscale - 0.5, yscale + 0.5)
    T(X(c0) + 50000 * sc + 1.0, yscale, "50 kb", size=FS)
    T(XR, yscale, f"LG04 {c0 / 1e6:.2f}\u2013{c1 / 1e6:.2f} Mb", size=FS, ha="right")
    # ---------------- zoom ----------------
    yz = yscale + 9.6; hz = 5.0                                             # zoom track centre / run height
    _rect(ax, Z(z0), Z(z1), yz - 1.0, yz + 1.0, BAND, z=1)
    for s, e, c in zip(cop.s, cop.e, base):
        if e > z0 and s < z1: ax.add_patch(Rectangle((Z(s), yz - hz / 2), Z(e) - Z(s), hz, facecolor=cls_col[c], edgecolor="white", linewidth=0.12, zorder=3))
    for r in H.itertuples():
        if r.hi > z0 and r.lo < z1: _rte(ax, Z(r.lo), Z(r.hi), yz, 2.9, forward=not r.plus, head=1.1)
    for r in INS.itertuples():
        if r.hi > z0 and r.lo < z1: _rect(ax, Z(r.lo) - 0.12, Z(r.hi) + 0.12, yz - 1.7, yz + 1.7, INS_COL, z=5)
    yt = yz + hz / 2 + 0.5; htz = 2.5
    for i, (g0, g1) in enumerate(gaps, 1):
        if i in types and zu0 <= i <= zu1:
            ga, gb = max(g0, z0), min(g1, z1); _rect(ax, Z(ga), Z(gb), yt, yt + htz, TYPE_COL[types[i]], z=3)
            if (gb - ga) * zsc > 3.0: T((Z(ga) + Z(gb)) / 2, yt + htz / 2 + 0.05, types[i], size=FS, ha="center", weight="bold")
    yu = yt + htz + 1.5                                                     # unit brackets
    for k in range(zu0, zu1 + 1):
        a, b = int(R.start0[k - 1]), int(R.start0[k]); xa, xb = Z(a) + 0.2, Z(b) - 0.2
        ax.plot([xa, xa, xb, xb], [yu - 0.6, yu, yu, yu - 0.6], color=ST.NAVY, lw=0.5, solid_capstyle="butt", zorder=6)
        T((xa + xb) / 2, yu + 1.3, "unit", size=FS, ha="center")
    ysz = yz - hz / 2 - 2.2                                                  # zoom scale bar: 5 kb, top right of the zoom
    C.hline(ax, XR - 5000 * zsc, XR, ysz); C.vline(ax, XR - 5000 * zsc, ysz - 0.5, ysz + 0.5); C.vline(ax, XR, ysz - 0.5, ysz + 0.5); T(XR - 5000 * zsc - 1.0, ysz, "5 kb", size=FS, ha="right")
    # ---------------- legends (text widths measured so that nothing overlaps) ----------------
    rend = fig.canvas.get_renderer(); inv = ax.transData.inverted()
    def put(x, y, s, **k):
        t = T(x, y, s, size=FS, **k); bb = t.get_window_extent(rend); return inv.transform((bb.x1, bb.y0))[0]
    GAP = 2.6
    ly = yu + 5.0; x = 2.0; x = put(x, ly, "monomer class") + 1.6
    for c in CLASS_ORDER + ["U"]:
        _rect(ax, x, x + 1.8, ly - 0.9, ly + 0.9, cls_col[c], z=3); x = put(x + 2.5, ly, "partial copy" if c == "U" else c) + GAP * 0.55
    ly2 = ly + 3.9; x = 2.0; nA = sum(t == "A" for t in types.values()); nB = sum(t == "B" for t in types.values()); nCc = sum(t == "C" for t in types.values())
    x = put(x, ly2, "spacer") + 1.6
    for t, lab in (("A", f"A  with ins ({nA})"), ("B", f"B  remainder ({nB})"), ("C", f"C  own variants ({nCc})")):
        _rect(ax, x, x + 1.8, ly2 - 0.9, ly2 + 0.9, TYPE_COL[t], z=3); x = put(x + 2.5, ly2, lab) + GAP
    ly3 = ly2 + 3.9; x = 2.0
    _rte(ax, x, x + 3.6, ly3, 1.9, True, 1.1); x = put(x + 4.4, ly3, "RTE, arrow 5′→3′") + GAP
    _rect(ax, x, x + 0.5, ly3 - 1.0, ly3 + 1.0, INS_COL, z=3); x = put(x + 1.5, ly3, "245-bp ins") + GAP
    _rect(ax, x, x + 3.0, ly3 - 0.9, ly3 + 0.9, BAND, z=2); x = put(x + 3.8, ly3, "candidate centromere") + GAP
    # ---------------- numbers ----------------
    kb = lambda a, b: round((b - a) / 1e3, 1)
    sat_bp = int((R.end - R.start0).sum()); gap_bp = int(sum(b - a for a, b in gaps))
    ncls = pd.Series(base).value_counts().to_dict()
    out = {"array": AID, "chrom": "C_mac_LG04", "array_start_bp": a0, "array_end_bp": a1, "array_span_kb": kb(a0, a1), "centromere_start": int(c0), "centromere_end": int(c1), "centromere_kb": kb(c0, c1),
           "runs": len(R), "interruptions": len(gaps), "satellite_in_runs_kb": round(sat_bp / 1e3, 1), "interruptions_kb": round(gap_bp / 1e3, 1), "copies_in_runs": len(cop),
           "RTE_arrows_total": len(H), "RTE_arrows_forward_(minus_hit)": int((~H.plus).sum()), "RTE_arrows_reverse_(plus_hit)": int(H.plus.sum()),
           "interruptions_with_RTE": int(H["int"].nunique()), "ins_245bp_copies": len(INS), "ins_interruptions": len(set(INS["int"])),
           "types_A": nA, "types_B": nB, "types_C": nCc, "untyped_interruption_1_truncated_end_piece_bp": int(gaps[0][1] - gaps[0][0]),
           "type_per_interruption": ",".join(f"{i}{t}" for i, t in sorted(types.items())),
           "unit_def": "run + the interruption that follows it (run start -> next run start)", "zoom_units": f"{zu0}-{zu1}", "zoom_span_kb": kb(z0, z1), "overview_mm_per_kb": round(sc * 1e3, 4), "zoom_mm_per_kb": round(zsc * 1e3, 4),
           "zoom_types": "".join(types[i] for i in range(zu0, zu1 + 1)),
           "unit_bp_median": int(np.median(U.period)), "unit_bp_min": int(U.period.min()), "unit_bp_max": int(U.period.max()),
           "evidence_level_spacer_types": "L2"}
    for c in CLASS_ORDER + ["U"]: out[f"copies_{c}"] = int(ncls.get(c, 0))
    for th in (70, 80, 90): out[f"sweep_ident{th}_partial_copies"], out[f"sweep_ident{th}_copies_changed_vs80"] = sw[f"ident{th}"]
    for ml in (150, 500, 1000, 2000): out[f"sweep_RTE_minHSP{ml}_arrows"], out[f"sweep_RTE_minHSP{ml}_interruptions"] = sw[f"rte{ml}"]
    out["ins_15mer_share_min_A_set"] = round(min(share[i] for i in A_set if i in share), 3); out["ins_15mer_share_max_non_A"] = round(max(v for i, v in share.items() if i not in A_set), 3)
    for th, v in A_sweep.items(): out[f"sweep_A_kmer_share_ge{th}_n"] = len(v); out[f"sweep_A_kmer_share_ge{th}_equals_ins_set"] = bool(set(v) == (A_set - {1}))
    for (tag, mm), (nc, sizes, ok) in cut.items(): out[f"sweep_BvsC_cut_{tag}_minor{mm}"] = f"{nc} columns, clusters {sizes[0]}+{sizes[1]}, = C set: {ok}"
    out["sweep_BvsC_all_six_recover_C"] = all(v[2] for v in cut.values())
    return out


if __name__ == "__main__":
    Wm, Hm = 100.0, 46.0
    res = C.standalone("fig2_panel_f", draw, Wm, Hm, "f")
    C.write_numbers(HERE / "fig2_panel_f_numbers.tsv", [(k, v, "") for k, v in res.items()], header=("item", "value", "note"))
    print(res)
