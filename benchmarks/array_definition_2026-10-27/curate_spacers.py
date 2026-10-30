"""Manual-curation sheet for the spacers (interruptions) of one array (2026-10-27, v2). Sequence only, no annotation calls.
Spacers = gaps between consecutive catalogue runs (>= 200 bp), flipped to the family strand when the array is '-' in centromeres/array_strand.tsv
(k-mer-based strand). Representative = median-length spacer, drawn with F = 300 bp of the neighbouring satellite on each side (shaded) so the
junctions show.
Dot plots: 10-mer matches kept only when they continue a diagonal for >= 3 consecutive positions (sensitive to ~80 % homology, no random dots);
microsatellites (1-6 bp motif >= 16 bp) masked as pale-blue bands (microsatellite class colour); slate = same strand, magenta = reverse complement.
  A representative (with flanks) vs itself       B vs the LG04 RTE-type spacer int15       C vs the RTE source consensus
  D monomer similarity along the representative: 50-bp windows (step 10) aligned locally to the array's main monomer x 2 (both strands);
    identity of the best alignment (degraded satellite at the junctions shows as a rising curve)
  E-F vs the representative spacers of two other arrays
  G all-vs-all of the array's spacers (shared 15-mers, Jaccard), array order, with lengths
usage: python -P curate_spacers.py <array_id> <png> [other_array_id,other_array_id]"""
import sys, re, csv, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "F:/data/satellite_FISH/runs_of_record_v4.6_2026-10-27/code"); import palette_v46 as PAL
CM_ID = matplotlib.colors.LinearSegmentedColormap.from_list("magma_r_soft", plt.get_cmap("magma_r")(__import__("numpy").linspace(0, 0.82, 256)))   # identity map that stops before black
from Bio import Align
AID, PNG = sys.argv[1:3]; OTH = sys.argv[3].split(",") if len(sys.argv) > 3 else []
R0 = "F:/data/fClaHyb_func_annotation_v1.0"; DG = "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22"; FL = 300; BS = chr(92)
RCT = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan"); rc = lambda s: s.translate(RCT)[::-1]
STRAND = {r["array_id"]: r["strand_vs_family"] for r in csv.DictReader(open(f"{DG}/figure/diagnostics/centromeres/array_strand.tsv"), delimiter="\t")}
LOC = Align.PairwiseAligner(mode="local", match_score=1, mismatch_score=-1, open_gap_score=-2, extend_gap_score=-1)
class Fa:
    def __init__(self, p):
        self.fh = open(p, "rb"); self.idx = {l.split()[0]: tuple(map(int, l.split()[1:5])) for l in open(p + ".fai")}
    def get(self, c, a, b):
        L, off, lb, lby = self.idx[c]; a = max(0, a); b = min(L, b); self.fh.seek(off + (a // lb) * lby + a % lb)
        return self.fh.read((b - a) + ((b - a) // lb + 2) * (lby - lb)).replace(b"\n", b"").replace(b"\r", b"")[:b - a].upper()
def spacers(aid):
    sp = aid.split("_")[0]; C = pd.read_csv(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_satelliteDNA/fClaHyb_{sp}.satDNA_v46.catalogue.tsv", sep="\t")
    x = C[C.array_id == aid].sort_values("start0"); fa = Fa(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.fa"); ch = x.chrom.iloc[0]; rev = STRAND.get(aid) == "-"
    co = [(int(a), int(b)) for a, b in zip(x.end.values[:-1], x.start0.values[1:]) if b - a >= 200]
    if rev: co = co[::-1]
    seqs = [fa.get(ch, a, b) for a, b in co]; seqs = [rc(s) for s in seqs] if rev else seqs
    m = sorted(range(len(co)), key=lambda i: len(seqs[i]))[len(co) // 2] if co else None
    repf = fa.get(ch, co[m][0] - FL, co[m][1] + FL) if co else b""
    PO = f"F:/data/satellite_FISH/pipeline_v4_5_fixed/out/fClaHyb_{sp}_satellites/03_array_consensus.tsv"
    cons = max(((int(r["n_copies_in_class"] or 0), r["consensus_sequence"].upper().encode()) for r in csv.DictReader(open(PO), delimiter="\t") if r["array_id"] == aid), default=(0, b""))[1]
    if rev: repf, cons = rc(repf), rc(cons)
    return seqs, (seqs[m] if co else b""), repf, cons, ch
def microsats(s):
    t = s.decode(); out = []
    for k in range(1, 7):
        for m in re.finditer(("(([ACGT]{%d})" + BS + "2{%d,})") % (k, max(2, 12 // k)), t):
            if m.end() - m.start() >= 16: out.append((m.start(), m.end()))
    return out
def masked(s):
    m = np.zeros(len(s), dtype=bool)
    for a, b in microsats(s): m[a:b] = True
    return m
def dots(ax, a, b, title, k=10, run=3, shade_a=None):
    ma, mb = masked(a), masked(b)
    for x0, x1 in microsats(a): ax.axvspan(x0 / 1e3, x1 / 1e3, color=PAL.CLASS["microsatellite"], lw=0)
    for y0, y1 in microsats(b): ax.axhspan(y0 / 1e3, y1 / 1e3, color=PAL.CLASS["microsatellite"], lw=0)
    if shade_a: [ax.axvspan(x0 / 1e3, x1 / 1e3, color=PAL.shade(PAL.CLASS["satellite"], 0.85), lw=0, zorder=0) for x0, x1 in shade_a]
    idx, idr = {}, {}
    for j in range(len(b) - k + 1):
        if not mb[j:j + k].any(): idx.setdefault(b[j:j + k], []).append(j)
    rb = rc(b)
    for j in range(len(rb) - k + 1):
        jj = len(b) - k - j
        if not mb[jj:jj + k].any(): idr.setdefault(rb[j:j + k], []).append(jj)
    P = set(); Q = set()
    for i in range(len(a) - k + 1):
        w = a[i:i + k]
        if b"N" in w or ma[i:i + k].any(): continue
        for j in idx.get(w, ()): P.add((i, j))
        for j in idr.get(w, ()): Q.add((i, j))
    keepP = [(i, j) for i, j in P if all((i + d, j + d) in P for d in range(1, run)) or all((i - d, j - d) in P for d in range(1, run))]
    keepQ = [(i, j) for i, j in Q if all((i + d, j - d) in Q for d in range(1, run)) or all((i - d, j + d) in Q for d in range(1, run))]
    if keepP: x, y = zip(*keepP); ax.scatter(np.array(x) / 1e3, np.array(y) / 1e3, s=0.4, color="#3d4f66")
    if keepQ: x, y = zip(*keepQ); ax.scatter(np.array(x) / 1e3, np.array(y) / 1e3, s=0.4, color="#b5468a")
    ax.set_xlim(0, len(a) / 1e3); ax.set_ylim(0, len(b) / 1e3); ax.set_title(title, fontsize=9.5, loc="left"); ax.tick_params(labelsize=7)
    ax.set_xlabel("representative spacer with flanks (kb)", fontsize=8); ax.set_ylabel("kb", fontsize=8)
SP, REP, REPF, CONS, CH = spacers(AID)
INT15 = "".join(l.strip() for l in open(f"{DG}/split/spacer_content/int15.fa") if not l.startswith(">")).upper().encode()
RTE = "".join(l.strip() for l in open(f"{DG}/seq/Gar_full_consensus.fa") if not l.startswith(">")).upper().encode()
shade = [(0, FL), (len(REPF) - FL, len(REPF))]
fig = plt.figure(figsize=(18, 17)); fig.patch.set_facecolor("white")
pos = lambda i: [0.05 + (i % 3) * 0.32, 0.66 - (i // 3) * 0.30, 0.26, 0.24]
dots(fig.add_axes(pos(0)), REPF, REPF, f"A  representative vs itself ({len(REP):,} bp + {FL} bp satellite each side, shaded)", shade_a=shade)
dots(fig.add_axes(pos(1)), REPF, INT15, "B  vs LG04 RTE-type spacer int15", shade_a=shade)
dots(fig.add_axes(pos(2)), REPF, RTE, "C  vs RTE source consensus", shade_a=shade)
ax = fig.add_axes(pos(3)); ref = CONS * 2; xs, ys = [], []
for w0 in range(0, len(REPF) - 50, 10):
    w = REPF[w0:w0 + 50]; best = max(LOC.score(w, ref), LOC.score(rc(w), ref)) if b"N" not in w else 0
    xs.append((w0 + 25) / 1e3); ys.append(max(0.0, (best + 50) / 100))                 # score (m - x) over 50 bp -> identity ~ (score + 50) / 100
ax.axvspan(0, FL / 1e3, color=PAL.shade(PAL.CLASS["satellite"], 0.85), lw=0); ax.axvspan((len(REPF) - FL) / 1e3, len(REPF) / 1e3, color=PAL.shade(PAL.CLASS["satellite"], 0.85), lw=0)
ax.plot(xs, ys, color="#3d4f66", lw=1); ax.axhline(0.75, color=PAL.INK2, lw=0.8, ls=":"); ax.set_ylim(0.4, 1.02); ax.set_xlim(0, len(REPF) / 1e3)
ax.set_title("D  monomer similarity along the representative (50-bp windows; dotted = 0.75)", fontsize=9.5, loc="left"); ax.tick_params(labelsize=7)
ax.set_xlabel("representative spacer with flanks (kb)", fontsize=8); ax.set_ylabel("identity to the monomer", fontsize=8)
for n_, o in enumerate(OTH[:2]):
    _, r2, rf2, _, c2 = spacers(o); dots(fig.add_axes(pos(4 + n_)), REPF, rf2, f"{'EF'[n_]}  vs representative spacer of {o} ({c2.split('_')[-1]}, with flanks)", shade_a=shade)
k = 15; S = [set(s[j:j + k] for j in range(len(s) - k + 1) if b"N" not in s[j:j + k]) for s in SP]
J = np.array([[len(a & b) / max(1, len(a | b)) for b in S] for a in S]) if len(S) > 1 else np.ones((1, 1))
ax = fig.add_axes([0.05, 0.04, 0.3, 0.27]); im = ax.imshow(J, cmap=CM_ID, vmin=0, vmax=1); fig.colorbar(im, ax=ax, fraction=0.046, label="shared 15-mers (Jaccard)")
ax.set_title(f"G  all-vs-all of the {len(SP)} spacers (>= 200 bp), array order", fontsize=9.5, loc="left", pad=8)
ax2 = fig.add_axes([0.42, 0.04, 0.18, 0.27]); ax2.barh(range(len(SP)), [len(s) / 1e3 for s in SP], color="#888"); ax2.invert_yaxis(); ax2.set_xlabel("spacer length (kb)", fontsize=8); ax2.set_yticks([])
ax2.set_title("spacer lengths", fontsize=9.5, loc="left", pad=8)
fig.text(0.02, 0.975, f"{AID} ({CH.split('_')[-1]}, strand {STRAND.get(AID, '?')} vs family): spacer curation sheet, median spacer {len(REP):,} bp", fontsize=13, fontweight="bold")
fig.text(0.02, 0.955, "dot plots: 10-mers kept when >= 3 consecutive along a diagonal; slate = same strand, magenta = reverse complement; pale-blue bands = microsatellites (masked); "
         "pale-teal shading = 300 bp of satellite flanking the representative", fontsize=9.5)
fig.savefig(PNG, dpi=450); fig.savefig(PNG.replace(".png", ".pdf")); print("wrote", PNG, "| spacers:", len(SP), "| median", len(REP), "bp")