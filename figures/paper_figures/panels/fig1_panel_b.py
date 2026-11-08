#!/usr/bin/env python3
"""Figure 1 panel b (88 x 56 mm): nucleotide view of monomer alignments, house style.  Centromeric group (GarSat03, MacSat03: the RTE-carrying centromeric families, same stripes in
both species) and 5S-related group (GarSat04 = 5S-derived satellite; MacSat06 = the C. macrocephalus 5S rDNA array, a gene array: a 5S-related segment in every monomer, bracketed).
2 x 2 sub-panels, rows = consecutive monomers of the example run of that family in the SatChain v4.7 runs of record (rule changed 2026-11-01: v4.7 has no rows for Gar LG23 / LG27, the
example chromosomes of the 2026-10-27 version; the example run is now, among the 5 largest S1 runs of the family, the one with the most homogeneous 30-monomer stack), cut at the repeat period, each row
phase-corrected to the row above (shift +-SHIFT bp); Mac consensus rotated onto the Gar one; 5S-related arrays oriented so the 5S gene reads forward; the segment = local alignment of a 5S rRNA gene to the row consensus.
Pipeline copied from thesis_defense_slides_tpl/code_R34/plot_slide_G07c_nucleotide_view.py (same parameters NROW 30, SHIFT 15, SKIP 3000, same period source = v4.5.3 unit_bp of the array).
usage: draw(fig, left_mm, top_mm, w_mm, h_mm, letter) -> dict ; python fig1_panel_b.py renders it alone and writes numbers + sweep tables."""
import sys, pathlib, collections, csv, random, json
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "code"))
import satellite_style as ST
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle
from Bio import SeqIO, Align
from Bio.Seq import Seq

A_ = pathlib.Path("F:/data/fClaHyb_func_annotation_v1.0"); R = ST.RUNS
R_OLD = pathlib.Path("F:/data/satellite_FISH/runs_of_record_v4.5.3_2026-10-18")
FIVE_S_FA = "F:/data/satellite_FISH/share_LITE_v4.5_2026-10-15/phylogeny/lineageB_5S_derived/lineageB_5S_derived_copies.fa"
NROW0, SHIFT0, SKIP0 = 30, 15, 3000
PANELS = [("cen", "Gar", "GarSat03"), ("cen", "Mac", "MacSat03"), ("5S", "Gar", "GarSat04"), ("5S", "Mac", "MacSat06")]     # (group, haplotype, family); chromosome = largest v4.7 S1 run
SPNAME = {"Gar": "C. gariepinus", "Mac": "C. macrocephalus"}
NUCCOL = [ST.NUC[c] for c in "ACGT"] + ["#D9DCDE"]
CMAP = ListedColormap(NUCCOL); IDX = {c: i for i, c in enumerate("ACGT")}
AL = Align.PairwiseAligner(mode="local", match_score=1, mismatch_score=-1, open_gap_score=-2, extend_gap_score=-1)
FIVE_S = next(str(r.seq) for r in SeqIO.parse(FIVE_S_FA, "fasta") if "|5S_rRNA|" in r.id)
_UNIT = None
# panel-g scoring (2026-10-28): medoid 115-bp 5S gene (fig1_panel_g_data/summary.json), local alignment match 2 / mismatch -3 / gap -5 / -2, both strands, circular monomer.  The 5S-related segment (box) and the
# orientation of the 5S-related stacks are decided with THIS scoring, so that panels b and g agree (the earlier default scoring 1/-1/-2/-1 gave a 60 / 80 bp box and the opposite strand for GarSat04)
GENE = json.load(open(str(HERE / "fig1_panel_g_data" / "summary.json")))["ref_seq"]
ALG = Align.PairwiseAligner(mode="local", match_score=2, mismatch_score=-3, open_gap_score=-5, extend_gap_score=-2)
rc = lambda x: str(Seq(x).reverse_complement())


def five_hit(cons):
    """(start, end, score_plus, score_minus): medoid gene vs the circular consensus (+ strand); end may exceed len(cons) = the segment wraps the seam"""
    q = cons + cons[:len(GENE)]; a = ALG.align(q, GENE)[0]
    qm = rc(cons); qm = qm + qm[:len(GENE)]
    return int(a.aligned[0][0][0]), int(a.aligned[0][-1][1]), float(a.score), float(ALG.score(qm, GENE))


def unit_old():
    """stack period of every (species, family, chrom) = unit_bp of the largest S1 run in the v4.5.3 run table (the v4.6 / v4.7 unit_bp of an array is not always the repeat period of that array)"""
    global _UNIT
    if _UNIT is None:
        _UNIT = {}
        for h in ("Gar", "Mac"):
            o = pd.read_csv(R_OLD / h / "runs_of_record.tsv", sep="\t"); o = o[o.tier == "S1"]
            for (f, c), g in o.groupby(["family", "chrom"]): _UNIT[(h, f, c)] = int(g.sort_values("run_bp", ascending=False).unit_bp.iloc[0])
    return _UNIT


def fetch(h, chrom, a, b):
    fa = A_ / f"fClaHyb_{h}/fClaHyb_{h}_genome/fClaHyb_{h}_LG.fa"
    for l in open(str(fa) + ".fai"):
        n, ln, off, lb, lby = l.split("\t")[:5]
        if n == chrom: ln, off, lb, lby = int(ln), int(off), int(lb), int(lby); break
    with open(fa, "rb") as fh:
        s0 = off + (a // lb) * lby + a % lb; s1 = off + (b // lb) * lby + b % lb
        fh.seek(s0); return fh.read(s1 - s0).replace(b"\n", b"").replace(b"\r", b"").decode().upper()


def stack(seq, P, nrow, shift):
    rows, pos = [], 0
    for _ in range(nrow):
        if not rows: rows.append(seq[pos:pos + P]); pos += P; continue
        prev = rows[-1]; best = min(range(-shift, shift + 1), key=lambda d: sum(x != y for x, y in zip(prev, seq[pos + d:pos + d + P])))
        pos += best; rows.append(seq[pos:pos + P]); pos += P
    return rows


def consensus(rows): return "".join(collections.Counter(c).most_common(1)[0][0] for c in zip(*rows))


def best_rot(ref, q):
    out = []
    for strand, s in (("+", q), ("-", str(Seq(q).reverse_complement()))):
        for k in range(len(s)): out.append((sum(x == y for x, y in zip(ref, s[k:] + s[:k])), k, strand))
    return max(out)


NCAND = 5
_PICK = {}


def pick_run(h, f, nrow=NROW0, shift=SHIFT0, skip=SKIP0):
    """the example run of a family = among the NCAND largest S1 runs of the family in the v4.7 runs of record, the one whose 30-monomer stack is the most homogeneous (highest mean
    row identity to the stack consensus): a well-homogenised, typical stretch.  All candidates are written to the numbers table (rows 'candidate')."""
    key = (h, f)
    if key in _PICK: return _PICK[key][0]
    r = ST.runs(h); r = r[(r.family == f) & (r.tier == "S1")].sort_values("run_bp", ascending=False).head(NCAND)
    cand = []
    for x in r.itertuples():
        c = str(x.chrom); P = unit_old().get((h, f, c), int(x.unit_bp)); a = int(x.start0) + min(skip, int(x.run_bp) // 4)
        rows = stack(fetch(h, c, a, a + P * (nrow + 12)), P, nrow, shift); cons = consensus(rows)
        rid = float(np.mean([sum(p == q for p, q in zip(rr, cons)) / len(cons) for rr in rows]))
        cand.append((rid, x, dict(family=f, species=h, chrom=c, run_bp=int(x.run_bp), run_start0=int(x.start0), period_bp=P, row_identity_pct=round(100 * rid, 1), centromeric_L1=x.centromeric_L1, cen_candidate_ror=x.cen_candidate_ror)))
    best = max(cand, key=lambda t: t[0]); _PICK[key] = (best[1], [t[2] for t in cand]); return best[1]


def candidates():
    return [dict(chosen=int(d["chrom"] == str(_PICK[k][0].chrom) and d["run_start0"] == int(_PICK[k][0].start0)), **d) for k in _PICK for d in _PICK[k][1]]


def build(nrow=NROW0, shift=SHIFT0, skip=SKIP0):
    data = {}
    for grp, h, f in PANELS:
        r = pick_run(h, f); c = str(r.chrom)
        P = unit_old().get((h, f, c), int(r.unit_bp)); a = int(r.start0) + min(skip, int(r.run_bp) // 4); seq = fetch(h, c, a, a + P * (nrow + 12))
        if grp == "5S":
            f_ = AL.score(FIVE_S, seq[:4 * P]); r_ = AL.score(FIVE_S, str(Seq(seq[:4 * P]).reverse_complement()))
            if r_ > f_: seq = str(Seq(seq).reverse_complement())
        data[f] = dict(h=h, c=c, a=a, P=P, run_bp=int(r.run_bp), run_start0=int(r.start0), v47_unit=int(r.unit_bp), period_from_v453=((h, f, c) in unit_old()), unit_class=str(r.unit_class),
                       element_class=str(r.element_class), seq_cluster=str(r.seq_cluster), seq=seq, rows=stack(seq, P, nrow, shift), grp=grp)
    for grp, gf, mf in (("cen", "GarSat03", "MacSat03"), ("5S", "GarSat04", "MacSat06")):
        G = data[gf]
        if grp == "5S":
            cg = consensus(G["rows"]); aln = AL.align(FIVE_S, cg + cg)[0]; st = aln.aligned[1][0][0] % len(cg)
            k = (st - 20) % len(cg); G["rows"] = [r_[k:] + r_[:k] for r_ in G["rows"]]
        cg = consensus(G["rows"]); M = data[mf]; cm = consensus(M["rows"])
        sc, k, strand = best_rot(cg, cm)
        rows = M["rows"] if strand == "+" else [str(Seq(r_).reverse_complement()) for r_ in M["rows"]]
        M["rows"] = [r_[k:] + r_[:k] for r_ in rows]; M["phase_id"] = sc / min(len(cg), len(cm))
    G, M = data["GarSat04"], data["MacSat06"]; cg = consensus(G["rows"]); _, _, sp_, sm_ = five_hit(cg)
    G["flipped"] = M["flipped"] = bool(sm_ > sp_)
    if G["flipped"]:        # panel-g scoring finds the gene on the other strand: reverse-complement both 5S-related stacks, rotate so the segment starts at column 20, re-phase Mac onto Gar
        G["rows"] = [rc(r_) for r_ in G["rows"]]; s0 = five_hit(consensus(G["rows"]))[0]; k = (s0 - 20) % G["P"]; G["rows"] = [r_[k:] + r_[:k] for r_ in G["rows"]]
        rows = [rc(r_) for r_ in M["rows"]]; sc, k, strand = best_rot(consensus(G["rows"]), consensus(rows))
        rows = rows if strand == "+" else [rc(r_) for r_ in rows]; M["rows"] = [r_[k:] + r_[:k] for r_ in rows]; M["phase_id"] = sc / min(G["P"], M["P"])
    for f, d in data.items():
        cons = consensus(d["rows"]); d["cons"] = cons
        d["row_id"] = float(np.mean([sum(x == y for x, y in zip(r_, cons)) / len(cons) for r_ in d["rows"]]))
        s0, s1, scp, scm = five_hit(cons); d["five_raw"] = (s0 % len(cons), s0 % len(cons) + (s1 - s0), int(scp)); d["five_strand_scores"] = (scp, scm)
    return data


def five_s_null(d, nshuf=200, seed=1):
    """for every displayed row: local alignment score of the 5S gene vs the row, vs the 99th percentile of the same score on mononucleotide-shuffled copies of the row (n shuffles).  Returns fraction of rows above the shuffled 99th percentile"""
    rnd = random.Random(seed); hit = 0
    for r_ in d["rows"]:
        s = AL.score(FIVE_S, r_); s_rc = AL.score(FIVE_S, str(Seq(r_).reverse_complement())); s = max(s, s_rc)
        sh = []
        for _ in range(nshuf):
            l = list(r_); rnd.shuffle(l); l = "".join(l); sh.append(max(AL.score(FIVE_S, l), AL.score(FIVE_S, str(Seq(l).reverse_complement()))))
        hit += s > np.percentile(sh, 99)
    return hit / len(d["rows"])


# ---------------------------------------------------------------- drawing
def _w(fig, t): return t.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi * 25.4


def draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None, data=None, order=("cen", "5S")):
    data = data or build()
    ax = ST.ax_mm(fig, left_mm, top_mm, w_mm, h_mm); ax.set_xlim(0, w_mm); ax.set_ylim(h_mm, 0); ax.axis("off")
    if letter: ST.letter(fig, left_mm, top_mm, letter)
    gap = 2.0; colw = (w_mm - gap) / 2; Wmax = max(d["P"] for d in data.values())
    LH, STRIP, GS, MH, BR = 5.6, 1.9, 0.5, 10.5, 3.4          # label block, consensus strip, gap, matrix height, bracket zone (all mm)
    y_head = 5.4; y_top = 8.0; rowh = LH + STRIP + GS + MH + BR; rows_out = []
    grp_name = {"cen": "Centromeric families", "5S": "5S-related families"}
    fams = {"cen": ("GarSat03", "MacSat03"), "5S": ("GarSat04", "MacSat06")}
    for k, grp in enumerate(order):
        x0 = k * (colw + gap)
        ax.text(x0 + 0.3, y_head, grp_name[grp], fontsize=7, fontweight="bold", color=ST.GROUP[grp], ha="left", va="center")
        for j, f in enumerate(fams[grp]):
            d = data[f]; y = y_top + j * (rowh + 0.4); gene = d["unit_class"] == "gene_array"
            t_f = ax.text(x0 + 0.3, y + 1.4, f, fontsize=7, fontweight="bold", color=ST.family_colour(f), ha="left", va="center")
            nx = x0 + 0.3 + _w(fig, t_f) + 1.2
            t_sp = ax.text(nx, y + 1.5, SPNAME[d["h"]], fontsize=6, style="italic", ha="left", va="center")
            ax.text(nx + _w(fig, t_sp) + 1.0, y + 1.5, f"LG{d['c'][-2:]}", fontsize=6, ha="left", va="center")
            what = "5S rDNA array" if gene else ("5S-derived satellite" if grp == "5S" else "")
            sub = (f"{what} \u00b7 " if what else "") + f"{d['P']} bp \u00b7 rows {d['row_id'] * 100:.0f} % identical to consensus"
            t_s = ax.text(x0 + 0.3, y + 4.3, sub, fontsize=6, color=ST.NAVY, ha="left", va="center")
            if _w(fig, t_s) > colw - 0.6: t_s.set_text(sub.replace(" to consensus", ""))
            if _w(fig, t_s) > colw - 0.6: t_s.set_text(sub.replace(" to consensus", "").replace("5S-derived satellite", "5S-derived"))
            wd = (colw - 0.6) * d["P"] / Wmax; xs = x0 + 0.3
            cons = np.array([[IDX.get(c, 4) for c in d["cons"]]])
            ax.imshow(cons, cmap=CMAP, vmin=0, vmax=4, extent=(xs, xs + wd, y + LH + STRIP, y + LH), aspect="auto", interpolation="nearest", zorder=2, rasterized=True)
            M = np.array([[IDX.get(c, 4) for c in r_.ljust(d["P"], "N")[:d["P"]]] for r_ in d["rows"]])
            ym = y + LH + STRIP + GS
            ax.imshow(M, cmap=CMAP, vmin=0, vmax=4, extent=(xs, xs + wd, ym + MH, ym), aspect="auto", interpolation="nearest", zorder=2, rasterized=True)
            ax.set_xlim(0, w_mm); ax.set_ylim(h_mm, 0)
            if grp == "5S":
                a5, b5 = d["five_raw"][:2]; P = d["P"]
                segs = [(a5, min(b5, P))] + ([(0, b5 - P)] if b5 > P else [])
                yb = ym + MH + 0.9
                for (s0, s1) in segs:
                    xa, xb = xs + wd * s0 / P, xs + wd * s1 / P
                    ax.plot([xa, xa, xb, xb], [yb - 0.7, yb, yb, yb - 0.7], color=ST.NAVY, lw=0.5, solid_capstyle="butt", zorder=5)
                    ax.add_patch(Rectangle((xa, y + LH - 0.15), xb - xa, STRIP + GS + MH + 0.3 + 0.15 + 0.0, fill=False, edgecolor=ST.RTE, linewidth=1.0, zorder=6, joinstyle="miter"))   # yellow outline box around the 5S-related segment (consensus strip + all rows)
                s0, s1 = max(segs, key=lambda s: s[1] - s[0]); xc = xs + wd * (s0 + s1) / 2 / P
                tw = 1.25 * len(f"5S-related segment, {b5 - a5} bp"); xt = min(max(xc, x0 + 0.3 + tw / 2), x0 + colw - tw / 2)
                ax.text(xt, yb + 1.6, f"5S-related segment, {b5 - a5} bp", fontsize=6, ha="center", va="center")
            rows_out.append(dict(panel="b", group=grp, family=f, species=d["h"], chrom=d["c"], unit_class_v47=d["unit_class"], element_class_v47=d["element_class"], seq_cluster_v47=d["seq_cluster"],
                                 monomer_bp_stack_period=d["P"], period_from_v453_table=int(d["period_from_v453"]), v47_unit_bp=d["v47_unit"],
                                 run_bp=d["run_bp"], run_start0=d["run_start0"], array_start0_shown=d["a"], rows_shown=len(d["rows"]),
                                 row_identity_to_consensus_pct=round(d["row_id"] * 100, 1), phase_identity_to_Gar_consensus_pct=round(d.get("phase_id", 1) * 100, 1),
                                 five_S_segment_start_bp=d["five_raw"][0] if grp == "5S" else "", five_S_segment_end_bp=d["five_raw"][1] if grp == "5S" else "",
                                 five_S_local_score_on_2x_consensus=d["five_raw"][2], stack_reverse_complemented_vs_earlier_panel=int(d.get("flipped", False)) if grp == "5S" else ""))
    # nucleotide key
    ky = h_mm - 1.6; kx = 0.3
    for c in "ACGT":
        ax.add_patch(Rectangle((kx, ky - 0.9), 1.8, 1.8, facecolor=ST.NUC[c], edgecolor="none", zorder=3))
        ax.text(kx + 2.5, ky, c, fontsize=6, ha="left", va="center"); kx += 5.2
    ax.text(kx + 1.0, ky, "rows = consecutive monomers; columns = aligned positions", fontsize=6, ha="left", va="center", color=ST.NAVY)
    return dict(rows=rows_out)


def write_tables(num, sweep=None, stem="fig1_panel_b"):
    with open(HERE / f"{stem}_numbers.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(num["rows"][0]), delimiter="\t"); w.writeheader(); w.writerows(num["rows"])
    if sweep:
        with open(HERE / f"{stem}_sweep.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(sweep[0]), delimiter="\t"); w.writeheader(); w.writerows(sweep)


def run_sweep():
    """one-at-a-time sweep of the three choices of the view: rows shown (20 / 30 / 40), phase-correction window (+-10 / 15 / 20 bp), start offset into the run (0 / 3000 / 6000 bp)"""
    out = []
    for nrow, shift, skip in [(30, 15, 3000), (20, 15, 3000), (40, 15, 3000), (30, 10, 3000), (30, 20, 3000), (30, 15, 0), (30, 15, 6000)]:
        d = build(nrow, shift, skip)
        for f in ("GarSat03", "MacSat03", "GarSat04", "MacSat06"):
            x = d[f]; out.append(dict(rows=nrow, shift_bp=shift, start_offset_bp=skip, family=f, chrom=x["c"], row_identity_pct=round(x["row_id"] * 100, 1), phase_identity_to_Gar_pct=round(x.get("phase_id", 1) * 100, 1),
                                      five_S_start=x["five_raw"][0] if x["grp"] == "5S" else "", five_S_end=x["five_raw"][1] if x["grp"] == "5S" else ""))
    return out


if __name__ == "__main__":
    ST.house_style(); W, H = 88, 56
    data = build(); fig = ST.new_figure(W, H); num = draw(fig, 0, 0, W, H, "b", data=data)
    for r_ in num["rows"]:
        f = r_["family"]; r_["rows_with_5S_hit_above_shuffled_p99_pct"] = round(100 * five_s_null(data[f]), 0)
    ST.save(fig, str(HERE / "_standalone" / "fig1_panel_b"))
    sweep = run_sweep() if "--nosweep" not in sys.argv else None
    write_tables(num, sweep)
    with open(HERE / "fig1_panel_b_candidates.tsv", "w", newline="") as fh:
        cd = candidates(); w = csv.DictWriter(fh, fieldnames=list(cd[0]), delimiter="\t"); w.writeheader(); w.writerows(cd)
    for r_ in num["rows"]: print(r_)
