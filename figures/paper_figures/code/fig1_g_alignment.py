"""Figure 1 panel g: align the 5S rRNA gene to the consensus monomers of the 5S-related families (GarSat04 = 5S-derived satellite; GarSat49, MacSat06 = 5S rDNA gene arrays) and,
as the control, the centromeric GarSat03 / MacSat03 group.  (Variable names LB / LA = the two sets; the word "lineage" is no longer used in any figure text, SatChain v4.7, 2026-11-01.)
Biopython PairwiseAligner, local, both strands of the monomer, circular monomer (a hit may wrap over the seam).  No HMMER, no WSL.  Windows python, 1 process, small memory.
 inputs : 16 annotated 5S genes (barrnap, rRNA.sorted.gff3.gz = the 16 '5S rRNA reference' sequences of the lineage B tree, taken UNROTATED from the genome),
          family consensus monomers (v4.6 FISH_package families/*_consensus.fa; GarSat49 and the v4.5 set from pipeline_v4_5_fixed/out/*/04_satellite_family_consensuses.fa)
 outputs: panels/fig1_panel_g_data/  refs_5S.fa, monomer_consensus.fa, best_alignments.tsv, sweep.tsv, null.tsv, copies_check.tsv, summary.json
usage: python -I fig1_g_alignment.py"""
import gzip, json, random, itertools, pathlib, statistics as stt
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
OUT = pathlib.Path(__file__).resolve().parent.parent / "panels" / "fig1_panel_g_data"
SF = pathlib.Path(r"F:\data\satellite_FISH")
GEN = pathlib.Path(r"F:\data\fClaHyb_func_annotation_v1.0")
rc = lambda s: str(Seq(s).reverse_complement())
LB = ["GarSat04", "GarSat49", "MacSat06"]; LA = ["GarSat03", "MacSat03"]
FAMDIR = {"GarSat04": "GarSat04_CgaSat001-216", "MacSat06": "MacSat06_CmaSat002-214", "GarSat03": "GarSat03_CgaSat002-226", "MacSat03": "MacSat03_CmaSat001-225"}


def fasta(p):
    d, k = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"): k = l[1:]; d[k] = ""
        elif k: d[k] += l.upper()
    return d


def genome_seq(sp, chrom, s, e):
    p = GEN / f"fClaHyb_{sp}" / f"fClaHyb_{sp}_genome" / f"fClaHyb_{sp}_LG.fa"
    f = {l.split()[0]: l.split()[1:] for l in open(str(p) + ".fai")}
    L, off, bl, bw = map(int, f[chrom])
    with open(p, "rb") as h:
        a = off + s // bl * bw + s % bl; b = off + (e - 1) // bl * bw + (e - 1) % bl + 1
        h.seek(a); return h.read(b - a).decode().replace("\n", "").upper()


def refs_5S():
    names = [l[1:].strip() for l in open(SF / "share_LITE_v4.5_2026-10-15/phylogeny/lineageB_5S_derived/lineageB_5S_derived_copies.fa") if l.startswith(">") and "5S_rRNA" in l]
    rows = {}
    for sp in ("Gar", "Mac"):
        with gzip.open(SF / f"annotation_overlays/{sp}/fClaHyb_{sp}.rRNA.sorted.gff3.gz", "rt") as fh:
            rows[sp] = [l.rstrip("\n").split("\t") for l in fh if "5S" in l and not l.startswith("#")]
    out = {}
    for n in names:
        sp, _, lg, st = n.split("|")
        r = [x for x in rows[sp] if x[0].endswith(lg) and int(x[3]) == int(st)][0]
        s = genome_seq(sp, r[0], int(r[3]) - 1, int(r[4]))
        out[f"{sp}|{lg}|{r[3]}"] = s if r[6] == "+" else rc(s)
    return out


def consensus_sets():
    fam = {}
    for f, d in FAMDIR.items():
        s = fasta(SF / f"runs_of_record_v4.6_2026-10-27/FISH_package/families/{d}/{f}_consensus.fa")
        fam[f] = s[[x for x in s if "family_consensus" in x][0]]
    v45 = {}
    for sp in ("Gar", "Mac"):
        s = fasta(SF / f"pipeline_v4_5_fixed/out/fClaHyb_{sp}_satellites/04_satellite_family_consensuses.fa")
        for k, v in s.items():
            if k.split()[0].split("|")[0] in LB + LA: v45[k.split()[0].split("|")[0]] = v
    fam["GarSat49"] = v45["GarSat49"]          # GarSat49 is not in the v4.6 FISH package: v4.5 consensus (209 bp)
    return fam, v45


def mk(match=2, mism=-3, go=-5, ge=-2):
    a = PairwiseAligner(mode="local"); a.match_score = match; a.mismatch_score = mism; a.open_gap_score = go; a.extend_gap_score = ge
    return a


def best_hit(al, gene, mono, want_cols=False):
    """best local alignment of the gene to the monomer, both strands, monomer circular (mono + its first len(gene) bases).
    Returns coordinates on the + strand of the consensus as given: m0, m1 (0-based, half-open, m1 may exceed len when the hit wraps), gene g0, g1."""
    L = len(mono); best = None
    for strand, q in (("+", mono), ("-", rc(mono))):
        qq = q + q[:len(gene)]
        a = al.align(qq, gene)[0]
        if best is None or a.score > best[0].score: best = (a, strand, qq)
    a, strand, qq = best
    mb, gb = a.aligned[0], a.aligned[1]                 # target = monomer (first arg), query = gene
    m0, m1 = int(mb[0][0]), int(mb[-1][1]); g0, g1 = int(gb[0][0]), int(gb[-1][1])
    cols = []; mt = 0; n = 0
    for (ma, mbb), (ga, gbb) in zip(mb, gb):
        for i in range(mbb - ma):
            ok = qq[ma + i] == gene[ga + i]; cols.append((ma + i, ga + i, ok)); mt += ok; n += 1
    for k in range(1, len(mb)):
        n += (mb[k][0] - mb[k - 1][1]) + (gb[k][0] - gb[k - 1][1])
    wrap = m1 > L
    if strand == "+": ms, me = m0, m1
    else: ms, me = L - m1, L - m0                       # (may be negative when wrapped)
    if strand == "-":
        cols = [((L - 1 - p) % L, g, ok) for p, g, ok in cols]
    else:
        cols = [(p % L, g, ok) for p, g, ok in cols]
    return dict(score=float(a.score), strand=strand, m0=ms % L, m1=ms % L + (m1 - m0), span_m=m1 - m0, g0=g0, g1=g1, span_g=g1 - g0, aln_len=int(n), matches=int(mt),
                ident=float(mt / n if n else 0), wrap=wrap, cols=cols if want_cols else None)


def shuffled(s, rnd):
    l = list(s); rnd.shuffle(l); return "".join(l)


def main():
    OUT.mkdir(exist_ok=True)
    refs = refs_5S()
    with open(OUT / "refs_5S.fa", "w") as fh:
        for k, v in refs.items(): fh.write(f">{k} len={len(v)}\n{v}\n")
    fam, v45 = consensus_sets()
    with open(OUT / "monomer_consensus.fa", "w") as fh:
        for k, v in fam.items(): fh.write(f">{k}|set_v46{'_(GarSat49_from_v4.5)' if k == 'GarSat49' else ''} len={len(v)}\n{v}\n")
        for k, v in v45.items(): fh.write(f">{k}|set_v45 len={len(v)}\n{v}\n")
    al0 = mk(); names = list(refs)
    mean = {a: stt.mean(al0.score(refs[a], refs[b]) for b in names if b != a) for a in names}
    ref_name = max(mean, key=mean.get); ref = refs[ref_name]
    idm = []
    for a, b in itertools.combinations(names, 2):
        idm.append(best_hit(al0, refs[a], refs[b])["ident"])
    summ = dict(n_refs=len(refs), ref_lengths=[len(v) for v in refs.values()], ref_medoid=ref_name, ref_len=len(ref), ref_seq=ref,
                refs_pairwise_ident_median=stt.median(idm), refs_pairwise_ident_min=min(idm))
    rows = []
    for setname, cs in (("v46", fam), ("v45", v45)):
        for f, m in cs.items():
            for rn, rs in refs.items():
                h = best_hit(al0, rs, m)
                rows.append(dict(set=setname, family=f, mono_len=len(m), ref=rn, ref_len=len(rs), score=h["score"], strand=h["strand"], mono_start=h["m0"], mono_end=h["m1"], mono_span=h["span_m"],
                                 gene_start=h["g0"], gene_end=h["g1"], gene_span=h["span_g"], aln_len=h["aln_len"], matches=h["matches"], identity=round(h["ident"], 4), wraps_seam=h["wrap"]))
    pd.DataFrame(rows).to_csv(OUT / "best_alignments.tsv", sep="\t", index=False)
    det = {}
    for f in LB + LA:
        h = best_hit(al0, ref, fam[f], True)
        det[f] = dict(strand=h["strand"], m0=h["m0"], m1=h["m1"], g0=h["g0"], g1=h["g1"], score=h["score"], ident=h["ident"], aln_len=h["aln_len"], wrap=h["wrap"], mono_len=len(fam[f]),
                      cols=[(int(a), int(b), bool(c)) for a, b, c in h["cols"]])
    summ["detail"] = det
    sw = []
    for (ma, mi, go, ge) in itertools.product((1, 2, 3), (-1, -2, -3, -4), (-3, -5, -8), (-1, -2)):
        al = mk(ma, mi, go, ge)
        for f in LB + LA:
            h = best_hit(al, ref, fam[f])
            sw.append(dict(match=ma, mismatch=mi, gap_open=go, gap_extend=ge, family=f, score=h["score"], strand=h["strand"], mono_start=h["m0"], mono_end=h["m1"], mono_span=h["span_m"],
                           gene_start=h["g0"], gene_end=h["g1"], aln_len=h["aln_len"], identity=round(h["ident"], 4)))
    pd.DataFrame(sw).to_csv(OUT / "sweep.tsv", sep="\t", index=False)
    rnd = random.Random(20261027); nul = []
    for f in LB + LA:
        m = fam[f]; real = best_hit(al0, ref, m)["score"]
        sc = [best_hit(al0, ref, shuffled(m, rnd))["score"] for _ in range(1000)]
        mu = stt.mean(sc); sd = stt.pstdev(sc)
        nul.append(dict(family=f, mono_len=len(m), real_score=real, null_mean=round(mu, 2), null_sd=round(sd, 2), null_max=max(sc), null_q99=sorted(sc)[989], n_ge_real=sum(x >= real for x in sc), n_shuffles=1000,
                        z=round((real - mu) / sd, 2) if sd else None, p_emp=(sum(x >= real for x in sc) + 1) / 1001))
    pd.DataFrame(nul).to_csv(OUT / "null.tsv", sep="\t", index=False)
    cp = []
    for f, d in FAMDIR.items():
        mons = fasta(SF / f"runs_of_record_v4.6_2026-10-27/FISH_package/families/{d}/{f}_monomers.fa")
        keys = list(mons); random.Random(1).shuffle(keys)
        for k in keys[:200]:
            h = best_hit(al0, ref, mons[k])
            cp.append(dict(family=f, monomer=k, mono_len=len(mons[k]), score=h["score"], aln_len=h["aln_len"], identity=round(h["ident"], 4), span_m=h["span_m"], strand=h["strand"], mono_start=h["m0"], mono_end=h["m1"]))
    pd.DataFrame(cp).to_csv(OUT / "copies_check.tsv", sep="\t", index=False)
    json.dump(summ, open(OUT / "summary.json", "w"), indent=1)
    print("done", ref_name, len(ref))


if __name__ == "__main__":
    main()
