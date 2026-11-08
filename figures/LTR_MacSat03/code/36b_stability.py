"""Step 36b (2026-10-24): resolvability / stability score per MacSat03 centromere (SatChain v4.6 runs + interruptions).
1 copy confidence: share of step-03 copies inside v4.6 MacSat03 runs with identity >= 75 (confident), 60-75 (weak), < 60 (junk).
2 boundary stability: per junction, distance between the satellite run edge and the RTE fragment end (split/v46_fragments.tsv);
  share of units within +-tol bp of the modal distance (tol 10/25/50), left and right junction.
3 positional identity: copies indexed 1..P from each junction of a run (P = 6), in run orientation. Identity proxy = global alignment score
  (match 1, mismatch 0, gap -1) / max length. Same-position identity (across units) vs adjacent-position identity (within a unit);
  leave-one-out nearest-position classification accuracy (chance = 1/P); 200 label permutations for a p-value."""
import pathlib, collections, random, pandas as pd, numpy as np
from Bio import Align
exec(open(pathlib.Path(__file__).parent / "32e_extract.py").read().split("random.seed")[0])
random.seed(7); fa = Fa(G["Mac"]); D = pathlib.Path("F:/data/fClaHyb_func_annotation_v1.0/fClaHyb_Mac/fClaHyb_Mac_satelliteDNA")
C = pd.read_csv(D / "fClaHyb_Mac.satDNA_v46.catalogue.tsv", sep="\t"); C = C[C.family == "MacSat03"]
cop = pd.read_csv("F:/data/satellite_FISH/pipeline_v4_5_fixed/out/fClaHyb_Mac_satellites/03_monomer_copies.tsv", sep="\t",
                  usecols=["chrom", "genomic_start_0based", "genomic_end", "strand", "best_alignment_identity", "observed_length_bp"])
cop = cop[cop.chrom.isin(C.chrom.unique())]
F = pd.read_csv(W / "split" / "v46_fragments.tsv", sep="\t"); F = F[F.aligned_bp >= 100]
GA = Align.PairwiseAligner(mode="global", match_score=1, mismatch_score=0, open_gap_score=-1, extend_gap_score=-1)
P = 6
def ident(a, b): return GA.score(a, b) / max(len(a), len(b))
rows, posrows = [], []
for ch, runs in C.groupby("chrom"):
    lab = ch.replace("C_mac_", ""); runs = runs.sort_values("start0")
    cc = cop[cop.chrom == ch]
    # 1 copy confidence inside runs
    ins = np.zeros(len(cc), bool)
    for r in runs.itertuples(): ins |= ((cc.genomic_start_0based >= r.start0) & (cc.genomic_end <= r.end)).values
    ci = cc[ins]; idt = ci.best_alignment_identity
    rec = dict(chrom=lab, runs=len(runs), copies=len(ci), conf=(idt >= 75).mean(), weak=((idt >= 60) & (idt < 75)).mean(), junk=(idt < 60).mean())
    # 2 boundary stability
    f = F[F.chrom == lab]; rec["units_with_RTE"] = len(f)
    for side in ("left_flank_bp", "right_flank_bp"):
        x = f[side]
        for tol in (10, 25, 50):
            rec[f"{side.split('_')[0]}_junction_pm{tol}"] = ((x - x.mode().iloc[0]).abs() <= tol).mean() if len(x) else np.nan
    # 3 positional identity: copies of runs with >= 2P copies, positions from both ends
    S = {}
    sat = [r for r in runs.itertuples() if r.n_conf >= 2 * P]
    plus = (ci.strand == "+").mean() >= 0.5
    for u, r in enumerate(sat):
        c = ci[(ci.genomic_start_0based >= r.start0) & (ci.genomic_end <= r.end)].sort_values("genomic_start_0based")
        if len(c) < 2 * P: continue
        seqs = [fa.get(ch, a, b) for a, b in zip(c.genomic_start_0based, c.genomic_end)]
        if not plus: seqs = [rc(s) for s in seqs[::-1]]
        for p in range(P): S[(u, "L", p + 1)] = seqs[p]; S[(u, "R", p + 1)] = seqs[-1 - p]
    units = sorted({k[0] for k in S}); rec["units_for_position_test"] = len(units)
    if len(units) >= 6:
        units = random.sample(units, min(20, len(units)))
        for end in ("L", "R"):
            keys = [(u, end, p) for u in units for p in range(1, P + 1)]
            M = np.array([[ident(S[a], S[b]) if a != b else np.nan for b in keys] for a in keys])
            lab_p = np.array([k[2] for k in keys]); lab_u = np.array([k[0] for k in keys])
            same = np.nanmean([M[i, j] for i in range(len(keys)) for j in range(len(keys)) if i != j and lab_p[i] == lab_p[j] and lab_u[i] != lab_u[j]])
            adj = np.nanmean([M[i, j] for i in range(len(keys)) for j in range(len(keys)) if lab_u[i] == lab_u[j] and abs(lab_p[i] - lab_p[j]) == 1])
            def acc(lp):
                ok = 0
                for i in range(len(keys)):
                    sc = {p: np.nanmean([M[i, j] for j in range(len(keys)) if lp[j] == p and lab_u[j] != lab_u[i]]) for p in range(1, P + 1)}
                    ok += max(sc, key=sc.get) == lp[i]
                return ok / len(keys)
            a0 = acc(lab_p)
            null = []
            for _ in range(200):
                lp = lab_p.copy()
                for u in set(lab_u): idx = np.where(lab_u == u)[0]; lp[idx] = np.random.permutation(lp[idx])
                null.append(acc(lp))
            rec[f"{end}_same_pos_id"] = same; rec[f"{end}_adjacent_id"] = adj; rec[f"{end}_pos_accuracy"] = a0
            rec[f"{end}_pos_p"] = (1 + sum(n >= a0 for n in null)) / 201
            for p in range(1, P + 1):
                ii = np.where(lab_p == p)[0]
                posrows.append(dict(chrom=lab, end=end, position=p, same_pos_id=np.nanmean(M[np.ix_(ii, ii)])))
    rows.append(rec); print(lab, "done", flush=True)
R = pd.DataFrame(rows); R.to_csv(W / "split" / "stability_scores.tsv", sep="\t", index=False); pd.DataFrame(posrows).to_csv(W / "split" / "stability_positions.tsv", sep="\t", index=False)
pd.set_option("display.width", 300); print(R.round(2).to_string(index=False))
