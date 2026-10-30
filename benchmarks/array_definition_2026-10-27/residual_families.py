"""Unannotated residuals of the major centromeres: extract, cluster into families, profile (2026-10-27). Windows Python only (no WSL).
1 Segments = parts of each major centromere span that are neither satellite run nor covered by RTE / SINE / TE whole element / tRNA / rRNA /
  Rfam ncRNA / exon / other satellite family (definitions of residuals.py); >= 100 bp kept, >= 300 bp clustered.
2 Families: canonical 15-mer MinHash sketch (s = 300) per segment; containment = shared / min(sketch sizes); single-linkage components at
  containment >= c (c swept 0.2 / 0.3 / 0.5; 0.3 reported).
3 Profile per family (representative = member with the most links, ties -> longer): members, centromeres, length median, GC, microsatellite
  share (1-6 bp motifs >= 16 bp), terminal direct / inverted repeats (>= 12-mer shared between the first and last 200 bp), longest ORF (aa) over
  6 frames, genome copies = loci where >= 30 % of the representative's non-microsatellite 21-mers fall within a window of its length
  (21-mer scan of both genome FASTAs), split into inside / outside the 54 major centromere spans.
Outputs: residual_segments.tsv / .fa, residual_families.tsv, residual_family_members.tsv, residual_families.png
usage: python -P residual_families.py"""
import os, re, sys, collections, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
B = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, B)
exec(open(os.path.join(B, "residuals.py"), encoding="utf-8").read().split("rows, feats = [], []")[0])   # loaders: find, gff, union, clip, bp, I, RTE, R0
from kmer_homogeneity import Fasta, kcodes
BS = chr(92)
def microsat_mask(s):
    m = np.zeros(len(s), dtype=bool); t = s.decode()
    for k in range(1, 7):
        for x in re.finditer(("(([ACGT]{%d})" + BS + "2{%d,})") % (k, max(2, 12 // k)), t):
            if x.end() - x.start() >= 16: m[x.start():x.end()] = True
    return m
segs = []
for sp in ("Mac", "Gar"):
    C = pd.read_csv(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_satelliteDNA/fClaHyb_{sp}.satDNA_v46.catalogue.tsv", sep="\t"); fa = Fasta(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.fa")
    src = {k: gff(find(sp, *v), t) for k, (v, t) in {"tRNA": (("tRNA_tRNAscan", "*.sorted.gff3.gz"), {"tRNA"}), "rRNA": (("rRNA_barrnap", "*.rRNA.sorted.gff3.gz"), {"rRNA"}),
           "ncRNA": (("ncRNA_rfam", "*.rfam.all_classes.sorted.gff3.gz"), None), "exon": (("genes_lifton", "*.lifton.sorted.gff3.gz"), {"exon", "CDS"}),
           "SINE": (("repeats_AnnoSINE", "*.SINE.AnnoSINE.annotations.gff3.gz"), None)}.items()}
    te = gff(find(sp, "repeats_EDTA", f"fClaHyb_{sp}.TEanno.sorted.gff3.gz"))
    for r in I[I.species == sp].itertuples():
        x = C[C.array_id == r.array_id]; c = r.chrom; lo, hi = int(x.start0.min()), int(x.end.max())
        runs = union(list(zip(x.start0, x.end))); non = []; cur = lo
        for a, b in runs:
            if a > cur: non.append((cur, a))
            cur = max(cur, b)
        cov = []
        for d in src.values(): cov += [(a, b) for a, b, *_ in d.get(c, []) if b > lo and a < hi]
        cov += [(a, b) for a, b in zip(RTE[RTE.chrom == c].a, RTE[RTE.chrom == c].b) if b > lo and a < hi]
        cov += [(a, b) for a, b, t, *_ in te.get(c, []) if b > lo and a < hi and b - a >= 500 and bp(clip([(a, b)], non)) >= 0.8 * (b - a) and t not in ("repeat_region", "target_site_duplication", "long_terminal_repeat")]
        oth = C[(C.chrom == c) & (C.array_id != r.array_id) & (C.end > lo) & (C.start0 < hi)]; cov += list(zip(oth.start0, oth.end))
        cov = clip(cov, non)
        for a, b in non:                                                       # non-satellite minus covered
            cur = a
            for ca, cb in [iv for iv in cov if iv[1] > a and iv[0] < b]:
                if ca > cur: segs.append((sp, c, r.array_id, r.family, cur, ca))
                cur = max(cur, cb)
            if b > cur: segs.append((sp, c, r.array_id, r.family, cur, b))
    seqmap = {}
S = pd.DataFrame(segs, columns=["species", "chrom", "array_id", "family", "start", "end"]); S["bp"] = S.end - S.start; S = S[S.bp >= 100].reset_index(drop=True)
FA = {sp: Fasta(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_genome/fClaHyb_{sp}_LG.fa") for sp in ("Mac", "Gar")}
seqs = [FA[r.species].get(r.chrom, int(r.start), int(r.end)).upper() for r in S.itertuples()]   # genome is soft-masked
S["GC"] = [round(sum(ch in b"GC" for ch in q) / max(1, len(q)) * 100, 1) for q in seqs]; S["microsat_share"] = [round(float(microsat_mask(q).mean()), 3) for q in seqs]
S["seg_id"] = [f"R{i:05d}" for i in range(len(S))]; S.to_csv(f"{B}/residual_segments.tsv", sep="\t", index=False)
with open(f"{B}/residual_segments.fa", "w", newline="\n") as fh:
    for r, q in zip(S.itertuples(), seqs): fh.write(f">{r.seg_id} {r.species} {r.chrom}:{r.start}-{r.end} {r.family}\n{q.decode()}\n")
print("segments >= 100 bp:", len(S), "| bp:", int(S.bp.sum()), "| by species:", S.groupby("species").bp.sum().to_dict(), flush=True)
# sketches and containment
big = S.index[S.bp >= 300].tolist(); SK = 300; sk = {}
for i in big:
    q = seqs[i]; m = microsat_mask(q); q2 = bytes(b if not mm else 78 for b, mm in zip(q, m))
    c = np.unique(kcodes(q2, 15)); h = np.sort((c * np.uint64(0x9E3779B97F4A7C15)) ^ (c >> np.uint64(29)))[:SK]
    if len(h) >= 30: sk[i] = h
inv = collections.defaultdict(list)
for i, h in sk.items():
    for v in h.tolist(): inv[v].append(i)
shared = collections.Counter()
for lst in inv.values():
    if 1 < len(lst) <= 2000:
        for a in range(len(lst)):
            for b in range(a + 1, len(lst)): shared[(lst[a], lst[b])] += 1
def components(cut):
    par = {i: i for i in sk}
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    deg = collections.Counter()
    for (a, b), n in shared.items():
        if n / min(len(sk[a]), len(sk[b])) >= cut: par[f(a)] = f(b); deg[a] += 1; deg[b] += 1
    comp = collections.defaultdict(list)
    for i in sk: comp[f(i)].append(i)
    return [v for v in comp.values() if len(v) >= 2], deg
SW = {}
for cut in (0.2, 0.3, 0.5):
    comps, _ = components(cut); SW[cut] = comps
    print(f"containment >= {cut}: {len(comps)} families with >= 2 members, covering {sum(len(v) for v in comps)} of {len(sk)} segments >= 300 bp", flush=True)
comps, deg = components(0.3); comps = sorted(comps, key=lambda v: -sum(S.bp[v]))
def orf(q):
    RC = bytes.maketrans(b"ACGTN", b"TGCAN"); best = 0
    for s_ in (q, q.translate(RC)[::-1]):
        for f in range(3):
            run = 0
            for i in range(f, len(s_) - 2, 3):
                if s_[i:i + 3] in (b"TAA", b"TAG", b"TGA"): best = max(best, run); run = 0
                else: run += 1
            best = max(best, run)
    return best
def terminal(q):
    RC = bytes.maketrans(b"ACGTN", b"TGCAN"); a, b = q[:200], q[-200:]; ka = {a[i:i + 12] for i in range(len(a) - 11)}
    direct = any(b[i:i + 12] in ka for i in range(len(b) - 11)); rb = b.translate(RC)[::-1]; inv_ = any(rb[i:i + 12] in ka for i in range(len(rb) - 11))
    return ("direct " if direct else "") + ("inverted" if inv_ else "")
fam_rows, mem_rows = [], []
cen_spans = collections.defaultdict(list)
for sp in ("Mac", "Gar"):
    C = pd.read_csv(f"{R0}/fClaHyb_{sp}/fClaHyb_{sp}_satelliteDNA/fClaHyb_{sp}.satDNA_v46.catalogue.tsv", sep="\t")
    for r in I[I.species == sp].itertuples():
        x = C[C.array_id == r.array_id]; cen_spans[r.chrom].append((int(x.start0.min()), int(x.end.max())))
reps = {}
for k, v in enumerate(comps[:40]):
    rep = max(v, key=lambda i: (deg[i], S.bp[i])); reps[f"F{k + 1:02d}"] = rep
    for i in v: mem_rows.append(dict(residual_family=f"F{k + 1:02d}", **S.loc[i, ["seg_id", "species", "chrom", "array_id", "family", "start", "end", "bp", "GC", "microsat_share"]].rename({"family": "sat_family"}).to_dict()))
# genome copies: 21-mer scan
famk = {}
for f, i in reps.items():
    q = seqs[i]; m = microsat_mask(q); q2 = bytes(b if not mm else 78 for b, mm in zip(q, m)); famk[f] = np.unique(kcodes(q2, 21))
allk = np.unique(np.concatenate(list(famk.values()))); owner = {}
for f, ks in famk.items():
    for v in ks.tolist(): owner.setdefault(v, set()).add(f)
copies = {f: collections.Counter() for f in famk}
for sp in ("Mac", "Gar"):
    fa = FA[sp]
    for c in sorted(fa.idx):
        seq = fa.get(c, 0, fa.idx[c][0]); s_ = np.frombuffer(seq, dtype=np.uint8)
        # positions of k-mers present in any family (chunked to bound memory)
        for off in range(0, len(seq), 20_000_000):
            sub = seq[off:off + 20_000_000 + 20]; codes = kcodes(sub, 21)
            hit = np.where(np.isin(codes, allk))[0]
            if not len(hit): continue
            for f, ks in famk.items():
                hf = hit[np.isin(codes[hit], ks)] + off
                if len(hf) < 5: continue
                L = len(seqs[reps[f]]); need = max(5, int(0.3 * len(ks))); j0 = 0
                loci = []
                for j in range(len(hf)):                                       # sliding window of the representative's length
                    while hf[j] - hf[j0] > L: j0 += 1
                    if j - j0 + 1 >= need and (not loci or hf[j] - loci[-1] > L): loci.append(hf[j])
                for p in loci:
                    inside = any(a - 2000 <= p <= b + 2000 for a, b in cen_spans.get(c, []))
                    copies[f][(sp, "major centromere" if inside else "elsewhere")] += 1
    print("scanned", sp, flush=True)
for k, v in enumerate(comps[:40]):
    f = f"F{k + 1:02d}"; i = reps[f]; d = S.loc[v]
    fam_rows.append(dict(family=f, members=len(v), centromeres=d.species.str.cat(d.chrom.str.split("_").str[-1], sep=" ").nunique(), cen_list=", ".join(sorted(set(d.species + " " + d.chrom.str.split("_").str[-1])))[:120],
                         bp_total=int(d.bp.sum()), len_median=int(d.bp.median()), rep=S.seg_id[i], rep_bp=int(S.bp[i]), GC=float(d.GC.median()), microsat=float(d.microsat_share.median()),
                         terminal_repeats=terminal(seqs[i]) or "none", longest_orf_aa=orf(seqs[i]),
                         **{f"copies_{sp}_{w.replace(' ', '_')}": copies[f][(sp, w)] for sp in ("Mac", "Gar") for w in ("major centromere", "elsewhere")}))
F = pd.DataFrame(fam_rows); F.to_csv(f"{B}/residual_families.tsv", sep="\t", index=False); pd.DataFrame(mem_rows).to_csv(f"{B}/residual_family_members.tsv", sep="\t", index=False)
pd.set_option("display.width", 320); pd.set_option("display.max_colwidth", 60); print(F.to_string(index=False))