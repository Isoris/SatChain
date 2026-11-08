"""Consensus of the dispersed element, its boundaries (conservation profile), ORFs, TSDs and the terminal sequence (2026-10-24)."""
import pathlib, collections, re, numpy as np
from Bio.Seq import Seq
W = pathlib.Path("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/seq")
def read(p):
    d = {}; k = None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"): k = l[1:]; d[k] = []
        else: d[k].append(l.upper())
    return {k: "".join(v) for k, v in d.items()}
for name in ("Gar_full", "Mac_full", "all_full"):
    A = read(W / f"{name}.aln.fa"); M = np.array([list(s) for s in A.values()]); n, L = M.shape
    occ = (M != "-").mean(0)
    cons, ident = [], []
    for j in range(L):
        c = collections.Counter(x for x in M[:, j] if x != "-")
        b, f = c.most_common(1)[0] if c else ("-", 0); ident.append(f / n)
        if occ[j] >= 0.5: cons.append(b)
    cons = "".join(cons); ident = np.array(ident)
    # body = columns where >= 60% of copies share the base; boundaries = first/last 50-col window with mean identity >= 0.6
    w = np.convolve(ident, np.ones(50) / 50, "same"); body = np.where(w >= 0.6)[0]
    c0, c1 = body[0], body[-1]
    cons_body = "".join(collections.Counter(x for x in M[:, j] if x != "-").most_common(1)[0][0] for j in range(c0, c1 + 1) if occ[j] >= 0.5)
    open(W / f"{name}_consensus.fa", "w").write(f">{name}_body\n{cons_body}\n")
    print(f"\n== {name}: {n} copies, alignment {L} cols, body cols {c0}-{c1}, consensus body {len(cons_body)} bp")
    for strand, s in (("+", cons_body), ("-", str(Seq(cons_body).reverse_complement()))):
        orfs = []
        for fr in range(3):
            p = str(Seq(s[fr:len(s) - (len(s) - fr) % 3]).translate())
            for m in re.finditer(r"[^*]{100,}", p): orfs.append((len(m.group()), fr, m.start() * 3 + fr, m.group()))
        orfs.sort(reverse=True)
        print(f"  strand {strand}: ORFs >= 100 aa:", [(a, f"{b}:{c}") for a, b, c, _ in orfs[:5]])
        with open(W / f"{name}_orfs_{'plus' if strand == '+' else 'minus'}.faa", "w") as fh:
            for a, b, c, p in orfs: fh.write(f">{name}_{strand}_f{b}_{c}_{a}aa\n{p}\n")
    print("  first 60 bp of body (consensus orientation):", cons_body[:60])
    print("  last 60 bp of body:", cons_body[-60:])
    # TSDs per copy: unaligned genomic flanks around the body columns
    tsd = []
    for k, s in A.items():
        left = "".join(x for x in s[:c0] if x != "-")[-25:]; right = "".join(x for x in s[c1 + 1:] if x != "-")[:25]
        best = ""
        for ln in range(25, 4, -1):
            for i in range(0, 25 - ln + 1):
                sub = left[i:i + ln]
                if len(sub) == ln and sub in right[:ln + 10] and i >= 25 - ln - 10: best = sub; break
            if best: break
        tsd.append(len(best))
    t = np.array(tsd); print("  TSD-like direct repeats at both ends (>=5 bp, within 10 bp of the boundaries): copies", int((t >= 5).sum()), "of", len(t), "| lengths", collections.Counter(t[t >= 5]).most_common(6))
