"""Shared satellite-block definition for steps 21-22 (2026-10-23 fix).
Block = longest run of MacSat03 monomers in the unit with gaps <= 50 bp, excluding monomers that overlap element hits (the element's built-in
~2-monomer satellite stretch); hits whose query lies only inside that stretch (2880-3400) are ignored, otherwise they mask the whole block. Edge segments are positional, not monomer calls: EDGE bp right after the block start and right before the block end,
taken in unit orientation (block start = side next to the preceding element's 3' end)."""
EDGE = 448
# element consensus carries its own MacSat03 stretch; hits confined to it are satellite matches, not element (fix 2026-10-23b)
BUILTIN_SAT = (2880, 3400)


def block_and_edges(mono_u, he_u, strand, get, chrom, rc):
    m = mono_u.sort_values("genomic_start_0based")
    if len(he_u):
        qlo = he_u[["q0", "q1"]].min(axis=1); qhi = he_u[["q0", "q1"]].max(axis=1)
        he_u = he_u[~((qlo >= BUILTIN_SAT[0]) & (qhi <= BUILTIN_SAT[1]))]
    if len(he_u):
        iv = list(zip(he_u.lo, he_u.hi))
        m = m[[not any(a < e and b > s for s, e in iv) for a, b in zip(m.genomic_start_0based, m.genomic_end)]]
    if m.empty: return "", "", "", 0
    runs, cur = [], [m.iloc[0]]
    for r in m.iloc[1:].itertuples():
        if r.genomic_start_0based - cur[-1].genomic_end <= 50: cur.append(r)
        else: runs.append(cur); cur = [r]
    runs.append(cur)
    best = max(runs, key=lambda x: x[-1].genomic_end - x[0].genomic_start_0based)
    b0, b1 = best[0].genomic_start_0based, best[-1].genomic_end
    seq = get(chrom, b0, b1); seq = seq if strand == "+" else rc(seq)
    return seq, seq[:EDGE], seq[-EDGE:], len(best)


def window_hits(h, w0, w1, min_ov=20, sat_filter=True):
    """Element hits used for one consensus window (fix 2026-10-24): drop hits confined to the built-in satellite stretch (they are
    satellite monomers), then keep non-overlapping hits in query space, best bit score first. Before this, all hits were concatenated,
    so windows near 2,900-3,400 became multi-kb chimeras of satellite monomers and internal-repeat hits."""
    q_lo = h[["q0", "q1"]].min(axis=1); q_hi = h[["q0", "q1"]].max(axis=1)
    if sat_filter: h = h[~((q_lo >= BUILTIN_SAT[0]) & (q_hi <= BUILTIN_SAT[1]))]
    taken, out = [], []
    for r in h.sort_values("bits", ascending=False).itertuples():
        a, b = max(r.q0, w0), min(r.q1, w1)
        if b - a < min_ov: continue
        if any(a < y - 10 and b > x + 10 for x, y in taken): continue
        taken.append((a, b)); out.append(r)
    return out
