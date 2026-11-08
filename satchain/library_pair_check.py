#!/usr/bin/env python3
"""library_pair_check.py  (SatChain v4.7, 2026-11-08): is a 'shared' family pair of the two haplotypes really the same family?

Purpose (declared): QC of the 59h library call. 59h calls two families homologous from ONE alignment (the unit consensus of A vs the
consensus of B at >= 0.80 identity over >= 50 % of the unit). This script checks every named pair at four levels, each against a null,
and draws them side by side; it decides nothing downstream, it tells the author whether the library call holds.

Levels, per pair (A in one haplotype, B in the other):
  1 monomer     unit A vs unit B (rotation-aware local alignment, parasail SW; identity, aligned columns, gaps);
                copy level: up to --n-copies confident copies of A aligned to consensus A and to consensus B, and the same for B
                (A->A, A->B, B->B, B->A identity distributions); null = the same copies aligned to --n-null random consensuses of the
                other haplotype with a unit within +-25 % of the length; second-best: identity of A to its second-best other-haplotype
                family from the 59g hits (the gap between best and second-best says how unique the match is); reciprocal: is A the best
                hit of B?
  2 composition GC, CpG o/e, Jensen-Shannon distance between the canonical 4-mer profiles of the two units; null = JSD of random
                cross-haplotype unit pairs of similar length (percentile reported)
  3 organisation from 59h: arrays, array length median / max, bp share in arrays >= 10 kb, chromosomes, loci, largest-locus share,
                candidate-centromere share, copies, divergence, strand coherence; drawn as array-length histograms and bp per chromosome
  4 topology    53b of the largest run of each family (existing-source runs only): n_states, entropy rate, memory excess, structure
                index, order-2 excess; transition matrices drawn side by side (states rank-named V1.., so the two are comparable)
Verdict per pair: same_family (consensus >= 0.80 AND copies of A are as close to B as to their own consensus within 10 points AND
null percentile >= 99 AND reciprocal), sister_family (consensus 0.70-0.80 or copy gap > 10 points), not_supported (otherwise).
Outputs <dest>/pairs_check.tsv, <dest>/pair_<A>_<B>.png/.pdf (450 dpi), <dest>/log.txt.
usage: library_pair_check.py <final dir> <v4.5 out root> <new03 pattern {G}> <dest> [--pairs A:B,A:B,...] [--top 12] [--n-copies 150] [--n-null 20] [--seed 7]
"""
import sys, os, csv, math, collections, argparse, random, time
import numpy as np, parasail
ap = argparse.ArgumentParser(); ap.add_argument('final'); ap.add_argument('out45'); ap.add_argument('new03'); ap.add_argument('dest')
ap.add_argument('--pairs', default=''); ap.add_argument('--top', type=int, default=12); ap.add_argument('--n-copies', type=int, default=150); ap.add_argument('--n-null', type=int, default=20); ap.add_argument('--seed', type=int, default=7)
a = ap.parse_args(); os.makedirs(a.dest, exist_ok=True); random.seed(a.seed); LOG = open(f'{a.dest}/log.txt', 'a', encoding='utf-8')
try: sys.stdout.reconfigure(encoding='utf-8')     # Windows console is cp1252; the log uses arrows
except Exception: pass
np.seterr(divide='ignore', invalid='ignore')
def P(*x):
    s = ' '.join(str(v) for v in x); print(s, flush=True); LOG.write(s + '\n'); LOG.flush()
MAT = parasail.matrix_create('ACGT', 2, -3); COMP = str.maketrans('ACGTacgt', 'TGCAtgca'); rc = lambda s: s.translate(COMP)[::-1]
def sw(q, t):
    r = parasail.sw_trace_striped_16(q, t, 5, 2, MAT); tb = r.traceback; cols = len(tb.query); m = tb.comp.count('|')
    return dict(score=int(r.score), ident=m / cols if cols else 0.0, cols=cols, gaps=tb.query.count('-') + tb.ref.count('-'), q=tb.query, c=tb.comp, t=tb.ref, tcov=len(tb.ref.replace('-', '')) / len(t))
def best_orient(q, t):
    """unit vs unit, rotation-aware: query doubled, both strands of the target."""
    h1, h2 = sw(q * 2, t), sw(q * 2, rc(t)); h = h1 if h1['score'] >= h2['score'] else h2; h['strand'] = '+' if h is h1 else '-'; return h
def copy_vs_cons(copy, cons):
    h1, h2 = sw(cons * 2, copy), sw(cons * 2, rc(copy)); h = h1 if h1['score'] >= h2['score'] else h2
    return h['ident'] if h['tcov'] >= 0.6 else h['ident'] * h['tcov']     # a copy that aligns over < 60 % of its length is penalised
def kmer_profile(s, k=4):
    s = s.upper(); c = collections.Counter()
    for i in range(len(s) - k + 1):
        w = s[i:i + k]
        if set(w) <= set('ACGT'): c[min(w, rc(w))] += 1
    tot = sum(c.values()) or 1; return {w: v / tot for w, v in c.items()}
def jsd(p, q):
    keys = set(p) | set(q); pa = np.array([p.get(k, 0) for k in keys]); qa = np.array([q.get(k, 0) for k in keys]); m = (pa + qa) / 2
    kl = lambda x, y: float(np.sum(np.where(x > 0, x * np.log2(x / np.where(y > 0, y, 1)), 0)))
    return math.sqrt(0.5 * kl(pa, m) + 0.5 * kl(qa, m))
# ---------------- load
U = {}; MET = {}; HITS = collections.defaultdict(list); CAT = {}; TOPO = {}; TRANS = collections.defaultdict(lambda: collections.Counter())
for G in ('Gar', 'Mac'):
    for r in csv.DictReader(open(f'{a.final}/unit_content_2026-11-07/59g_unit_content_{G}.tsv'), delimiter='\t'): U[(G, r['family'])] = r['unit_seq'].upper();
    for r in csv.DictReader(open(f'{a.final}/family_metrics_2026-11-07/59h_family_metrics_{G}.tsv'), delimiter='\t'): MET[(G, r['family'])] = r
    for h in csv.DictReader(open(f'{a.final}/unit_content_2026-11-07/59g_hits_{G}.tsv'), delimiter='\t'):
        if h['ref_group'] == 'SAT' and float(h['unit_frac']) >= 0.5 and float(h['ref_cov']) >= 0.5: HITS[(G, h['family'])].append((h['ref'], float(h['identity'])))
    CAT[G] = list(csv.DictReader(open(f'{a.final}/cat59/59_catalogue_{G}.tsv'), delimiter='\t'))
    tp = f'{a.final}/cat53/{G}_existing/53b/53b_topologies.tsv'
    if os.path.exists(tp):
        for r in csv.DictReader(open(tp), delimiter='\t'):
            if r['unit'] == 'array': TOPO[r['array_id']] = r
        for r in csv.DictReader(open(f'{a.final}/cat53/{G}_existing/53b/53b_transitions.tsv'), delimiter='\t'):
            if r['unit'] == 'array': TRANS[r['array_id']][(r['from'], r['to'])] += int(r['count'])
lib = list(csv.DictReader(open(f'{a.final}/family_metrics_2026-11-07/59h_library.tsv'), delimiter='\t'))
pairs = []
if a.pairs:
    for p in a.pairs.split(','):
        x, y = p.split(':'); pairs.append((x, y))
else:
    seen = set()
    for r in sorted([r for r in lib if r['status'] == 'shared'], key=lambda r: -int(r['bp'])):
        k = tuple(sorted((r['family'], r['homolog'])))
        if k in seen: continue
        seen.add(k); pairs.append((r['family'], r['homolog']))
        if len(pairs) >= a.top: break
    for x, y in (('GarSat04', 'MacSat06'), ('GarSat09', 'MacSat03')):
        if tuple(sorted((x, y))) not in seen: pairs.append((x, y)); seen.add(tuple(sorted((x, y))))
hap = lambda f: 'Gar' if f.startswith(('GarSat', 'NewSatGar')) else 'Mac'
P('# pairs', pairs)
# ---------------- copies of the wanted families (streamed once per haplotype from the 03 tables)
want = collections.defaultdict(set)
for x, y in pairs:
    for f in (x, y): want[hap(f)].add(f)
COPIES = collections.defaultdict(list)
for G in ('Gar', 'Mac'):
    runs = collections.defaultdict(list); arrays = {}
    for r in CAT[G]:
        f = r['family_v46'] or r['family']
        if f in want[G]:
            aid = (r['array_id_v46'] or r['array_id']).split('~')[0]; runs[aid].append((int(r['start0']), int(r['end']), f))
    for p in (f'{a.out45}/fClaHyb_{G}_satellites/03_monomer_copies.tsv', f'{a.new03.format(G=G)}/03_monomer_copies.tsv'):
        if not os.path.exists(p): continue
        with open(p) as fh:
            rd = csv.reader(fh, delimiter='\t'); hdr = next(rd); ix = {h: i for i, h in enumerate(hdr)}
            for row in rd:
                aid = row[ix['array_id']]
                if aid not in runs: continue
                if row[ix['complete_or_truncated']] != 'complete' or float(row[ix['best_alignment_identity']] or 0) < 80: continue
                s, e = int(row[ix['genomic_start_0based']]), int(row[ix['genomic_end']])
                for rs, re_, f in runs[aid]:
                    if s >= rs and e <= re_: COPIES[f].append(row[ix['sequence']].upper()); break
    P(G, 'copies collected', {f: len(COPIES[f]) for f in want[G]})
# null consensus pools per haplotype (units of satellite families with a run of record)
POOL = {G: [(f, s) for (g, f), s in U.items() if g == G and (g, f) in MET and MET[(g, f)]['element_class_59b'] == 'satellite'] for G in ('Gar', 'Mac')}
def null_pool(G, L, n):
    c = [s for f, s in POOL[G] if 0.75 * L <= len(s) <= 1.25 * L] or [s for f, s in POOL[G]]; return random.sample(c, min(n, len(c)))
# composition null: random cross-haplotype unit pairs of similar length
def jsd_null(L, n=300):
    out = []
    for _ in range(n):
        s1 = random.choice([s for f, s in POOL['Gar'] if 0.75 * L <= len(s) <= 1.25 * L] or [s for f, s in POOL['Gar']]); s2 = random.choice([s for f, s in POOL['Mac'] if 0.75 * L <= len(s) <= 1.25 * L] or [s for f, s in POOL['Mac']])
        out.append(jsd(kmer_profile(s1), kmer_profile(s2)))
    return np.array(out)
# ---------------- per pair
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, textwrap
W = lambda t, n=78: textwrap.fill(t, n)
out = open(f'{a.dest}/pairs_check.tsv', 'w')
out.write('A\tB\tlibrary_identity_59h\tconsensus_identity\tconsensus_cols\tconsensus_gaps\tstrand\treciprocal\tsecond_best_B_side\tsecond_best_identity\tbest_minus_second\t'
          'A_to_A_median\tA_to_B_median\tB_to_B_median\tB_to_A_median\tnull_A_to_other_median\tnull_percentile_A_to_B\tnull_percentile_B_to_A\tcopy_gap_A\tcopy_gap_B\t'
          'GC_A\tGC_B\tjsd_4mer\tjsd_null_percentile\tunit_A\tunit_B\tarrays_A\tarrays_B\tloci_A\tloci_B\tlargest_locus_share_A\tlargest_locus_share_B\tcen_share_A\tcen_share_B\tdivergence_A\tdivergence_B\t'
          'topo_run_A\ttopo_run_B\tn_states_A\tn_states_B\tentropy_rate_A\tentropy_rate_B\tmemory_excess_A\tmemory_excess_B\tstructure_index_A\tstructure_index_B\tverdict\n')
for A, B in pairs:
    t0 = time.time(); gA, gB = hap(A), hap(B); uA, uB = U.get((gA, A)), U.get((gB, B))
    if not uA or not uB: P('skip (no unit)', A, B); continue
    h = best_orient(uA, uB)
    libid = next((r['homolog_identity'] for r in lib if r['family'] == A and r['homolog'] == B), '')
    recip = 'yes' if next((r['homolog'] for r in lib if r['family'] == B), '') == A else 'no'
    hb = sorted([(ref, idt) for ref, idt in HITS.get((gA, A), []) if ref.startswith('MacSat' if gA == 'Gar' else 'GarSat') and ref != B], key=lambda x: -x[1])
    sb, sbi = (hb[0] if hb else ('', 0.0))
    cA = random.sample(COPIES[A], min(a.n_copies, len(COPIES[A]))); cB = random.sample(COPIES[B], min(a.n_copies, len(COPIES[B])))
    AA = np.array([copy_vs_cons(c, uA) for c in cA]); AB = np.array([copy_vs_cons(c, uB) for c in cA]); BB = np.array([copy_vs_cons(c, uB) for c in cB]); BA = np.array([copy_vs_cons(c, uA) for c in cB])
    nullA = [np.median([copy_vs_cons(c, s) for c in cA[:40]]) for s in null_pool(gB, len(uA), a.n_null)]; nullB = [np.median([copy_vs_cons(c, s) for c in cB[:40]]) for s in null_pool(gA, len(uB), a.n_null)]
    pctA = 100 * np.mean(np.array(nullA) < np.median(AB)) if len(nullA) else float('nan'); pctB = 100 * np.mean(np.array(nullB) < np.median(BA)) if len(nullB) else float('nan')
    gapA = 100 * (np.median(AA) - np.median(AB)); gapB = 100 * (np.median(BB) - np.median(BA))
    pA, pB = kmer_profile(uA), kmer_profile(uB); j = jsd(pA, pB); jn = jsd_null((len(uA) + len(uB)) // 2); jpct = 100 * np.mean(jn > j)
    gc = lambda s: (s.count('G') + s.count('C')) / len(s)
    mA, mB = MET.get((gA, A), {}), MET.get((gB, B), {})
    # topology of the largest existing-source run
    def topo_run(G, f):
        rs = [r for r in CAT[G] if (r['family_v46'] or r['family']) == f and r['run_id'] in TOPO]
        return max(rs, key=lambda r: int(r['bp']))['run_id'] if rs else ''
    tA, tB = topo_run(gA, A), topo_run(gB, B); TA, TB = TOPO.get(tA, {}), TOPO.get(tB, {})
    verdict = 'same_family' if h['ident'] >= 0.80 and gapA <= 10 and gapB <= 10 and min(pctA, pctB) >= 99 and recip == 'yes' else ('sister_family' if h['ident'] >= 0.70 and min(pctA, pctB) >= 95 else 'not_supported')
    out.write('\t'.join(map(str, [A, B, libid, round(h['ident'], 3), h['cols'], h['gaps'], h['strand'], recip, sb, round(sbi, 3), round(h['ident'] - sbi, 3),
              round(np.median(AA), 3), round(np.median(AB), 3), round(np.median(BB), 3), round(np.median(BA), 3), round(float(np.median(nullA)), 3) if nullA else '', round(pctA, 1), round(pctB, 1), round(gapA, 1), round(gapB, 1),
              round(gc(uA), 3), round(gc(uB), 3), round(j, 3), round(jpct, 1), len(uA), len(uB), mA.get('n_arrays', ''), mB.get('n_arrays', ''), mA.get('n_loci', ''), mB.get('n_loci', ''), mA.get('largest_locus_share', ''), mB.get('largest_locus_share', ''),
              mA.get('bp_share_in_cen_interval', ''), mB.get('bp_share_in_cen_interval', ''), mA.get('divergence_pct', ''), mB.get('divergence_pct', ''), tA, tB, TA.get('n_states', ''), TB.get('n_states', ''), TA.get('entropy_rate_bits', ''), TB.get('entropy_rate_bits', ''),
              TA.get('memory_excess_bits', ''), TB.get('memory_excess_bits', ''), TA.get('structure_index', ''), TB.get('structure_index', ''), verdict])) + '\n'); out.flush()
    # ---------------- figure: 2 x 3
    fig, axs = plt.subplots(2, 3, figsize=(16.5, 10), constrained_layout=True); fig.patch.set_facecolor('white'); cA_, cB_ = '#8a1c2b', '#1f2d4d'
    ax = axs[0, 0]; q, c_, t = h['q'], h['c'], h['t']; n = len(q); x = np.arange(n)
    col = ['#2a9d4a' if ch == '|' else ('#c0392b' if q[i] != '-' and t[i] != '-' else '#999999') for i, ch in enumerate(c_)]
    ax.bar(x, np.ones(n), width=1, color=col, linewidth=0); ax.set_xlim(0, n); ax.set_yticks([]); ax.set_xlabel('alignment column'); ax.set_title(f"1  monomer: unit {A} ({len(uA)} bp) vs unit {B} ({len(uB)} bp)\n" + W(f"identity {h['ident']:.3f} over {h['cols']} columns ({h['gaps']} gaps), strand {h['strand']}; library call {libid}; reciprocal {recip}; 2nd-best other-haplotype family {sb or '-'} {sbi:.2f}"), fontsize=8, loc='left')
    ax.text(0.01, 0.5, 'green match · red mismatch · grey gap', transform=ax.transAxes, fontsize=7.5, va='center', color='white', bbox=dict(fc='black', alpha=0.4, lw=0))
    ax = axs[0, 1]; data = [AA, AB, BB, BA, np.array(nullA), np.array(nullB)]; labels = ['A→A', 'A→B', 'B→B', 'B→A', 'A→null', 'B→null']
    vp = ax.violinplot([d * 100 for d in data if len(d)], showmedians=True, widths=0.8)
    for body, c0 in zip(vp['bodies'], [cA_, cA_, cB_, cB_, '#999999', '#999999']): body.set_facecolor(c0); body.set_alpha(0.35)
    ax.set_xticks(range(1, 7)); ax.set_xticklabels(labels, fontsize=8); ax.set_ylabel('identity (%)'); ax.set_xlabel(f'A = {A} copies, B = {B} copies; →X = aligned to unit X; null = random units of the other haplotype', fontsize=7.5)
    ax.set_title(f'1  copy level (n = {len(cA)} / {len(cB)} confident copies)\n' + W(f"medians: A→A {np.median(AA):.2f}, A→B {np.median(AB):.2f} (gap {gapA:.1f} pts); B→B {np.median(BB):.2f}, B→A {np.median(BA):.2f} (gap {gapB:.1f}); A→B above {pctA:.0f} % and B→A above {pctB:.0f} % of the null"), fontsize=8, loc='left')
    ax = axs[0, 2]; keys = sorted(set(pA) | set(pB)); ax.scatter([pA.get(k, 0) * 100 for k in keys], [pB.get(k, 0) * 100 for k in keys], s=10, color='#555555', alpha=0.7); lim = max(max(pA.values()), max(pB.values())) * 100
    ax.plot([0, lim], [0, lim], color='#bbbbbb', lw=0.8); ax.set_xlabel(f'4-mer frequency (%) in unit {A}'); ax.set_ylabel(f'in unit {B}'); ax.set_title('2  composition of the two units\n' + W(f"GC {gc(uA):.2f} / {gc(uB):.2f}; 4-mer JSD {j:.3f}: more similar than {jpct:.0f} % of random cross-haplotype unit pairs of this length"), fontsize=8, loc='left')
    ax = axs[1, 0]
    for f, G, c0 in ((A, gA, cA_), (B, gB, cB_)):
        bp = np.array([int(r['bp']) for r in CAT[G] if (r['family_v46'] or r['family']) == f]); ax.hist(np.log10(bp), bins=np.arange(2, 6.2, 0.2), alpha=0.45, color=c0, label=f'{f}: {len(bp)} arrays, median {np.median(bp) / 1e3:.1f} kb, max {bp.max() / 1e3:.0f} kb')
    ax.set_xlabel('log10 array length (bp)'); ax.set_ylabel('arrays'); ax.legend(fontsize=7, frameon=False, loc='upper left'); ax.set_title('3  genomic organisation (A / B)\n' + W(f"loci {mA.get('n_loci', '?')} / {mB.get('n_loci', '?')}; largest-locus share {mA.get('largest_locus_share', '?')} / {mB.get('largest_locus_share', '?')}; bp in candidate centromere {mA.get('bp_share_in_cen_interval', '?')} / {mB.get('bp_share_in_cen_interval', '?')}; divergence {mA.get('divergence_pct', '?')} / {mB.get('divergence_pct', '?')} %; chromosomes {mA.get('n_chromosomes', '?')} / {mB.get('n_chromosomes', '?')}"), fontsize=8, loc='left')
    for k, (f, G, c0, tr, T) in enumerate(((A, gA, cA_, tA, TA), (B, gB, cB_, tB, TB))):
        ax = axs[1, 1 + k]
        if tr and TRANS.get(tr):
            st = sorted({s for pr in TRANS[tr] for s in pr}, key=lambda s: (not s.startswith('V'), s)); M = np.zeros((len(st), len(st)))
            for (u, v), n_ in TRANS[tr].items(): M[st.index(u), st.index(v)] = n_
            Mn = M / np.maximum(M.sum(1, keepdims=True), 1); im = ax.imshow(Mn, cmap='Reds' if k == 0 else 'Blues', vmin=0, vmax=1); ax.set_xticks(range(len(st))); ax.set_yticks(range(len(st))); ax.set_xticklabels(st, fontsize=6, rotation=90); ax.set_yticklabels(st, fontsize=6)
            ax.set_xlabel('next state'); ax.set_ylabel('state'); plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03, label='P(next | state)')
            ax.set_title(f"4  topology, largest run of {f} ({tr})\n" + W(f"{T.get('n_states', '?')} states; entropy rate {float(T.get('entropy_rate_bits', 0)):.2f} bits; memory excess {float(T.get('memory_excess_bits', 0)):.3f}; structure index {float(T.get('structure_index', 0)):.2f}; order-2 excess {T.get('order2_excess_bits', '?')}", 60), fontsize=8, loc='left')
        else: ax.text(0.5, 0.5, f'no 53b topology for {f}\n(the final chain ran 53/53b on satellite-class runs only;\ngene-array families and new-source runs have none)', ha='center', va='center', fontsize=9); ax.set_axis_off()
    fig.suptitle(f'Library pair check: A = {A} ({gA}) vs B = {B} ({gB}); verdict {verdict}', fontsize=12, x=0.01, ha='left')
    fig.savefig(f'{a.dest}/pair_{A}_{B}.png', dpi=450); fig.savefig(f'{a.dest}/pair_{A}_{B}.pdf'); plt.close(fig)
    P(f'{A} vs {B}: consensus {h["ident"]:.3f}, copies A→A {np.median(AA):.2f} A→B {np.median(AB):.2f} B→B {np.median(BB):.2f} B→A {np.median(BA):.2f}, null pct {pctA:.0f}/{pctB:.0f}, JSD pct {jpct:.0f}, recip {recip}, verdict {verdict} ({time.time() - t0:.0f}s)')
out.close(); P('# done')
