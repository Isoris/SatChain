#!/usr/bin/env python3
"""59g_unit_content.py  (SatChain v4.7, 2026-11-07; what is INSIDE the repeat unit of every catalogue family)

Why: 59b labels a family 'gene_array' from annotation COORDINATES only (overlap ratio against one swept cut-off). It never looks at the
monomer sequence, so a family at ratio 0.71 vs cut 0.7 (MacSat06) flips with the cut, a 5S-derived satellite (GarSat04) and a 5S gene
array (GarSat49) are told apart only by whether barrnap happened to call the gene, and composite units (gene + spacer, satellite + TE
fragment, two satellite monomers) are invisible.

What it does, per haplotype G and per catalogue family (family_v46):
  1. unit = the step-03 consensus of the dominant class of the family's largest run (existing out/ table or the new-candidate table);
     the unit is doubled so a circular alignment can start anywhere.
  2. references extracted from the genome by the annotation GFFs, one class each: 5S (barrnap), rRNA_other (5.8S/18S/28S), tRNA
     (tRNAscan), U1..U6 (Rfam), miR-430, miRNA_other (curated), Histone3, UnaL2 (Rfam; a LINE-3'-end family), rfam_other;
     plus the EDTA TE library (TElib.fa) and every family consensus of BOTH haplotypes (04_satellite_family_consensuses.fa).
  3. k-mer screen (k = 13 for genes, 15 for TEs and satellites; >= 3 / 4 shared k-mers, both strands) -> candidates -> local
     alignment (parasail Smith-Waterman, match 2 / mismatch -3 / gap 5,2) -> best hit per class: identity, reference coverage,
     fraction of the unit covered, position on the unit.
  4. segment accounting of the unit: gene / TE / other-satellite / unassigned (spacer), greedy by alignment score, no overlaps.
     A TE consensus that covers >= 85 % of the unit at >= 85 % identity and is <= 3 unit lengths long is flagged
     'telib_is_this_satellite' (EDTA built a 'TE' from the satellite) and NOT counted as TE content.
  5. composition of the unit (length, GC, CpG o/e, entropies k = 1..3) and the catalogue context (bp, runs, mean copy identity,
     strand coherence, EDTA overlap, 59b ratio and unit share, centromere flags, tier).
  6. THE DECISION (one per family; 2026-11-07 review): unit_content_rule in {gene_array, gene_derived_satellite, satellite, unresolved}
     from three axes kept as separate columns:
       sequence_axis      intact_gene (identity >= 0.90, reference coverage >= 0.80) / gene_derived (>= 0.70, >= 0.5) / none
       annotation (59b)   element_class_59b, ratio_59b (threshold-sensitive when |ratio - 0.7| < 0.1)
       unit_organization  gene_plus_spacer (unit <= 2.5 x gene) / gene_in_long_unit / no_gene
     gene_array = intact gene AND 59b agrees AND not threshold-sensitive; gene_derived_satellite = gene-derived unit, no annotation call;
     satellite = no gene homology, no annotation call; unresolved = the axes disagree or the ratio is threshold-sensitive.
     TE content, other-family content and composite flags are reported (te_frac, othersat_frac, composite_sat_pair, segments) but do
     not enter the decision. The sweep (gene identity 0.80 / 0.90 x coverage 0.6 / 0.8) is kept for the record.
  7. EXPLORATORY ONLY (no decision depends on it): pooled PCA of the feature table, k-means / Ward with silhouette, cluster x label
     tables and the PCA scatter. Purpose: show whether the labelled groups are separable in the data; nothing downstream reads it.
Outputs (<dest>/): 59g_unit_content_<G>.tsv, 59g_refs_<G>.fa, 59g_hits_<G>.tsv (every aligned hit kept), 59g_rule_sweep.tsv,
59g_pca_scores.tsv, 59g_pca_loadings.tsv, 59g_clusters.tsv, 59g_cluster_vs_rule.tsv, 59g_cluster_vs_59b.tsv, 59g_pca.png/.pdf, 59g_log.txt.
usage: 59g_unit_content.py <final cat59 dir> <annotation root> <v4.5 out root> <new03 dir pattern with {G}> <dest> [--tags Gar,Mac]
Evidence level L1 (sequence); the labels are descriptions of unit content, not function.
"""
import sys, os, csv, gzip, re, math, collections, time, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from faidx_lite import Fai
import parasail

ap = argparse.ArgumentParser()
ap.add_argument('cat59'); ap.add_argument('ann'); ap.add_argument('out45'); ap.add_argument('new03'); ap.add_argument('dest')
ap.add_argument('--tags', default='Gar,Mac'); ap.add_argument('--kgene', type=int, default=13); ap.add_argument('--krep', type=int, default=15)
ap.add_argument('--min-kmers-gene', type=int, default=3); ap.add_argument('--min-kmers-rep', type=int, default=4)
ap.add_argument('--max-cand', type=int, default=6, help='candidates per class sent to alignment'); ap.add_argument('--seed', type=int, default=7)
a = ap.parse_args(); os.makedirs(a.dest, exist_ok=True)
LOG = open(os.path.join(a.dest, '59g_log.txt'), 'a')
def P(*x):
    s = ' '.join(str(v) for v in x); print(s, flush=True); LOG.write(s + '\n'); LOG.flush()
P('# 59g start', time.strftime('%Y-%m-%d %H:%M'), 'args', vars(a))
COMP = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
rc = lambda s: s.translate(COMP)[::-1]
MAT = parasail.matrix_create('ACGT', 2, -3)
ENC = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

def kmers(seq, k):
    """2-bit integer k-mers of a sequence (skips k-mers with non-ACGT)."""
    out = []; v = 0; n = 0; mask = (1 << (2 * k)) - 1
    for ch in seq.upper():
        c = ENC.get(ch)
        if c is None: n = 0; v = 0; continue
        v = ((v << 2) | c) & mask; n += 1
        if n >= k: out.append(v)
    return out

def sw(query, ref):
    """Local alignment query (doubled unit) vs ref. Returns dict with score, identity, aligned columns, query span, ref span."""
    r = parasail.sw_trace_striped_16(query, ref, 5, 2, MAT)
    tb = r.traceback; q, c, t = tb.query, tb.comp, tb.ref
    cols = len(q); matches = sum(1 for x in c if x == '|'); ql = len(q.replace('-', '')); tl = len(t.replace('-', ''))
    qe = r.end_query + 1; te = r.end_ref + 1
    return dict(score=int(r.score), ident=matches / cols if cols else 0.0, cols=cols, q0=qe - ql, q1=qe, t0=te - tl, t1=te, ref_len=len(ref))

def read_gff(path, want):
    """Yield (chrom, start0, end, strand, class) for GFF lines whose class, given by want(line fields), is not None."""
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        for l in f:
            if l[:1] == '#' or not l.strip(): continue
            c = l.rstrip('\n').split('\t')
            if len(c) < 9: continue
            cl = want(c)
            if cl: yield c[0], int(c[3]) - 1, int(c[4]), c[6], cl

def attr(s): return dict(x.split('=', 1) for x in s.split(';') if '=' in x)

def gene_refs(G, fai):
    """Reference sequences per gene class from the host's own annotation (sequence taken from the genome, strand applied)."""
    A = os.path.join(a.ann, f'fClaHyb_{G}')
    feats = []
    feats += list(read_gff(f'{A}/fClaHyb_{G}_rRNA_barrnap/fClaHyb_{G}.rRNA.sorted.gff3.gz', lambda c: '5S' if '5S_rRNA' in c[8] else 'rRNA_other'))
    feats += list(read_gff(f'{A}/fClaHyb_{G}_tRNA_tRNAscan/fClaHyb_{G}.tRNA.sorted.gff3.gz', lambda c: 'tRNA' if c[2] in ('tRNA', 'gene', 'pseudogene') else None))
    def rf(c):
        n = attr(c[8]).get('Name', '')
        if n in ('U1', 'U2', 'U3', 'U4', 'U5', 'U6'): return 'snRNA_' + n
        if n == '5S_rRNA': return '5S'
        if n == 'tRNA': return 'tRNA'
        if n == 'miR-430' or n.lower().startswith('mir-430'): return 'miR-430'
        if n in ('Histone3', 'UnaL2'): return n
        if n.startswith(('LSU', 'SSU', '5_8S')): return 'rRNA_other'
        return 'rfam_other'
    feats += list(read_gff(f'{A}/fClaHyb_{G}_ncRNA_rfam/fClaHyb_{G}.rfam.all_classes.sorted.gff3.gz', rf))
    feats += list(read_gff(f'{A}/fClaHyb_{G}_miRNA/fClaHyb_{G}.miRNA.curated.sorted.gff3.gz', lambda c: 'miR-430' if 'mir-430' in attr(c[8]).get('Name', '').lower() else 'miRNA_other'))
    # one sequence per distinct (class, chrom, start): dedupe tracks that report the same gene twice
    seen = set(); refs = []; per = collections.Counter(); cap = {'tRNA': 3000, 'rRNA_other': 300, 'rfam_other': 3000, 'UnaL2': 2000}
    cache = {}   # each chromosome read once (fai.get per feature re-read the whole chromosome: hours instead of seconds)
    for chrom, s, e, strand, cl in sorted(feats):
        if (cl, chrom, s) in seen or e - s < 40 or e - s > 6000: continue
        seen.add((cl, chrom, s))
        if per[cl] >= cap.get(cl, 5000): continue
        try:
            if chrom not in cache: cache[chrom] = fai.get(chrom)
            seq = cache[chrom][s:e].upper()
        except Exception: continue
        if strand == '-': seq = rc(seq)
        per[cl] += 1; refs.append((f'{cl}|{chrom}:{s}-{e}', cl, seq))
    P(G, 'gene references per class', dict(per))
    return refs

def read_fa(path):
    name = None; buf = []
    with open(path) as f:
        for l in f:
            if l.startswith('>'):
                if name: yield name, ''.join(buf)
                name = l[1:].split()[0]; buf = []
            else: buf.append(l.strip())
    if name: yield name, ''.join(buf)

def consensus_tables(G):
    """(array_id, class) -> consensus sequence, from the v4.5 out table and the new-candidate table."""
    d = {}
    for p in (os.path.join(a.out45, f'fClaHyb_{G}_satellites', '03_array_consensus.tsv'), os.path.join(a.new03.format(G=G), '03_array_consensus.tsv')):
        if not os.path.exists(p): P('missing consensus table', p); continue
        for r in csv.DictReader(open(p), delimiter='\t'):
            d[(r['array_id'], r['monomer_class'])] = r['consensus_sequence'].upper()
    return d

def entropy(seq, k):
    c = collections.Counter(seq[i:i + k] for i in range(len(seq) - k + 1)); n = sum(c.values())
    return -sum(v / n * math.log2(v / n) for v in c.values()) / (2 * k) if n else 0.0

# ------------------------------------------------------------------ main, per haplotype
tags = a.tags.split(','); allrows = []
satfa = {}  # family consensus of both haplotypes (other-satellite references)
for G in tags:
    p = os.path.join(a.out45, f'fClaHyb_{G}_satellites', '04_satellite_family_consensuses.fa')
    if os.path.exists(p):
        for n, s in read_fa(p): satfa[n] = s.upper()
P('family consensus references (both haplotypes):', len(satfa))
sweep_rows = []
for G in tags:
    t0 = time.time()
    A = os.path.join(a.ann, f'fClaHyb_{G}'); fai = Fai(f'{A}/fClaHyb_{G}_genome/fClaHyb_{G}_LG.fa')
    cat = list(csv.DictReader(open(os.path.join(a.cat59, f'59_catalogue_{G}.tsv')), delimiter='\t'))
    b59 = {r['family']: r for r in csv.DictReader(open(os.path.join(a.cat59, f'59b_gene_arrays_{G}.tsv')), delimiter='\t')} if os.path.exists(os.path.join(a.cat59, f'59b_gene_arrays_{G}.tsv')) else {}
    cons = consensus_tables(G)
    fam = collections.defaultdict(list)
    for r in cat: fam[r['family_v46'] or r['family']].append(r)
    # ---- units
    units = {}; nounit = 0
    for f, rs in fam.items():
        big = max(rs, key=lambda r: int(r['bp']))
        seq = cons.get((big['array_id_v46'] or big['array_id'], big['dominant_class'])) or cons.get((big['array_id'], big['dominant_class']))
        if not seq or len(seq) < 7: nounit += 1; continue
        units[f] = dict(unit=seq, run=big)
    P(G, 'families', len(fam), 'with a unit consensus', len(units), 'without', nounit)
    # ---- references
    grefs = gene_refs(G, fai)
    terefs = [(n, 'TE', s.upper()) for n, s in read_fa(f'{A}/fClaHyb_{G}_repeats_EDTA/fClaHyb_{G}.TElib.fa')]
    srefs = [(n, 'SAT', s) for n, s in satfa.items()]
    with open(os.path.join(a.dest, f'59g_refs_{G}.fa'), 'w') as fh:
        for n, cl, s in grefs: fh.write(f'>{n}\n{s}\n')
    # ---- k-mer index of the doubled units (forward strand; references are screened on both strands)
    idx = {a.kgene: collections.defaultdict(list), a.krep: collections.defaultdict(list)}
    fams = sorted(units)
    for fi, f in enumerate(fams):
        dd = units[f]['unit'] * 2
        for k in idx:
            for km in set(kmers(dd, k)): idx[k][km].append(fi)
    def screen(refs, k, mink):
        """shared k-mer counts per (family, reference), both strands of the reference."""
        hits = collections.defaultdict(collections.Counter)
        for ri, (n, cl, s) in enumerate(refs):
            for strand, seq in (('+', s), ('-', rc(s))):
                cnt = collections.Counter()
                for km in set(kmers(seq, k)):
                    for fi in idx[k].get(km, ()): cnt[fi] += 1
                for fi, c in cnt.items():
                    if c >= mink and c > hits[fi][(ri, strand)]: hits[fi][(ri, strand)] = c
        return hits
    hg = screen(grefs, a.kgene, a.min_kmers_gene); ht = screen(terefs, a.krep, a.min_kmers_rep); hs = screen(srefs, a.krep, a.min_kmers_rep)
    P(G, 'k-mer screen done', f'{time.time() - t0:.0f}s', 'families with gene cand', len(hg), 'TE cand', len(ht), 'sat cand', len(hs))
    # ---- alignments
    hitw = open(os.path.join(a.dest, f'59g_hits_{G}.tsv'), 'w'); hitw.write('family\tref_group\tref_class\tref\tstrand\tshared_kmers\tscore\tidentity\tcols\tunit_start\tunit_end\tunit_frac\tref_cov\tref_len\n')
    rows = []
    for fi, f in enumerate(fams):
        u = units[f]['unit']; L = len(u); dd = u * 2; big = units[f]['run']
        best = {}   # class -> best hit dict
        segs = []   # (score, start, end, kind, label) on the unit (0..L)
        def consider(group, cl, name, strand, kc, s):
            h = sw(dd, s if strand == '+' else rc(s))
            if h['cols'] < 20: return
            ufrac = min(1.0, (h['q1'] - h['q0']) / L); rcov = (h['t1'] - h['t0']) / h['ref_len']
            hitw.write(f"{f}\t{group}\t{cl}\t{name}\t{strand}\t{kc}\t{h['score']}\t{h['ident']:.3f}\t{h['cols']}\t{h['q0'] % L}\t{h['q0'] % L + (h['q1'] - h['q0'])}\t{ufrac:.3f}\t{rcov:.3f}\t{h['ref_len']}\n")
            rec = dict(group=group, cls=cl, name=name, strand=strand, kc=kc, ident=h['ident'], ufrac=ufrac, rcov=rcov, score=h['score'], q0=h['q0'] % L, span=h['q1'] - h['q0'], ref_len=h['ref_len'])
            key = cl if group == 'gene' else group
            if key not in best or rec['score'] > best[key]['score']: best[key] = rec
            return rec
        # genes: top candidates per class
        bycl = collections.defaultdict(list)
        for (ri, strand), kc in hg.get(fi, {}).items(): bycl[grefs[ri][1]].append((kc, ri, strand))
        for cl, lst in bycl.items():
            for kc, ri, strand in sorted(lst, reverse=True)[:a.max_cand]:
                r_ = consider('gene', cl, grefs[ri][0], strand, kc, grefs[ri][2])
                if r_: segs.append((r_['score'], r_['q0'], r_['q0'] + r_['span'], 'gene', cl))
        # TEs
        te_self = ''
        for kc, (ri, strand) in sorted(((kc, k) for k, kc in ht.get(fi, {}).items()), reverse=True)[:a.max_cand]:
            r_ = consider('TE', 'TE', terefs[ri][0], strand, kc, terefs[ri][2])
            if not r_: continue
            if r_['ufrac'] >= 0.85 and r_['ident'] >= 0.85 and r_['ref_len'] <= 3 * L:
                te_self = te_self or terefs[ri][0]; continue       # EDTA consensus built from this satellite: not TE content
            segs.append((r_['score'], r_['q0'], r_['q0'] + r_['span'], 'TE', terefs[ri][0]))
        # other satellite families (never the family itself)
        for kc, (ri, strand) in sorted(((kc, k) for k, kc in hs.get(fi, {}).items()), reverse=True)[:a.max_cand]:
            if srefs[ri][0] == f: continue
            r_ = consider('SAT', 'SAT', srefs[ri][0], strand, kc, srefs[ri][2])
            if r_ and r_['ident'] >= 0.80: segs.append((r_['score'], r_['q0'], r_['q0'] + r_['span'], 'SAT', srefs[ri][0]))
        # greedy non-overlapping accounting of the unit (circular: positions mod L)
        cover = np.zeros(L, dtype=np.int8); kinds = {'gene': 1, 'TE': 2, 'SAT': 3}; placed = []
        for sc, s0, s1, kind, lab in sorted(segs, reverse=True):
            pos = np.arange(s0, min(s1, s0 + L)) % L
            if cover[pos].any(): continue
            cover[pos] = kinds[kind]; placed.append((kind, lab, len(pos)))
        gene_bp = int((cover == 1).sum()); te_bp = int((cover == 2).sum()); sat_bp = int((cover == 3).sum()); spacer_bp = L - gene_bp - te_bp - sat_bp
        # composite_sat: two different other families on disjoint parts, < 0.70 identical to each other
        satfams = [(lab, n) for kind, lab, n in placed if kind == 'SAT' and n >= 40]
        comp_sat = 0
        if len(satfams) >= 2:
            s1_, s2_ = satfams[0][0], satfams[1][0]
            if s1_ in satfa and s2_ in satfa:
                hh = sw(satfa[s1_] * 2, satfa[s2_]); comp_sat = int(hh['ident'] < 0.70 or hh['cols'] < 0.5 * min(len(satfa[s1_]), len(satfa[s2_])))
        gb = max((v for k, v in best.items() if k not in ('TE', 'SAT')), key=lambda v: v['score'], default=None)
        te = best.get('TE'); sa = best.get('SAT')
        # catalogue context
        rs = fam[f]; bp = sum(int(r['bp']) for r in rs)
        wmed = lambda key: float(np.median([float(r[key]) for r in rs if r[key] not in ('', 'NA')] or [0]))
        ident_mean = float(np.average([float(r['mean_identity_conf']) for r in rs], weights=[int(r['bp']) for r in rs]))
        strand_coh = float(np.average([abs(float(r['reverse_fraction']) - 0.5) * 2 for r in rs], weights=[int(r['bp']) for r in rs]))
        edta = float(np.average([float(r['EDTA_overlap_fraction']) for r in rs], weights=[int(r['bp']) for r in rs]))
        b = b59.get(f, {})
        row = dict(haplotype=G, family=f, name_v46=big.get('name_v46', ''), unit_len=L, unit_GC=round((u.count('G') + u.count('C')) / L, 3),
                   unit_CpG_oe=round((u.count('CG') * L) / max(1, u.count('C') * u.count('G')), 3), ent1=round(entropy(u, 1), 3), ent2=round(entropy(u, 2), 3), ent3=round(entropy(u, 3), 3),
                   n_runs=len(rs), family_bp=bp, unit_bp_catalogue=wmed('unit_bp'), copy_identity=round(ident_mean, 1), strand_coherence=round(strand_coh, 3), edta_overlap=round(edta, 3),
                   element_class_59b=rs[0]['element_class'], gene_flag_59b=rs[0]['gene_flag'], ratio_59b=b.get('observed_over_expected', ''), unit_share_59b=b.get('unit_share', ''), top_class_59b=b.get('top_class', ''),
                   tier=rs[0]['tier'], seq_cluster=rs[0]['seq_cluster'], cen_L1_runs=sum(1 for r in rs if r['putative_centromeric_array_L1'] == 'yes'), cen_ror_runs=sum(1 for r in rs if r['cen_candidate_ror'] == 'yes'),
                   gene_best_class=gb['cls'] if gb else '', gene_best_ident=round(gb['ident'], 3) if gb else 0.0, gene_best_refcov=round(gb['rcov'], 3) if gb else 0.0, gene_best_unitfrac=round(gb['ufrac'], 3) if gb else 0.0, gene_best_ref=gb['name'] if gb else '', gene_best_reflen=gb['ref_len'] if gb else 0,
                   **{f'ident_{c}': round(best[c]['ident'], 3) if c in best else 0.0 for c in ('5S', 'tRNA', 'snRNA_U1', 'snRNA_U2', 'snRNA_U4', 'snRNA_U5', 'snRNA_U6', 'miR-430', 'miRNA_other', 'Histone3', 'UnaL2', 'rRNA_other', 'rfam_other')},
                   te_best=te['name'] if te else '', te_best_ident=round(te['ident'], 3) if te else 0.0, te_best_unitfrac=round(te['ufrac'], 3) if te else 0.0, te_best_len=te['ref_len'] if te else 0, telib_is_this_satellite=te_self,
                   sat_best=sa['name'] if sa else '', sat_best_ident=round(sa['ident'], 3) if sa else 0.0, sat_best_unitfrac=round(sa['ufrac'], 3) if sa else 0.0,
                   gene_bp=gene_bp, te_bp=te_bp, othersat_bp=sat_bp, spacer_bp=spacer_bp, gene_frac=round(gene_bp / L, 3), te_frac=round(te_bp / L, 3), othersat_frac=round(sat_bp / L, 3), spacer_frac=round(spacer_bp / L, 3),
                   composite_sat_pair=';'.join(x[0] for x in satfams[:2]) if comp_sat else '', segments=';'.join(f'{k}:{lab}:{n}' for k, lab, n in placed), unit_seq=u)
        rows.append(row)
    hitw.close()
    # ---- rule label with sweep
    # ---- the decision (2026-11-07, after review): ONE decision, outcomes gene_array / gene_derived_satellite / satellite / unresolved.
    # Axis 1 (sequence identity): what the unit resembles = best gene identity x reference coverage (this step).
    # Axis 2 (annotation): 59b observed/expected ratio against its cut (0.7) and whether it sits within 0.1 of the cut.
    # Axis 3 (organization): unit_organization = gene_plus_spacer (unit <= 2.5 x gene) / gene_in_long_unit / no_gene; TE and
    # other-family content are reported as columns (te_frac, othersat_frac, composite_sat_pair) and do NOT enter the decision.
    # unresolved = the axes disagree (annotation says gene array but the unit holds no intact gene, or the unit holds an intact gene
    # the annotation ratio does not support) or the 59b ratio is threshold-sensitive (|ratio - 0.7| < 0.1).
    def seq_axis(r, idc, covc):
        gi, gc = r['gene_best_ident'], r['gene_best_refcov']
        if gi >= idc and gc >= covc: return 'intact_gene'
        if gi >= 0.70 and gc >= 0.5: return 'gene_derived'
        return 'none'
    def organization(r):
        if r['gene_best_class'] == '' or seq_axis(r, 0.70, 0.5) == 'none': return 'no_gene'
        return 'gene_plus_spacer' if r['unit_len'] <= 2.5 * max(1, r['gene_best_reflen']) else 'gene_in_long_unit'
    def rule(r, idc, covc):
        s = seq_axis(r, idc, covc); b59 = r['element_class_59b'].startswith('gene_array')
        ratio = num(r['ratio_59b'], -1.0); sensitive = ratio >= 0 and abs(ratio - 0.7) < 0.1
        if s == 'intact_gene' and b59 and not sensitive: return 'gene_array'
        if s == 'gene_derived' and not b59: return 'gene_derived_satellite'
        if s == 'none' and not b59: return 'satellite'
        return 'unresolved'
    def num(v, d=0.0):
        try: return float(v)
        except Exception: return d
    for r in rows:
        r['sequence_axis'] = seq_axis(r, 0.90, 0.80); r['unit_organization'] = organization(r)
        r['unit_content_rule'] = rule(r, 0.90, 0.80)
        for idc in (0.80, 0.90):
            for covc in (0.6, 0.8): r[f'rule_id{idc}_cov{covc}'] = rule(r, idc, covc)
    for idc in (0.80, 0.90):
        for covc in (0.6, 0.8):
            c = collections.Counter(r[f'rule_id{idc}_cov{covc}'] for r in rows); bpc = collections.Counter()
            for r in rows: bpc[r[f'rule_id{idc}_cov{covc}']] += r['family_bp']
            for lab in sorted(c): sweep_rows.append((G, idc, covc, lab, c[lab], round(bpc[lab] / 1e6, 3)))
    cols = [k for k in rows[0] if k != 'unit_seq'] + ['unit_seq'] if rows else []
    with open(os.path.join(a.dest, f'59g_unit_content_{G}.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
    P(G, 'done', len(rows), 'families', f'{time.time() - t0:.0f}s', 'rule labels', dict(collections.Counter(r['unit_content_rule'] for r in rows)),
      '| 59b vs rule disagreements:', sum(1 for r in rows if (r['element_class_59b'].startswith('gene_array')) != (r['unit_content_rule'] == 'gene_array')))
    allrows += rows
with open(os.path.join(a.dest, '59g_rule_sweep.tsv'), 'w') as fh:
    fh.write('haplotype\tgene_ident_cut\tgene_cov_cut\tlabel\tfamilies\tMb\n')
    for r in sweep_rows: fh.write('\t'.join(map(str, r)) + '\n')

# ------------------------------------------------------------------ unsupervised view, pooled
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score
FEAT = ['log_unit_len', 'unit_GC', 'unit_CpG_oe', 'ent1', 'ent2', 'ent3', 'copy_identity', 'strand_coherence', 'edta_overlap', 'log_family_bp', 'log_n_runs',
        'gene_best_ident', 'gene_best_refcov', 'gene_best_unitfrac', 'gene_frac', 'te_best_ident', 'te_frac', 'sat_best_ident', 'othersat_frac', 'spacer_frac', 'ratio_59b_num', 'unit_share_59b_num',
        'ident_5S', 'ident_tRNA', 'ident_snRNA_max', 'ident_miR-430', 'ident_UnaL2']
def num(v, d=0.0):
    try: return float(v)
    except Exception: return d
for r in allrows:
    r['log_unit_len'] = math.log10(r['unit_len']); r['log_family_bp'] = math.log10(max(1, r['family_bp'])); r['log_n_runs'] = math.log10(max(1, r['n_runs']))
    r['ratio_59b_num'] = min(3.0, num(r['ratio_59b'])); r['unit_share_59b_num'] = num(r['unit_share_59b']); r['ident_snRNA_max'] = max(r[f'ident_snRNA_U{i}'] for i in (1, 2, 4, 5, 6))
X = np.array([[num(r[c]) for c in FEAT] for r in allrows]); Xs = StandardScaler().fit_transform(X)
pca = PCA(n_components=min(8, Xs.shape[1])).fit(Xs); S = pca.transform(Xs)
P('PCA explained variance', [round(v, 3) for v in pca.explained_variance_ratio_])
with open(os.path.join(a.dest, '59g_pca_loadings.tsv'), 'w') as fh:
    fh.write('feature\t' + '\t'.join(f'PC{i + 1}' for i in range(S.shape[1])) + '\n')
    for j, c in enumerate(FEAT): fh.write(c + '\t' + '\t'.join(f'{pca.components_[i, j]:.3f}' for i in range(S.shape[1])) + '\n')
Z = S[:, :5]; best = None; ks = []
for K in range(2, 11):
    km = KMeans(n_clusters=K, n_init=25, random_state=a.seed).fit(Z); sil = silhouette_score(Z, km.labels_)
    ag = AgglomerativeClustering(n_clusters=K, linkage='ward').fit(Z); sila = silhouette_score(Z, ag.labels_)
    ks.append((K, sil, sila)); P('K', K, 'silhouette kmeans', round(sil, 3), 'ward', round(sila, 3))
    if best is None or sil > best[1]: best = (K, sil, km.labels_, ag.labels_)
K, sil, lab, laba = best; P('chosen K (max k-means silhouette on PC1-5)', K, round(sil, 3))
with open(os.path.join(a.dest, '59g_clusters.tsv'), 'w') as fh:
    fh.write('K\tsilhouette_kmeans\tsilhouette_ward\n')
    for K_, s1, s2 in ks: fh.write(f'{K_}\t{s1:.3f}\t{s2:.3f}\n')
with open(os.path.join(a.dest, '59g_pca_scores.tsv'), 'w') as fh:
    fh.write('haplotype\tfamily\tunit_content_rule\telement_class_59b\tgene_flag_59b\tkmeans_cluster\tward_cluster\t' + '\t'.join(f'PC{i + 1}' for i in range(S.shape[1])) + '\n')
    for r, s, c1, c2 in zip(allrows, S, lab, laba): fh.write(f"{r['haplotype']}\t{r['family']}\t{r['unit_content_rule']}\t{r['element_class_59b']}\t{r['gene_flag_59b']}\t{c1}\t{c2}\t" + '\t'.join(f'{v:.3f}' for v in s) + '\n')
def cross(key, fn):
    t = collections.Counter(); bpc = collections.Counter()
    for r, c in zip(allrows, lab): t[(c, r[key])] += 1; bpc[(c, r[key])] += r['family_bp']
    labs = sorted({r[key] for r in allrows})
    with open(os.path.join(a.dest, fn), 'w') as fh:
        fh.write('cluster\tn_families\tMb\t' + '\t'.join(labs) + '\tdominant\n')
        for c in range(K):
            n = sum(t[(c, l)] for l in labs); mb = sum(bpc[(c, l)] for l in labs) / 1e6
            dom = max(labs, key=lambda l: t[(c, l)])
            fh.write(f'{c}\t{n}\t{mb:.3f}\t' + '\t'.join(str(t[(c, l)]) for l in labs) + f'\t{dom}\n')
cross('unit_content_rule', '59g_cluster_vs_rule.tsv'); cross('element_class_59b', '59g_cluster_vs_59b.tsv')
# ---- figure: PC1 x PC2 coloured three ways
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, 3, figsize=(14, 4.6)); fig.patch.set_facecolor('white')
def scat(ax, labels, title, order=None):
    labs = order or sorted(set(labels)); cmap = plt.get_cmap('tab10')
    for i, l in enumerate(labs):
        m = np.array([x == l for x in labels]); ax.scatter(S[m, 0], S[m, 1], s=9, alpha=0.6, color=cmap(i % 10), label=f'{l} ({m.sum()})', linewidths=0)
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.0%})'); ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.0%})'); ax.set_title(title, fontsize=10, loc='left'); ax.legend(fontsize=6.5, frameon=False, markerscale=2)
scat(axs[0], [r['unit_content_rule'] for r in allrows], 'a  unit content, rule label (59g)')
scat(axs[1], [r['element_class_59b'].split(':')[0] + (':' + r['gene_flag_59b'].split(':')[0] if r['gene_flag_59b'] not in ('', 'none') else '') for r in allrows], 'b  annotation-overlap label (59b)')
scat(axs[2], [f'cluster {c}' for c in lab], f'c  k-means on PC1-5, K = {K} (silhouette {sil:.2f})')
fig.suptitle(f'Repeat-unit content of {len(allrows)} catalogue families (Gar + Mac), SatChain 59g, {time.strftime("%Y-%m-%d")}', fontsize=10, x=0.01, ha='left')
fig.tight_layout(); fig.savefig(os.path.join(a.dest, '59g_pca.png'), dpi=450); fig.savefig(os.path.join(a.dest, '59g_pca.pdf')); P('figure written'); P('# 59g end', time.strftime('%H:%M'))
