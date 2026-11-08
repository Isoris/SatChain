#!/usr/bin/env python3
"""regenerate_numbers.py (2026-11-01): every number quoted in README.md / docs/METHODS.md / docs/STEP_METHOD.md, computed from the shipped
tables of one SatChain run folder (audit HB01, HB04, HB19, D22, C05, C20: README numbers that no table reproduced). Prints a markdown
block and writes benchmarks/NUMBERS_<date>.tsv (key, Gar, Mac, source table). Run it after every rerun and paste; never type a number by hand.
usage: regenerate_numbers.py <run folder with cat59/ hmm56_<G>/ lin55/ benchmarks/ cat53/> [--tags Gar,Mac]"""
import sys, os, csv, collections, argparse, datetime
ap = argparse.ArgumentParser(); ap.add_argument('V'); ap.add_argument('--tags', default='Gar,Mac'); a = ap.parse_args(); V = a.V; TAGS = a.tags.split(',')
N = collections.OrderedDict()
def put(k, g, v, src): N.setdefault(k, {})[g] = v; N[k]['source'] = src
def mb(x): return f'{x / 1e6:.2f}'
def rows(p): return list(csv.DictReader(open(p), delimiter='\t')) if os.path.exists(p) else []
for G in TAGS:
    cat = rows(f'{V}/cat59/59_catalogue_{G}.tsv')
    if cat:
        sat = [r for r in cat if r['element_class'] == 'satellite']; ga = [r for r in cat if r['element_class'].startswith('gene_array')]
        put('runs (all)', G, len(cat), '59_catalogue'); put('satellite families Mb', G, mb(sum(int(r['bp']) for r in sat)), '59_catalogue element_class=satellite')
        put('satellite Mb, unit >= 100 bp', G, mb(sum(int(r['bp']) for r in sat if r['unit_class'] == 'satellite')), '59_catalogue unit_class=satellite')
        put('minisatellite Mb (7-99 bp)', G, mb(sum(int(r['bp']) for r in sat if r['unit_class'] == 'mini')), '59_catalogue')
        put('gene arrays Mb / runs', G, f"{mb(sum(int(r['bp']) for r in ga))} / {len(ga)}", '59_catalogue element_class=gene_array')
        put('gene-array families', G, ', '.join(sorted({r['family'] + ':' + r['element_class'].split(':')[1] for r in ga})), '59_catalogue')
        new = [r for r in cat if r['source'] == 'new']; put('new runs / Mb', G, f"{len(new)} / {mb(sum(int(r['bp']) for r in new))}", '59_catalogue source=new')
        put('new runs carried from v4.6 draft', G, sum(1 for r in cat if r.get('carried_from', 'none') != 'none'), '59_catalogue carried_from')
        put('runs passing G1+G2+G3 (%)', G, f"{100 * sum(r['gates_passed'] == '3' for r in cat) / len(cat):.0f}", '59_catalogue gates_passed (G2 >= 85 since 2026-11-01)')
        put('G2 fails (identity < 85)', G, sum(r['G2'] == 'FAIL' for r in cat), '59_catalogue')
        put('overlapping runs resolved', G, sum(r['overlap_note'] != 'none' for r in cat), '59_catalogue overlap_note')
        cen = sum(int(r['bp']) for r in sat if r['putative_centromeric_array_L1'] == 'yes'); ror = sum(int(r['bp']) for r in sat if r['cen_candidate_ror'] == 'yes'); tot = sum(int(r['bp']) for r in sat)
        put('satellite share in step-31 centromere intervals (%, all satellite families)', G, f'{100 * cen / tot:.0f}', '59_catalogue putative_centromeric_array_L1'); put('satellite share in largest run cluster per chromosome (%, all satellite families)', G, f'{100 * ror / tot:.0f}', '59_catalogue cen_candidate_ror')
        s100 = [r for r in sat if r['unit_class'] == 'satellite']; t100 = sum(int(r['bp']) for r in s100)
        put('satellite (unit >= 100 bp) share in step-31 centromere intervals / largest run cluster (%)', G, f"{100 * sum(int(r['bp']) for r in s100 if r['putative_centromeric_array_L1'] == 'yes') / t100:.0f} / {100 * sum(int(r['bp']) for r in s100 if r['cen_candidate_ror'] == 'yes') / t100:.0f}", '59_catalogue unit_class=satellite')
        fams = collections.Counter(r['family_v46'] if 'family_v46' in r else r['family'] for r in sat); put('satellite families', G, len(fams), '59_catalogue family_v46')
        cl = collections.Counter(r['seq_cluster'] for r in sat); put('runs by seq_cluster', G, '; '.join(f'{k}:{v}' for k, v in sorted(cl.items())), '59_catalogue seq_cluster')
        clf = collections.defaultdict(set)
        for r in sat: clf[r['seq_cluster']].add(r.get('family_v46', r['family']))
        put('families with a seq_cluster (not none/unassigned)', G, sum(len(v) for k, v in clf.items() if k not in ('none', 'unassigned', 'na')), '59_catalogue (catalogue families, C05)')
        tiers = collections.Counter(r['tier'] for r in sat); put('runs by tier', G, '; '.join(f'{k}:{v}' for k, v in sorted(tiers.items())), '59_catalogue tier (59d, family_v46)')
    for p, k in ((f'{V}/hmm56_{G}/56_summary.tsv', 'HMM Mb at posterior 0.5 / 0.8 / 0.95'),):
        s = rows(p)
        if s: put(k, G, ' / '.join(r['hmm_Mb'] for r in s), '56_summary'); put('HMM runs at posterior 0.95', G, s[-1]['n_hmm_runs'], '56_summary')
    hp = f'{V}/hmm56_{G}/56_hmm_params.tsv'
    if os.path.exists(hp):
        d = dict(l.rstrip('\n').split('\t')[:2] for l in open(hp) if '\t' in l); put('HMM P(A->A) / P(B->B)', G, f"{d.get('P(A->A)')} / {d.get('P(B->B)')}", '56_hmm_params')
    b = rows(f'{V}/hmm56_{G}/val56_benchmark_{G}.tsv')
    if b:
        put('val56 runs kept at q0.8 / 0.9 / 0.95', G, ' / '.join(f"{r['runs_kept']}/{r['runs_total']}" for r in b), 'val56_benchmark')
        put('val56 Mb added beyond v4.5 at q0.8 / 0.9 / 0.95', G, ' / '.join(r['added_beyond_36_Mb'] for r in b), 'val56_benchmark')
        put('val56 fraction of v4.5 bp kept at q0.8 / 0.9 / 0.95', G, ' / '.join(r['frac_36_bp_kept'] for r in b), 'val56_benchmark')
        put('val56 background zero share', G, b[0]['background_zero_share'], 'val56_benchmark (bins outside all candidate intervals)')
    gs = rows(f'{V}/cat59/59b_gene_arrays_{G}_sweep.tsv')
    if gs: put('gene-array Mb at ratio cut 0.5 / 0.7 / 0.8 (unit share 0.5)', G, ' / '.join(r.get('gene_array', '0') for r in gs if r['unit_cut'] == '0.5'), '59b_gene_arrays_sweep')
    gw = rows(f'{V}/cat59/59_gates_sweep_{G}.tsv')
    if gw: put('runs passing all gates, sweep min-max (%)', G, f"{100 * min(int(r['n_pass_all']) for r in gw) / int(gw[0]['n_runs']):.0f}-{100 * max(int(r['n_pass_all']) for r in gw) / int(gw[0]['n_runs']):.0f}", '59_gates_sweep')
    for rule, folder in (('copy', 'cat59'), ('min', 'cat59_minrule')):
        te = rows(f'{V}/{folder}/59c_te_{G}.tsv')
        if te:
            tb = sum(int(r['bp']) for r in te); ins = [sum(int(r[f'TE_insert_bp_{c}']) for r in te) for c in ('0.2', '0.3', '0.5')]
            put(f'TE insertions inside runs, intrinsic cut 0.2 / 0.3 / 0.5 (carry rule {rule}) Mb', G, ' / '.join(mb(x) for x in ins), f'59c_te ({folder})')
            put(f'TE insertions as % of run bp (carry rule {rule})', G, f'{100 * ins[0] / tb:.1f}-{100 * ins[2] / tb:.1f}', f'59c_te ({folder})')
    it = rows(f'{V}/cat59/59e_interruptions_{G}.tsv')
    if it:
        c = collections.Counter(); bpc = collections.Counter()
        for r in it: c[r['class']] += 1; bpc[r['class']] += int(r['bp'])
        put('interruptions (n / Mb) by class', G, '; '.join(f'{k}:{c[k]}/{mb(bpc[k])}' for k in sorted(c)), '59e_interruptions'); put('interruptions across array_ids', G, sum(r['array_relation'] == 'across_arrays' for r in it), '59e_interruptions')
        ucl = {r['run_id']: r['unit_class'] for r in cat}; c2 = collections.Counter(); b2 = collections.Counter()      # the same, satellite-class (unit >= 100 bp) runs only
        for r in it:
            if ucl.get(r['run_left']) == 'satellite': c2[r['class']] += 1; b2[r['class']] += int(r['bp'])
        put('interruptions between unit >= 100 bp runs (n / Mb) by class', G, '; '.join(f'{k}:{c2[k]}/{mb(b2[k])}' for k in sorted(c2)), '59e_interruptions x 59_catalogue unit_class')
    sw = rows(f'{V}/cat59/59e_sweep_{G}.tsv')
    if sw:
        te_mb = [float(r['Mb_class']) for r in sw if r['class'] == 'TE_insertion' and r['max_gap'] == '50000']; put('TE-filled interruptions Mb, sweep min-max (gap <= 50 kb)', G, f'{min(te_mb):.2f}-{max(te_mb):.2f}', '59e_sweep')
    bt = rows(f'{V}/benchmarks/threshold_band_test_{G}_summary.tsv')
    if bt:
        put('band test: tandem share above null p95 by lag-identity band', G, '; '.join(f"{r['band']}:{float(r['tandem_above_null_p95']):.2f}" for r in bt), 'threshold_band_test_summary')
    s53 = rows(f'{V}/cat53/{G}_existing/53_array_summary.tsv')
    if s53:
        put('53 runs: unsegmentable / copy-windowed / bp-windowed', G, f"{sum(r.get('segmentable') == 'no' for r in s53)} / {sum(r.get('window_mode', '').startswith('copy') for r in s53)} / {sum(r.get('window_mode') == 'bp' and r.get('segmentable') != 'no' for r in s53)}", '53_array_summary (bp-windowed excludes unsegmentable runs, fixed 2026-11-09)')
        put('53 runs with > 1 segment', G, sum(int(r.get('n_segments', 1)) > 1 for r in s53 if r.get('n_segments', '1').isdigit()), '53_array_summary')
kc = rows(f'{V}/lin55/55_K_choice.tsv')
if kc:
    ch = [r for r in kc if r.get('chosen') == 'yes']
    if ch: put('seq clusters: chosen K (rule / fallback), min ARI across k, min ARI restarts, silhouette', 'pooled', f"{ch[0]['K']} ({ch[0]['note']}), {ch[0]['min_ARI_across_k']}, {ch[0]['min_ARI_restarts']}, {ch[0]['mean_silhouette']}", '55_K_choice')
cs = rows(f'{V}/lin55/55_cluster_summary.tsv')
if cs and kc:
    K = ch[0]['K']; put('seq clusters at chosen K (k = 6): n, Gar, Mac, median len, median 5S frac', 'pooled', '; '.join(f"{r['cluster']}: {r['n']}, {r['Gar']}, {r['Mac']}, {r['median_len']}, {r['median_frac_5S']}" for r in cs if r['k'] == '6' and r['K'] == K), '55_cluster_summary')
hs = rows(f'{V}/lin55/55_haplotype_stability.tsv')
if hs and kc:
    put('seq clusters: ARI pooled vs own-haplotype at chosen K (k = 6)', 'pooled', '; '.join(f"{r['haplotype']}:{r['ARI_pooled_vs_own']}" for r in hs if r['k'] == '6' and r['K'] == K), '55_haplotype_stability')
fc = rows(f'{V}/lin55/55d_family_cluster.tsv')
if fc:
    for G in TAGS: put('families assigned a seq_cluster / unassigned (sampled families)', G, f"{sum(r['subgenome'] == G and r['seq_cluster'] != 'unassigned' for r in fc)} / {sum(r['subgenome'] == G and r['seq_cluster'] == 'unassigned' for r in fc)}", '55d_family_cluster')
date = datetime.date.today().isoformat(); out = os.path.join(os.path.dirname(os.path.abspath(__file__)), f'NUMBERS_{date}.tsv')
with open(out, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n'); w.writerow(['key'] + TAGS + ['pooled', 'source'])
    for k, d in N.items(): w.writerow([k] + [d.get(g, '') for g in TAGS] + [d.get('pooled', ''), d['source']])
print(f'| number | {" | ".join(TAGS)} | pooled | source |\n|---|{"---|" * (len(TAGS) + 2)}')
for k, d in N.items(): print(f'| {k} | {" | ".join(str(d.get(g, "")) for g in TAGS)} | {d.get("pooled", "")} | {d["source"]} |')
print('written', out)
