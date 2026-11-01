"""Synthetic tests for the array-definition stage: gaps.py, chains.py, 61, 61b, 62, 63 (+ plot63), 64, 65 (2026-10-26; rewritten 2026-11-01, D20).
Runs both ways:  python -P tests/test_gaps_jsd.py   (every test_* function, in order)   and   python -P -m pytest tests -q.
All data are synthetic (220-bp or 170-bp monomers, 1-2 % copy noise, random spacers) and written to a temp dir; nothing real is read.
Known truths: run edges, '-' array orientation, a 400-bp short-spacer chain (D01), a composite chain with homologous spacers (D06),
a true 4-mer HOR vs a two-class mixture (D03), i.i.d. random monomer strings (D04)."""
import os, sys, random, tempfile, subprocess, csv, importlib.util, statistics
HERE = os.path.dirname(os.path.abspath(__file__)); SC = os.path.join(HERE, '..', 'satchain'); NL = chr(10)   # Unix line ends: .fai assumes 61 bytes per line
def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(SC, name + '.py')); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
RC = str.maketrans('ACGT', 'TGCA')
def rc(s): return s.translate(RC)[::-1]
COPY_COLS = ['array_id', 'chrom', 'genomic_start_0based', 'genomic_end', 'monomer_class', 'strand', 'reverse_complemented', 'best_alignment_identity',
             'complete_or_truncated', 'sequence']
CAT_COLS = ['run_id', 'chrom', 'start0', 'end', 'bp', 'array_id', 'family', 'unit_bp', 'n_conf', 'reverse_fraction', 'G3', 'element_class', 'name_v46']


def write_genome(d, chroms):
    """chroms: list of (name, seq); writes g.fa (60 bp lines) + .fai; returns the FASTA path"""
    fa = os.path.join(d, 'g.fa')
    with open(fa, 'w', newline=NL) as fh:
        for n, s in chroms: fh.write(f'>{n}\n' + '\n'.join(s[i:i + 60] for i in range(0, len(s), 60)) + '\n')
    off = 0; fai = []
    for n, s in chroms:
        hdr = len(n) + 2; fai.append(f'{n}\t{len(s)}\t{off + hdr}\t60\t61\n'); off += hdr + len(s) + (len(s) + 59) // 60
    open(fa + '.fai', 'w', newline=NL).writelines(fai)
    return fa


def run(script, args):
    r = subprocess.run([sys.executable, '-P', os.path.join(SC, script)] + [str(x) for x in args], capture_output=True, text=True)
    assert r.returncode == 0, f'{script} failed:\n{r.stdout[-1500:]}\n{r.stderr[-3000:]}'
    return r.stdout


def tsv(p): return list(csv.DictReader(open(p), delimiter='\t'))


_D = {}
def data():
    """one synthetic world, built once: chains of 220-bp monomer runs (2 % noise) with known edges.
    chr1: 5 runs x 15 copies, 3-kb random spacers, '+'          chrR: the same array given as a '-' array (consensus rc, copies rc'd)
    chrS: 3 runs, 400-bp spacers (short spacer, D01)             chrC: 5 runs, the SAME 3-kb spacer in every interruption (composite, D06)
    chr2: 2 runs around 1000 bp + 100 N + 1000 bp (assembly gap: never chained)   chr3: one 30-copy run (null positions)"""
    if _D: return _D
    rng = random.Random(1)
    def mut(s, p): return ''.join(c if rng.random() > p else rng.choice('ACGT') for c in s)
    def rnd(n): return ''.join(rng.choice('ACGT') for _ in range(n))
    mono = rnd(220); A = 220 * 15
    def build(nrun, spacer_of):
        seq = ''; copies = []; bounds = []
        for k in range(nrun):
            if k: seq += spacer_of(k)
            s0 = len(seq)
            for i in range(15):
                c = mut(mono, 0.02); copies.append((len(seq), len(seq) + len(c), c)); seq += c
            bounds.append((s0, len(seq)))
        return seq, copies, bounds
    c1, cp1, b1 = build(5, lambda k: rnd(3000)); cS, cpS, bS = build(3, lambda k: rnd(400))
    te = rnd(3000); cC, cpC, bC = build(5, lambda k: mut(te, 0.02))
    c2 = ''.join(mut(mono, 0.02) for _ in range(15)) + rnd(1000) + 'N' * 100 + rnd(1000) + ''.join(mut(mono, 0.02) for _ in range(15))
    b2 = [(0, A), (A + 2100, 2 * A + 2100)]
    c3 = ''.join(mut(mono, 0.02) for _ in range(30))
    d = tempfile.mkdtemp(prefix='satchain_test_')
    chroms = [('chr1', c1, cp1, b1, '+'), ('chrR', c1, cp1, b1, '-'), ('chrS', cS, cpS, bS, '+'), ('chrC', cC, cpC, bC, '+')]
    fa = write_genome(d, [(n, s) for n, s, *_ in chroms] + [('chr2', c2), ('chr3', c3)])
    gb = os.path.join(d, 'gaps.bed'); open(gb, 'w', newline=NL).write(f'chr2\t{A + 1000}\t{A + 1100}\nchr3\t3300\t3300\nchr1\t0\t0\n')
    cat, cop, con = (os.path.join(d, x) for x in ('cat.tsv', 'copies.tsv', 'cons.tsv'))
    with open(cat, 'w', newline='') as fh, open(cop, 'w', newline='') as fc, open(con, 'w', newline='') as fo:
        w = csv.writer(fh, delimiter='\t'); w.writerow(CAT_COLS)
        wc = csv.writer(fc, delimiter='\t'); wc.writerow(COPY_COLS)
        wo = csv.writer(fo, delimiter='\t'); wo.writerow(['array_id', 'monomer_class', 'n_copies_in_class', 'consensus_sequence'])
        truth = {}
        for n, s, cps, bds, st in chroms:
            aid = 'A_' + n
            for i, (a, b) in enumerate(bds):
                w.writerow([f'{n}_r{i}', n, a, b, b - a, aid, 'F', 220, 15, 1.0 if st == '-' else 0.0, 'PASS', 'satellite', 'SatF'])
                if i > 0: truth[(n, f'{n}_r{i}', 'start')] = a
                if i < len(bds) - 1: truth[(n, f'{n}_r{i}', 'end')] = b
            for a, b, c in cps: wc.writerow([aid, n, a, b, 'C01', st, 'yes' if st == '-' else 'no', 98, 'complete', rc(c) if st == '-' else c])
            wo.writerow([aid, 'C01', len(cps), rc(mono) if st == '-' else mono])
        for i, (a, b) in enumerate(b2): w.writerow([f'chr2_r{i}', 'chr2', a, b, b - a, 'A_chr2', 'F', 220, 15, 0.0, 'PASS', 'satellite', 'SatF'])
    _D.update(dir=d, fa=fa, gaps=gb, cat=cat, copies=cop, cons=con, truth=truth, A=A, mono=mono)
    return _D


# ---------------------------------------------------------------- gaps.py + chains.py
def test_gaps():
    D = data(); gm = load('gaps'); A = D['A']
    G = gm.Gaps(D['gaps'])
    checks = [G.between('chr2', A, A + 2100), not G.between('chr2', 0, A), G.between('chr3', 3000, 3500), not G.between('chr3', 3300, 3500),
              not G.between('chr1', 0, 100), G.count() == 2]
    Gs = gm.Gaps(None, D['fa'], min_n=10); checks.append(Gs.between('chr2', A, A + 2100))
    print('gaps.py checks:', checks); assert all(checks)


def test_chains():
    D = data(); gm = load('gaps'); cm = load('chains')
    rows = tsv(D['cat']); G = gm.Gaps(D['gaps'])
    ch, jx, st = cm.build_chains(rows, max_gap=50000, min_runs=2, gaps=G)
    print('chains:', st)
    assert st['links_refused_at_gaps'] == 1 and all(c['array_id'] != 'A_chr2' for c in ch), 'chr2 must not be chained across the N gap'
    assert sorted(c['n_runs'] for c in ch) == [3, 5, 5, 5] and len(jx) == 2 + 4 * 3
    assert all(0 <= j['spacer_bp'] <= 3000 for j in jx) and {j['spacer_bp'] for j in jx if j['chrom'] == 'chrS'} == {400}
    # gap 0 links (D23); family_v46 preferred over family (D02); '~g' stripped (D12)
    r = [dict(run_id='a', chrom='c', start0='0', end='100', array_id='X~g1', family='F', family_v46='V'), dict(run_id='b', chrom='c', start0='100', end='200', array_id='X~g2', family='F', family_v46='V'),
         dict(run_id='c', chrom='c', start0='300', end='400', array_id='X~g2', family='F', family_v46='W')]
    ch2, jx2, st2 = cm.build_chains(r, min_runs=2)
    assert len(ch2) == 1 and ch2[0]['array_id'] == 'X' and ch2[0]['n_runs'] == 2 and st2['links_refused_family'] == 1 and jx2[0]['spacer_bp'] == 0
    ch3, _, _ = cm.build_chains(r, min_runs=4); assert ch3 == []


# ---------------------------------------------------------------- 61b
def test_61b_jsd():
    D = data(); out = os.path.join(D['dir'], 'o61b'); sw = os.path.join(D['dir'], 'sweep61b.tsv')
    print(run('61b_jsd_boundaries.py', ['--catalogue', D['cat'], '--genome', D['fa'], '--gaps', D['gaps'], '--dest', out, '--tag', 't', '--arrays', 'A_chr1,A_chrS,A_chr2',
                                        '--k', '3,4,5', '--w-units', '1,2', '--search', '300', '--null', '60', '--sweep', '--sweep-out', sw])[-700:])
    J = tsv(os.path.join(out, '61b_jsd_junctions_t.tsv')); S = tsv(os.path.join(out, '61b_jsd_summary_t.tsv'))[0]; ST = tsv(os.path.join(out, '61b_jsd_stability_t.tsv'))
    print('summary:', S)
    assert int(S['chain_links_refused_at_gaps']) == 1, 'chr2 runs must not be chained across the N gap'
    assert {x['array_id'] for x in J} == {'A_chr1', 'A_chrS'}, 'chr1 and chrS form chains, chr2 does not'
    assert all(x['family_left'] == 'F' and x['spacer_bp'] in ('3000', '400') for x in J)
    assert all(x['short_spacer'] == ('1' if x['chrom'] == 'chrS' else '0') for x in J), 'short_spacer = spacer < 2 w (400 < 440)'
    tr = D['truth']
    for chrom, lim in (('chr1', 60), ('chrS', 50)):           # w = 1 unit = 220 bp: edges within 0.3 unit (synthetic median ~20 bp)
        for wu in ('1.0', '2.0'):
            sub = [x for x in J if x['w_units'] == wu and x['chrom'] == chrom]; calls = [x['call'] for x in sub]
            errs = [abs(int(x['jsd_edge']) - tr[(chrom, x['run_id'], 'end' if x['type'] == 'OUT' else 'start')]) for x in sub if x['call'] == 'boundary']
            print(f'61b {chrom} w = {wu} unit(s): calls', {c: calls.count(c) for c in set(calls)}, '| |JSD edge - true| bp: median', statistics.median(errs), 'max', max(errs))
            if wu == '1.0': assert calls.count('boundary') == len(sub) and max(errs) <= lim, f'{chrom} one-unit windows: every boundary called within {lim} bp'
            assert all(x['ramp'].startswith('up' if x['type'] == 'OUT' else 'down') for x in sub), 'ramp direction fixed by junction type (D19)'
    assert all(x['primary_edge'] != '' and int(x['primary_settings_boundary']) == 3 for x in ST), 'primary edge (w = 1 unit) from all 3 k'
    assert os.path.exists(sw) and len(tsv(sw)) == 6, 'sweep table: 3 k x 2 w'


# ---------------------------------------------------------------- 61
def test_61_edges():
    D = data(); out = os.path.join(D['dir'], 'o61'); sw = os.path.join(D['dir'], 'sweep61.tsv')
    print(run('61_unit_boundaries.py', ['--catalogue', D['cat'], '--copies', D['copies'], '--consensus', D['cons'], '--genome', D['fa'], '--dest', out, '--tag', 't',
                                        '--gaps', D['gaps'], '--flank', '100,300', '--sweep', '--sweep-out', sw])[-600:])
    J = tsv(os.path.join(out, '61_junctions_t.tsv')); tr = D['truth']
    assert len(J) == 2 * 2 * (4 + 4 + 2 + 4) and all(x['refined'] for x in J), '2 flanks x 2 edges per junction, all refined'
    assert tsv(os.path.join(out, '61_skipped_chains_t.tsv')) == []
    for chrom in ('chr1', 'chrR', 'chrS', 'chrC'):
        sub = [x for x in J if x['chrom'] == chrom]
        errs = [int(x['refined']) - tr[(chrom, x['run_id'], 'end' if (x['type'] == 'OUT') == (chrom != 'chrR') else 'start')] for x in sub]
        print(f'61 {chrom}: refined - true (bp) {sorted(set(errs))}; types {sorted(set(x["type"] for x in sub))}')
        assert max(abs(e) for e in errs) <= 5, f'{chrom}: edges recovered within 5 bp (incl. the "-" array)'
    # homology null (D06): random spacers (chr1, chrR) are (almost) never homologous, the composite chain (chrC) always is; medoids excluded
    for chrom, lim in (('chr1', 0.15), ('chrR', 0.15), ('chrC', 1.0)):
        h = [int(x['homologous_null']) for x in J if x['chrom'] == chrom and x['homologous_null'] != '']; p = sorted({x['null_p99'] for x in J if x['chrom'] == chrom})
        print(f'61 {chrom}: homologous_null share {sum(h) / len(h):.2f} of {len(h)}, null_n {sorted({x["null_n"] for x in J if x["chrom"] == chrom})}, null_p99 {p}, '
              f'flank_identity median {statistics.median(float(x["flank_identity"]) for x in J if x["chrom"] == chrom and x["is_medoid"] == "0"):.2f}')
        assert len(h) == 2 * 2 * 3 and (sum(h) / len(h) <= lim if lim < 1 else sum(h) == len(h)), f'{chrom}: homologous_null share vs {lim}'
    U = [l.rstrip('\n').split('\t') for l in open(os.path.join(out, '61_units_t.bed'))]
    print('61 units:', [(u[0], u[1], u[2], u[5], u[6], u[7]) for u in U][:4], '...', len(U))
    assert len(U) == 3 + 3 + 1 + 3 and all(len(u) == 8 for u in U), 'units = IN -> next IN: n_runs - 2 per chain'
    assert all(abs(int(u[6]) - 3300) <= 10 for u in U) and all(abs(int(u[7]) - (400 if u[0] == 'chrS' else 3000)) <= 10 for u in U), 'satellite_bp / interruption_bp (D11)'
    assert {u[5] for u in U if u[0] == 'chrR'} == {'-'} and {u[5] for u in U if u[0] != 'chrR'} == {'+'}
    assert os.path.exists(sw) and len(tsv(sw)) > 0


# ---------------------------------------------------------------- 62 (+ 63, 65 share the HOR genome)
_H = {}
def hor_data():
    """chrH: a true 4-mer HOR (4 monomers 10 % apart, 120 copies, 1 % noise); chrM: a two-class mixture (2 classes 8 % apart, i.i.d., 120 copies);
    chrN: a monomeric run (120 copies, 1 % noise)."""
    if _H: return _H
    rng = random.Random(7)
    def mut(s, p): return ''.join(c if rng.random() > p else rng.choice('ACGT') for c in s)
    base = ''.join(rng.choice('ACGT') for _ in range(170))
    mons = [mut(base, 0.10) for _ in range(4)]; A_, B_ = mut(base, 0.08), mut(base, 0.08)
    sets = {'chrH': [mut(mons[i % 4], 0.01) for i in range(120)], 'chrM': [mut(rng.choice([A_, B_]), 0.01) for _ in range(120)], 'chrN': [mut(base, 0.01) for _ in range(120)]}
    d = tempfile.mkdtemp(prefix='satchain_hor_'); fa = write_genome(d, [(n, ''.join(v)) for n, v in sets.items()])
    cat, cop = os.path.join(d, 'cat.tsv'), os.path.join(d, 'copies.tsv')
    with open(cat, 'w', newline='') as fh, open(cop, 'w', newline='') as fc:
        w = csv.writer(fh, delimiter='\t'); w.writerow(CAT_COLS); wc = csv.writer(fc, delimiter='\t'); wc.writerow(COPY_COLS)
        for n, v in sets.items():
            L = sum(len(x) for x in v); w.writerow([f'{n}_r0', n, 0, L, L, 'A_' + n, 'F', 170, len(v), 0.5 if n == 'chrM' else 0.0, 'FAIL' if n == 'chrM' else 'PASS', 'satellite', 'SatH'])
            p = 0
            for i, x in enumerate(v): wc.writerow(['A_' + n, n, p, p + len(x), f'C{i % 4 + 1:02d}' if n == 'chrH' else 'C01', '+', 'no', 97, 'complete', x]); p += len(x)
    _H.update(dir=d, fa=fa, cat=cat, copies=cop, n=120)
    return _H


def test_62_hor_vs_mixture():
    H = hor_data(); out = os.path.join(H['dir'], 'o62'); sw = os.path.join(H['dir'], 'sweep62.tsv')
    print(run('62_hor_detect.py', ['--catalogue', H['cat'], '--copies', H['copies'], '--genome', H['fa'], '--dest', out, '--tag', 't', '--perm', '500', '--ks', '6,15',
                                   '--sweep', '--sweep-out', sw])[-900:])
    C = {r['run_id']: r for r in tsv(os.path.join(out, '62_hor_calls_t.tsv'))}; L = tsv(os.path.join(out, '62_lag_identity_t.tsv'))
    for k, r in C.items(): print(f'62 {k}: copies {r["copies"]} max_lag {r["max_lag"]} best_lag {r["best_lag"]} excess {r["best_excess"]} p_max {r["p_max"]} q {r["q_max"]} '
                                 f'candidate {r["hor_candidate"]!r} HOR {r["hor_monomers"]!r} flag {r["strand_flag"]!r} k15 interval {r["k15_top_interval"]}')
    assert C['chrH_r0']['hor_monomers'] == '4' and C['chrH_r0']['hor_bp_estimate'] == str(4 * 170), 'true 4-mer HOR must be called with n = 4'
    assert C['chrM_r0']['hor_monomers'] == '' and C['chrN_r0']['hor_monomers'] == '', 'two-class mixture and monomeric run: no HOR'
    assert C['chrM_r0']['strand_flag'] == 'mixed_strand' and C['chrM_r0']['G3'] == 'FAIL' and C['chrH_r0']['strand_flag'] == '', 'reverse_fraction 0.5 flagged (D07)'
    assert all(int(r['max_lag']) == 40 for r in C.values()) and {int(x['lag']) for x in L} == set(range(1, 41)) and all(int(x['pairs']) == min(150, 120 - int(x['lag'])) for x in L)
    assert float(C['chrM_r0']['q_max']) > 0.05 and float(C['chrH_r0']['q_max']) <= 0.05
    SW = tsv(sw); print('62 sweep rows', len(SW), 'HOR calls per setting', sorted({int(r['hor_calls_fdr']) for r in SW}))
    assert len(SW) == 3 * 2 * 2 * 3 and all(int(r['hor_calls_fdr']) == 1 for r in SW), 'every sweep setting calls exactly the one true HOR'


def test_63_tracks_and_plot():
    H = hor_data(); out = os.path.join(H['dir'], 'o63'); L = 120 * 170
    r = subprocess.run([sys.executable, '-P', os.path.join(SC, '63_expansion_tracks.py'), '--chrom', 'chrH', '--start', '0', '--end', str(L + 200), '--name', 'h', '--genome', H['fa'],
                        '--dest', out, '--mode', 'hor', '--copies', H['copies']], capture_output=True, text=True)
    assert r.returncode != 0 and '--hor' in r.stderr + r.stdout, '--mode hor without --hor must stop (D10)'
    print(run('63_expansion_tracks.py', ['--chrom', 'chrH', '--start', 0, '--end', L + 200, '--name', 'h', '--genome', H['fa'], '--dest', out, '--mode', 'hor', '--hor', 4,
                                         '--copies', H['copies'], '--bins', 20, '--sweep', '--sweep-out', os.path.join(H['dir'], 'sweep63.tsv')])[-400:])
    U = tsv(os.path.join(out, '63_units_h.tsv')); S = tsv(os.path.join(out, '63_summary_h.tsv'))[0]
    assert len(U) == 30 and all(int(u['length']) >= 4 * 160 for u in U), '120 adjacent copies -> 30 HOR blocks'
    assert [int(u['win_units']) for u in U][-3:] == [3, 2, 1] and all(int(u['exact_win']) <= int(u['win_units']) for u in U), 'truncated last windows marked'
    print('63 summary:', S); assert float(S['exact_rate_null']) <= float(S['exact_rate']) + 1e-9 and int(S['truncated_windows']) == 9
    png = os.path.join(out, 'h.png'); print(run('plot63_array.py', [out, png, 'h:test']))
    assert os.path.exists(png) and os.path.exists(os.path.join(out, 'h.pdf')), 'plot63 writes PNG + PDF (D26)'


def test_65_null_and_hor():
    H = hor_data(); base = ['--catalogue', H['cat'], '--copies', H['copies'], '--t-min', '0.90', '--t-max', '0.95', '--t-step', '0.01', '--null', '10']
    out = os.path.join(H['dir'], 'o65')
    print(run('65_hor_mining.py', base + ['--runs', 'chrM_r0', '--dest', out, '--name', 'mix'])[-500:])
    T = tsv(os.path.join(out, 'mix_thresholds.tsv')); hs = tsv(os.path.join(out, 'mix_hors.tsv'))
    print('65 mixture: coverage', [r['hor_coverage'] for r in T], 'null', [r['null_coverage'] for r in T], 'excess', [r['excess'] for r in T])
    assert all(abs(float(r['excess'])) < 0.08 for r in T), 'i.i.d. two-class string: coverage within 0.08 of the shuffled null (D04)'
    assert hs == [] and all(r['chosen'] == '0' for r in T), 'no threshold above its null -> no HOR reported'
    print(run('65_hor_mining.py', base + ['--runs', 'chrH_r0', '--dest', out, '--name', 'hor', '--sweep', '--sweep-out', os.path.join(H['dir'], 'sweep65.tsv')])[-500:])
    T = tsv(os.path.join(out, 'hor_thresholds.tsv')); hs = tsv(os.path.join(out, 'hor_hors.tsv')); ch = [r for r in T if r['chosen'] == '1']
    print('65 HOR: chosen', ch[0] if ch else None, '| top', hs[0] if hs else None)
    assert ch and hs and hs[0]['length'] == '4' and float(hs[0]['coverage']) > 0.9 and float(ch[0]['null_coverage']) < 0.2, '4-mer HOR found, far above its null'
    SW = tsv(os.path.join(H['dir'], 'sweep65.tsv')); assert {r['min_copies'] for r in SW} == {'2', '3', '4'}


# ---------------------------------------------------------------- 64 (smoke: block permutation scheme, gaps, BH)
def test_64_blocks():
    D = data(); out = os.path.join(D['dir'], 'o64')
    print(run('64_object_profiles.py', ['--catalogue', D['cat'], '--copies', D['copies'], '--genome', D['fa'], '--dest', out, '--gaps', D['gaps'], '--group', 'chrom,array_id',
                                        '--loo-block', 'run_id,run_id', '--ks', '3,5', '--perm', '99', '--sweep', '--sweep-out', os.path.join(D['dir'], 'sweep64.tsv')])[-700:])
    O = tsv(os.path.join(out, '64_objects.tsv')); S = tsv(os.path.join(out, '64_separation.tsv'))
    types = {r['object_type']: 0 for r in O}
    for r in O: types[r['object_type']] += 1
    print('64 objects:', types); assert types['run'] == 20 and types['gap'] == 14 and types['copy'] == 18 * 15, 'chr2 gap object skipped (assembly gap, D17)'
    assert {r['strand'] for r in O if r['chrom'] == 'chrR'} == {'-'} and {r['strand'] for r in O if r['chrom'] == 'chr1'} == {'+'}, 'run / gap objects oriented by majority strand'
    for r in S: print('64', r['grouping'], r['object_type'], r['k'], 'groups', r['groups'], 'blocks', r['blocks'], r['perm_scheme'], 'p', r['p_perm'], 'q', r['q_perm'], 'loo', r['loo_accuracy'])
    cp = [r for r in S if r['object_type'] == 'copy' and r['grouping'] == 'chrom']
    assert cp and all(r['perm_scheme'] == 'between_blocks' and r['groups'] == '4' and int(r['blocks']) == 18 and r['q_perm'] != '' for r in cp)
    assert all(float(r['p_perm']) > 0.05 for r in cp), 'chr1 / chrR / chrS / chrC copies are the same monomer: no separation by chromosome'


if __name__ == '__main__':
    for name, fn in [(k, v) for k, v in list(globals().items()) if k.startswith('test_') and callable(v)]:
        print(f'=== {name}'); fn()
    print('ALL TESTS PASSED')
