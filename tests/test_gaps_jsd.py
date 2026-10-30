"""Synthetic tests for gaps.py and 61b_jsd_boundaries.py (2026-10-26). usage: python tests/test_gaps_jsd.py  (writes into a temp dir)"""
import os, sys, random, tempfile, subprocess, csv, importlib.util
HERE = os.path.dirname(os.path.abspath(__file__)); SC = os.path.join(HERE, '..', 'satchain')
spec = importlib.util.spec_from_file_location('gaps', os.path.join(SC, 'gaps.py')); gm = importlib.util.module_from_spec(spec); spec.loader.exec_module(gm)
rng = random.Random(1); NL = chr(10)          # Unix line ends: the .fai offsets assume 61 bytes per line
def mut(s, p): return ''.join(c if rng.random() > p else rng.choice('ACGT') for c in s)
mono = ''.join(rng.choice('ACGT') for _ in range(220))
sat = lambda n: ''.join(mut(mono, 0.02) for _ in range(n))
rnd = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
# chr1: sat(15) | spacer 3000 | sat(15) | spacer 3000 | sat(15)      boundaries known
# chr2: sat(15) | 1000 random + 100 N + 1000 random | sat(15)        -> the two runs must NOT be chained (gap between them)
A, SP = 220 * 15, 3000
c1 = sat(15) + rnd(SP) + sat(15) + rnd(SP) + sat(15)
b1 = [(0, A), (A + SP, 2 * A + SP), (2 * A + 2 * SP, 3 * A + 2 * SP)]
c2 = sat(15) + rnd(1000) + 'N' * 100 + rnd(1000) + sat(15)
b2 = [(0, A), (A + 2100, 2 * A + 2100)]
d = tempfile.mkdtemp(); fa = os.path.join(d, 'g.fa')
with open(fa, 'w', newline=NL) as fh:
    for n, s in (('chr1', c1), ('chr2', c2), ('chr3', sat(30))):
        fh.write(f'>{n}\n' + '\n'.join(s[i:i + 60] for i in range(0, len(s), 60)) + '\n')
off = 0; fai = []
for n, s in (('chr1', c1), ('chr2', c2), ('chr3', sat(30))):
    hdr = len(n) + 2; fai.append(f'{n}\t{len(s)}\t{off + hdr}\t60\t61\n'); off += hdr + len(s) + (len(s) + 59) // 60
open(fa + '.fai', 'w', newline=NL).writelines(fai)
gb = os.path.join(d, 'gaps.bed'); open(gb, 'w', newline=NL).write(f'chr2\t{A + 1000}\t{A + 1100}\nchr3\t3300\t3300\nchr1\t0\t0\n')
# --- gaps.py
G = gm.Gaps(gb)
checks = [G.between('chr2', A, A + 2100), not G.between('chr2', 0, A), G.between('chr3', 3000, 3500), not G.between('chr3', 3300, 3500),
          not G.between('chr1', 0, 100), G.count() == 2]
Gs = gm.Gaps(None, fa, min_n=10); checks.append(Gs.between('chr2', A, A + 2100))
print('gaps.py checks:', checks); assert all(checks)
# --- 61b
cat = os.path.join(d, 'cat.tsv')
with open(cat, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t'); w.writerow(['run_id', 'chrom', 'start0', 'end', 'array_id', 'family', 'unit_bp'])
    for i, (s, e) in enumerate(b1): w.writerow([f'r1_{i}', 'chr1', s, e, 'A1', 'F', 220])
    for i, (s, e) in enumerate(b2): w.writerow([f'r2_{i}', 'chr2', s, e, 'A2', 'F', 220])
out = os.path.join(d, 'o')
r = subprocess.run([sys.executable, os.path.join(SC, '61b_jsd_boundaries.py'), '--catalogue', cat, '--genome', fa, '--gaps', gb, '--dest', out, '--tag', 't',
                    '--k', '3,4,5', '--w', '100,200', '--search', '300', '--null', 60], capture_output=True, text=True) if False else \
    subprocess.run([sys.executable, os.path.join(SC, '61b_jsd_boundaries.py'), '--catalogue', cat, '--genome', fa, '--gaps', gb, '--dest', out, '--tag', 't',
                    '--k', '3,4,5', '--w-units', '1,2', '--search', '300', '--null', '60'], capture_output=True, text=True)
print(r.stdout[-600:], r.stderr[-1500:])
J = list(csv.DictReader(open(os.path.join(out, '61b_jsd_junctions_t.tsv')), delimiter='\t'))
S = list(csv.DictReader(open(os.path.join(out, '61b_jsd_summary_t.tsv')), delimiter='\t'))[0]
print('summary:', S)
assert int(S['chain_links_refused_at_gaps']) == 1, 'chr2 runs must not be chained across the N gap'
assert all(x['array_id'] == 'A1' for x in J), 'only chr1 forms a chain'
true = {('r1_0', 'OUT'): A, ('r1_1', 'IN'): A + SP, ('r1_1', 'OUT'): 2 * A + SP, ('r1_2', 'IN'): 2 * A + 2 * SP}
for wu in ('1.0', '2.0'):
    sub = [x for x in J if x['w_units'] == wu]; calls = [x['call'] for x in sub]
    errs = [abs(int(x['jsd_edge']) - true[(x['run_id'], x['type'])]) for x in sub if x['call'] == 'boundary']
    print(f'w = {wu} unit(s): calls', {c: calls.count(c) for c in set(calls)}, '| |JSD edge - true| bp: median', sorted(errs)[len(errs) // 2], 'max', max(errs))
    if wu == '1.0': assert calls.count('boundary') == len(sub) and max(errs) <= 50, 'one-unit windows: every boundary called within 50 bp'
print('ALL TESTS PASSED')