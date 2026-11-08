"""Narrow HOR re-test (decision 2026-11-09): steps 62 (decides) and 65 (cross-check) on the ten HOR-like centromeric arrays of thesis
Table 3.9 (v4.5 self-similarity evidence) plus three 'plateau' arrays as real-data negative controls. Purpose: confirm or correct the thesis
claim that MacSat1272 on Mac LG20 is the one higher-order repeat. No change to 62 / 65 code; their defaults, with --sweep.
Array = runs of the dominant family (family_v46) of the v4.7 catalogue overlapping the candidate centromere interval of that chromosome
(final/runs_of_record/<G>/centromere_check.tsv). Copies = v4.7 step-03 tables (existing + new runs), pre-filtered to the selected
intervals so each call reads a small file (same columns; selection by coordinates only).
Output: satellite_FISH/satchain_v4.7_2026-11-01/final/hor_narrow_2026-11-09/"""
import os, sys, csv, subprocess, time
D = r'F:\data'; V = os.path.join(D, 'satellite_FISH', 'satchain_v4.7_2026-11-01'); FIN = os.path.join(V, 'final')
SC = os.path.join(D, 'SatChain', 'satchain'); OUT = os.path.join(FIN, 'hor_narrow_2026-11-09'); os.makedirs(OUT, exist_ok=True)
csv.field_size_limit(1 << 30)
ARRAYS = [  # (genome, chrom, family, role, thesis margin evidence)
    ('Mac', 'C_mac_LG20', 'MacSat1272', 'HOR-like', '0.850 (0.05 yes)'), ('Mac', 'C_mac_LG23', 'MacSat03', 'HOR-like', '0.092 (0.05 yes)'),
    ('Mac', 'C_mac_LG27', 'MacSat03', 'HOR-like', '0.045'), ('Gar', 'C_gar_LG10', 'GarSat09', 'HOR-like', '0.045'),
    ('Gar', 'C_gar_LG25', 'GarSat03', 'HOR-like', '0.042'), ('Mac', 'C_mac_LG09', 'MacSat03', 'HOR-like', '0.040'),
    ('Mac', 'C_mac_LG05', 'MacSat03', 'HOR-like', '0.032'), ('Mac', 'C_mac_LG02', 'MacSat03', 'HOR-like', '0.025'),
    ('Mac', 'C_mac_LG24', 'MacSat03', 'HOR-like', '0.024'), ('Mac', 'C_mac_LG04', 'MacSat03', 'HOR-like', '0.022'),
    ('Gar', 'C_gar_LG11', 'GarSat03', 'control_plateau', ''), ('Gar', 'C_gar_LG09', 'GarSat04', 'control_plateau', ''),
    ('Mac', 'C_mac_LG13', 'MacSat03', 'control_plateau', ''),
]
def log(m): print(time.strftime('%H:%M:%S ') + m, flush=True); open(os.path.join(OUT, 'driver.log'), 'a').write(time.strftime('%F %T ') + m + '\n')
def tsv(p): return list(csv.DictReader(open(p, encoding='utf-8'), delimiter='\t'))
def run(args, logname):
    t0 = time.time(); r = subprocess.run([sys.executable, '-P'] + args, capture_output=True, text=True)
    open(os.path.join(OUT, logname), 'w', encoding='utf-8').write(r.stdout + '\n--- stderr ---\n' + r.stderr)
    if r.returncode: sys.exit(f'FAILED {logname}: {r.stderr[-2000:]}')
    log(f'{logname} done in {time.time() - t0:.0f} s')

sel = []   # (genome, array_label, run_ids)
for G in ('Gar', 'Mac'):
    cen = {r['chrom']: (int(r['cen_start0']), int(r['cen_end'])) for r in tsv(os.path.join(FIN, 'runs_of_record', G, 'centromere_check.tsv'))}
    cat = tsv(os.path.join(FIN, 'cat59', f'59_catalogue_{G}.tsv'))
    for g, chrom, fam, role, ev in ARRAYS:
        if g != G: continue
        s, e = cen[chrom]
        ids = [r['run_id'] for r in cat if r['chrom'] == chrom and r['family_v46'] == fam and int(r['end']) > s and int(r['start0']) < e]
        if not ids: sys.exit(f'no runs for {chrom} {fam}')
        sel.append((G, f'{chrom}_{fam}', ids, role, ev, s, e))
        log(f'{G} {chrom} {fam} {role}: {len(ids)} runs in {s}-{e}')
with open(os.path.join(OUT, 'arrays.tsv'), 'w') as f:
    f.write('genome\tarray\trole\tthesis_inter_minus_intra\tcen_start0\tcen_end\tn_runs\trun_ids\n')
    for G, lab, ids, role, ev, s, e in sel: f.write(f'{G}\t{lab}\t{role}\t{ev}\t{s}\t{e}\t{len(ids)}\t{",".join(ids)}\n')

for G in ('Gar', 'Mac'):
    wins = [(chrom, s, e) for g, lab, ids, role, ev, s, e in sel if g == G for chrom in [lab.rsplit('_', 1)[0]]]
    small = os.path.join(OUT, f'copies_{G}.tsv')
    if not os.path.exists(small):
        srcs = [os.path.join(V, 'out', f'fClaHyb_{G}_satellites', '03_monomer_copies.tsv'), os.path.join(V, f'new03v47_{G}', '03_monomer_copies.tsv')]
        n = 0; header = None
        with open(small + '.part', 'w', encoding='utf-8', newline='') as fo:
            for src in srcs:
                with open(src, encoding='utf-8', newline='') as fi:
                    h = fi.readline()
                    if header is None: header = h; fo.write(h)
                    elif h != header: sys.exit(f'copy table columns differ: {src}')
                    cols = h.rstrip('\r\n').split('\t'); ic, i0, i1 = cols.index('chrom'), cols.index('genomic_start_0based'), cols.index('genomic_end')
                    for line in fi:
                        x = line.split('\t', max(ic, i0, i1) + 1)
                        c, a, b = x[ic], int(x[i0]), int(x[i1])
                        if any(c == wc and b > ws and a < we for wc, ws, we in wins): fo.write(line); n += 1
        os.replace(small + '.part', small); log(f'{G}: {n} copies in the selected intervals')
    fa = os.path.join(D, 'fClaHyb_func_annotation_v1.0', f'fClaHyb_{G}', f'fClaHyb_{G}_genome', f'fClaHyb_{G}_LG.fa')
    cat = os.path.join(FIN, 'cat59', f'59_catalogue_{G}.tsv')
    ids = ','.join(i for g, lab, rid, *_ in sel if g == G for i in rid)
    run([os.path.join(SC, '62_hor_detect.py'), '--catalogue', cat, '--copies', small, '--genome', fa, '--dest', os.path.join(OUT, 'o62'), '--tag', G,
         '--runs', ids, '--sweep', '--sweep-out', os.path.join(OUT, 'o62', f'62_sweep_{G}.tsv')], f'62_{G}.log')
    for g, lab, rid, role, *_ in sel:
        if g != G: continue
        run([os.path.join(SC, '65_hor_mining.py'), '--catalogue', cat, '--copies', small, '--runs', ','.join(rid), '--dest', os.path.join(OUT, 'o65'),
             '--name', lab, '--sweep', '--sweep-out', os.path.join(OUT, 'o65', f'65_sweep_{lab}.tsv')], f'65_{lab}.log')
log('ALL DONE')
