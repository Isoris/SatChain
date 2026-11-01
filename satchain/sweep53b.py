#!/usr/bin/env python3
"""One-at-a-time sweep of the 53b cutoffs on one step-53 run of the RUNS OF RECORD (--runs): major-edge probability / count,
intrinsic-annotation fraction. Per setting vs default: edges whose core/accessory role changes, runs whose intrinsic annotations change,
segments whose variant structure index changes sign, library types at the first split-half null cut.
2026-11-01 (CA27 / HB17 / HB02): each setting writes to <53 run dir>/53b_sweep_<setting> (53b --dest), so the production <run>/53b is
never overwritten; 53b's stderr is not captured (shown live) and a failing setting stops the sweep with 53b's own traceback; the
library column is the current type_at_nullq*_D* name; --runs / --gaps / --allow-missing are passed through to 53b.
usage: sweep53b.py <53 run dir> <out> <func> <overlays> <tag> <outdir> --runs <runs of record tsv> [--gaps <bed>] [--perms 50]"""
import sys, subprocess, csv, os, argparse
ap = argparse.ArgumentParser()
for p in ('run', 'out', 'func', 'overlays', 'tag', 'dest'): ap.add_argument(p)
ap.add_argument('--runs', required=True, help='catalogue runs of record (run mode; the unit of record)')
ap.add_argument('--gaps'); ap.add_argument('--perms', default='50'); ap.add_argument('--order2-sims', default='20'); ap.add_argument('--allow-missing')
a = ap.parse_args(); RUN, DEST = a.run, a.dest
HERE = os.path.dirname(os.path.abspath(__file__))
SET = {'default': [], 'major_p_0.10': ['--major-p', '0.10'], 'major_p_0.25': ['--major-p', '0.25'], 'major_n_3': ['--major-n', '3'],
       'major_n_10': ['--major-n', '10'], 'intrinsic_0.2': ['--intrinsic-frac', '0.2'], 'intrinsic_0.5': ['--intrinsic-frac', '0.5'],
       'uniform_0.9': ['--uniform-frac', '0.9'], 'uniform_0.99': ['--uniform-frac', '0.99'], 'copy_scale_1.5': ['--copy-scale', '1.5'], 'copy_scale_4': ['--copy-scale', '4']}
os.makedirs(DEST, exist_ok=True); res = {}
common = ['--out', a.out, '--run', RUN, '--func', a.func, '--overlays', a.overlays, '--genome-tag', a.tag, '--perms', a.perms, '--order2-sims', a.order2_sims, '--runs', a.runs]
if a.gaps: common += ['--gaps', a.gaps]
if a.allow_missing: common += ['--allow-missing', a.allow_missing]
for name, extra in SET.items():
    d = os.path.join(RUN, f'53b_sweep_{name}')
    r = subprocess.run([sys.executable, '-P', os.path.join(HERE, '53b_topology_compare.py')] + common + ['--dest', d] + extra)
    if r.returncode: sys.exit(f'sweep53b: setting {name} failed (53b exit {r.returncode}, traceback above); nothing overwritten, stopping')
    R = lambda f: list(csv.DictReader(open(os.path.join(d, f)), delimiter='\t'))
    edges = {(r['array_id'], r['edge']): r['role'] for r in R('53b_core_accessory_edges.tsv')}
    intr = {r['array_id']: r['intrinsic_annotations'] for r in R('53b_topologies.tsv') if r['unit'] == 'array'}   # every run (CA24)
    si = {(r['array_id'], r['unit']): r['var_structure_index'] for r in R('53b_topologies.tsv')}
    lib_rows = R('53b_topology_library.tsv'); libcol = next((c for c in (lib_rows[0].keys() if lib_rows else []) if c.startswith('type_at')), None)
    lib = {(r['array_id'], r['unit']): r[libcol] for r in lib_rows} if libcol else {}
    res[name] = (edges, intr, si, lib, libcol or 'NA'); print(name, 'done ->', d, flush=True)
ref = res['default']; rows = []
sign = lambda v: 'NA' if v in ('NA', '') else (1 if float(v) > 0.05 else -1 if float(v) < -0.05 else 0)
for name, (edges, intr, si, lib, libcol) in res.items():
    ch = sum(1 for k, v in edges.items() if k in ref[0] and ref[0][k] != v); new = sum(1 for k in edges if k not in ref[0]); gone = sum(1 for k in ref[0] if k not in edges)
    ich = sum(1 for k, v in intr.items() if ref[1].get(k) != v)
    sch = sum(1 for k, v in si.items() if k in ref[2] and sign(v) != sign(ref[2][k]))
    rows.append([name, len(edges), sum(1 for v in edges.values() if v == 'core'), ch, new, gone, ich, sch, len(set(lib.values())), libcol])
with open(os.path.join(DEST, 'sweep53b_summary.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['setting', 'n_major_edges', 'n_core', 'role_changed', 'edges_new', 'edges_lost', 'runs_intrinsic_changed', 'units_SI_sign_changed', 'library_types_first_null_cut', 'library_column'])
    w.writerows(rows)
for r in rows: print('\t'.join(map(str, r)))
