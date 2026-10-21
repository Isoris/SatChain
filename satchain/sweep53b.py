#!/usr/bin/env python3
"""One-at-a-time sweep of the 53b cutoffs on one step-53 run: major-edge probability / count, intrinsic-annotation fraction.
Per setting vs default: edges whose core/accessory role changes, arrays whose intrinsic annotations change, segments whose
variant structure index changes sign, library types at D 0.05. usage: sweep53b.py <53 run dir> <out> <func> <overlays> <tag> <outdir>"""
import sys, subprocess, csv, os, shutil, collections
RUN, OUT, FUNC, OVL, TAG, DEST = sys.argv[1:7]
HERE = os.path.dirname(os.path.abspath(__file__))
SET = {'default': [], 'major_p_0.10': ['--major-p', '0.10'], 'major_p_0.25': ['--major-p', '0.25'], 'major_n_3': ['--major-n', '3'],
       'major_n_10': ['--major-n', '10'], 'intrinsic_0.2': ['--intrinsic-frac', '0.2'], 'intrinsic_0.5': ['--intrinsic-frac', '0.5'],
       'uniform_0.9': ['--uniform-frac', '0.9'], 'uniform_0.99': ['--uniform-frac', '0.99'], 'copy_scale_1.5': ['--copy-scale', '1.5'], 'copy_scale_4': ['--copy-scale', '4']}
os.makedirs(DEST, exist_ok=True); res = {}
for name, extra in SET.items():
    subprocess.run([sys.executable, os.path.join(HERE, '53b_topology_compare.py'), '--out', OUT, '--run', RUN, '--func', FUNC, '--overlays', OVL,
                    '--genome-tag', TAG, '--perms', '50'] + extra, check=True, capture_output=True)
    d = os.path.join(DEST, name); shutil.rmtree(d, ignore_errors=True); shutil.copytree(os.path.join(RUN, '53b'), d)
    R = lambda f: list(csv.DictReader(open(os.path.join(d, f)), delimiter='\t'))
    edges = {(r['array_id'], r['edge']): r['role'] for r in R('53b_core_accessory_edges.tsv')}
    intr = {r['array_id']: r['intrinsic_annotations'] for r in R('53b_core_accessory_edges.tsv')}
    si = {(r['array_id'], r['unit']): r['var_structure_index'] for r in R('53b_topologies.tsv')}
    lib = {(r['array_id'], r['unit']): r['type_at_D0.05'] for r in R('53b_topology_library.tsv')}
    res[name] = (edges, intr, si, lib); print(name, 'done', flush=True)
ref = res['default']; rows = []
sign = lambda v: 'NA' if v in ('NA', '') else (1 if float(v) > 0.05 else -1 if float(v) < -0.05 else 0)
for name, (edges, intr, si, lib) in res.items():
    ch = sum(1 for k, v in edges.items() if k in ref[0] and ref[0][k] != v); new = sum(1 for k in edges if k not in ref[0]); gone = sum(1 for k in ref[0] if k not in edges)
    ich = sum(1 for k, v in intr.items() if ref[1].get(k) != v)
    sch = sum(1 for k, v in si.items() if k in ref[2] and sign(v) != sign(ref[2][k]))
    rows.append([name, len(edges), sum(1 for v in edges.values() if v == 'core'), ch, new, gone, ich, sch, len(set(lib.values()))])
with open(os.path.join(DEST, 'sweep53b_summary.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['setting', 'n_major_edges', 'n_core', 'role_changed', 'edges_new', 'edges_lost', 'arrays_intrinsic_changed', 'units_SI_sign_changed', 'library_types_D0.05'])
    w.writerows(rows)
for r in rows: print('\t'.join(map(str, r)))
