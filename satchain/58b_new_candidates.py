#!/usr/bin/env python3
"""58b: final list of NEW tandem candidates (scan 57, tier p99.9, outside every step-02 array, lag >= 7 bp).
Each candidate (first 20 kb tested): recurrence excess at its lag and 3-copy chain fraction, both against a dinucleotide-block shuffle;
  new_array   chain excess > 0.01 (>= 3 tandem copies)          new_pair  recurrence excess > 0.02 but no chain (duplicated pair)
  weak        neither
plus TRF support (same unit, >= 3 copies) and TE overlap share (EDTA). Output 58b_new_candidates.bed (BED9, coloured by tier) and a summary.
usage: 58b_new_candidates.py <tag> <cal58 dir> <func genome folder> <pipeline out> <genome fasta>"""
import sys, importlib.util, collections, gzip, os
import numpy as np
TAG, CAL, FUNC, OUT, FA = sys.argv[1:6]
sys.argv = [sys.argv[0], TAG, CAL, '0.999', FUNC, OUT, FA, '0']       # reuse val58_new helpers without running its report
spec = importlib.util.spec_from_file_location('v', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'val58_new.py'))
src = open(spec.origin, encoding='utf-8').read().split("print(f'{TAG} tier p{TIER}")[0]
V = {}; exec(compile(src, spec.origin, 'exec'), V)
te = collections.defaultdict(list)
for l in gzip.open(f'{FUNC}/_repeats_EDTA/fClaHyb_{TAG}.TEanno.sorted.gff3.gz', 'rt'):
    if l.startswith('#'): continue
    f = l.split('\t', 5)
    if len(f) > 4 and f[2] not in ('repeat_region', 'target_site_duplication', 'long_terminal_repeat'): te[f[0]].append((int(f[3]) - 1, int(f[4])))
for v in te.values(): v.sort()
COL = {'new_array': '31,119,180', 'new_pair': '235,104,52', 'weak': '180,178,169'}
tot = collections.Counter(); n = collections.Counter(); trf_ok = collections.Counter()
with open(f'{CAL}/58b_new_candidates.bed', 'w') as bed:
    for cls, v in V['classes'].items():
        for c, s, e, lag, out in v:
            sq = V['CODE'][np.frombuffer(V['fetch'](c, s, min(e, s + 20000)), dtype=np.uint8)]; sh = V['dinuc_shuffle'](sq)
            r = V['rec_at'](sq, lag) - V['rec_at'](sh, lag); ch = V['chain_at'](sq, lag) - V['chain_at'](sh, lag)
            tier = 'new_array' if ch > 0.01 else 'new_pair' if r > 0.02 else 'weak'
            tr = V['trf_support'](c, s, e, lag); tef = min(1.0, V['ov'](te, c, s, e) / (e - s))
            k = (cls, tier, 'TE' if tef > 0.5 else 'nonTE'); tot[k] += out; n[k] += 1; trf_ok[k] += tr
            bed.write(f'{c}\t{s}\t{e}\t{tier}|unit{lag}|{cls}|TE{tef:.2f}|TRF{int(tr)}\t0\t.\t{s}\t{e}\t{COL[tier]}\n')
print(f'{TAG}: class\ttier\tTE\tn\tMb\tTRF-supported')
for k in sorted(tot): print(f'{TAG}\t' + '\t'.join(k) + f'\t{n[k]}\t{tot[k] / 1e6:.2f}\t{trf_ok[k] / n[k]:.2f}')
arr = sum(v for k, v in tot.items() if k[1] == 'new_array'); arr_nte = sum(v for k, v in tot.items() if k[1] == 'new_array' and k[2] == 'nonTE')
print(f'{TAG}: new arrays (>= 3 copies) {arr / 1e6:.2f} Mb, of which non-TE {arr_nte / 1e6:.2f} Mb')
