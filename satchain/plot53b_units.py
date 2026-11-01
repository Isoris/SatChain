#!/usr/bin/env python3
"""Big topology (whole run) vs subtopologies (segments): transition-probability heatmaps on the run's states, plus memory,
structure index and divergence from the whole run. 2026-11-01 (HB11): the 53b unit named 'array' is one run of record in run mode;
--runs <runs tsv> labels it so (interval, bp); a step-02 interval needs --extents-ok. PNG at --dpi (default 450) plus a vector PDF.
usage: plot53b_units.py <53b dir> <run_id> <png> [--min-state-count 5] --runs <runs tsv> [--dpi 450]"""
import sys, csv, collections, argparse
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

ap = argparse.ArgumentParser(); ap.add_argument('d53b'); ap.add_argument('aid'); ap.add_argument('png')
ap.add_argument('--min-state-count', type=int, default=5); ap.add_argument('--runs'); ap.add_argument('--dpi', type=int, default=450); ap.add_argument('--extents-ok', action='store_true')
A = ap.parse_args(); D, AID, PNG, MINC = A.d53b, A.aid, A.png, A.min_state_count
if not A.runs and not A.extents_ok: ap.error('--runs <runs of record> (run_id) is the unit; a step-02 interval needs --extents-ok (CA43)')
rrow = next((r for r in csv.DictReader(open(A.runs), delimiter='\t') if r['run_id'] == AID), None) if A.runs else None
if A.runs and rrow is None: sys.exit(f'{AID} not in {A.runs}')
UNIT = 'run' if rrow else 'step-02 interval'
INK, INK2, GRID = '#1f1f1e', '#5f5e58', '#e4e2dc'
SEQ = LinearSegmentedColormap.from_list('blue', ['#ffffff', '#cde2fb', '#86b6ef', '#3987e5', '#1c5cab', '#0d366b'])
top = [r for r in csv.DictReader(open(f'{D}/53b_topologies.tsv'), delimiter='\t') if r['array_id'] == AID]
tr = collections.defaultdict(collections.Counter)
for r in csv.DictReader(open(f'{D}/53b_transitions.tsv'), delimiter='\t'):
    if r['array_id'] == AID: tr[r['unit']][(r['from'], r['to'])] += int(r['count'])
comp = collections.Counter()
for (x, y), v in tr['array'].items(): comp[x] += v; comp[y] += v
order = lambda s: (0 if s.startswith('V') else 1 if s in ('U', 'OF') else 2, s)
states = sorted([s for s, n in comp.items() if n >= MINC], key=order)
units = [r for r in top if int(r['n_trans'] or 0) >= 10]
plt.rcParams.update({'font.family': 'Aptos', 'font.size': 9.5, 'axes.edgecolor': GRID, 'xtick.color': INK2, 'ytick.color': INK2})
fig, axs = plt.subplots(2, len(units), figsize=(3.1 * len(units) + 1, 7.4), gridspec_kw={'height_ratios': [3.2, 1]}, squeeze=False)
for j, u in enumerate(units):
    t = tr[u['unit']]; row = collections.Counter()
    for (x, y), v in t.items(): row[x] += v
    M = [[t[(x, y)] / row[x] if row[x] else float('nan') for y in states] for x in states]
    ax = axs[0, j]; im = ax.imshow(M, cmap=SEQ, vmin=0, vmax=1)
    for a_, x in enumerate(states):
        for b_, y in enumerate(states):
            if row[x] and t[(x, y)]: ax.text(b_, a_, f'{M[a_][b_]:.2f}', ha='center', va='center', fontsize=7.5, color='white' if M[a_][b_] > 0.55 else INK)
        if not row[x]: ax.add_patch(plt.Rectangle((-0.5, a_ - 0.5), len(states), 1, color='#f2f1ed', lw=0))
    ax.set_xticks(range(len(states))); ax.set_xticklabels(states, rotation=55, ha='right', fontsize=8)
    ax.set_yticks(range(len(states))); ax.set_yticklabels(states if j == 0 else [''] * len(states), fontsize=8)
    lab = f'whole {UNIT} (big topology)' if u['unit'] == 'array' else f"{u['unit']}: {u['segment_profile'].replace('_', ' ')}"
    ax.set_title(f"{lab}\n{int(u['start0']) / 1e6:.3f}-{int(u['end']) / 1e6:.3f} Mb, {int(u['n_trans']):,} transitions", fontsize=9.5, color=INK)
    if j == 0: ax.set_ylabel('current state', color=INK2)
    ax.set_xlabel('next state', color=INK2)
    b = axs[1, j]; vals = [('memory (bits)', u['memory_MI_bits']), ('structure index', u['structure_index']), ('D vs array', u['D_trans_vs_array'])]
    for k, (nm, v) in enumerate(vals):
        v = float(v) if v not in ('NA', '') else None
        if v is not None: b.barh(k, v, color='#2a78d6' if nm != 'structure index' else ('#2a78d6' if v >= 0 else '#eb6834'), height=0.6)
        b.text((max(v, 0) if v is not None else 0) + 0.05, k, f'{nm}: ' + ('NA' if v is None else f'{v:.3f}'), ha='left', va='center', fontsize=8.5, color=INK)
    b.axvline(0, color=INK2, lw=0.8); b.set_xlim(-1, 2.2); b.set_ylim(2.6, -0.6); b.set_yticks([])
    for s_ in ('top', 'right', 'left'): b.spines[s_].set_visible(False)
span = f", {rrow['chrom']}:{int(rrow['start0']):,}-{int(rrow['end']):,}, {int(rrow['end']) - int(rrow['start0']):,} bp" if rrow else ''
fig.suptitle(f'{AID} ({UNIT}{span}): big topology vs subtopologies (P(next state | current state); grey row = state absent)', x=0.01, ha='left', fontsize=12.5, fontweight='bold', color=INK)
fig.text(0.01, 0.925, 'states: V1, V2 = monomer variants by rank (after class merge); .SINE = copy carrying the intrinsic AnnoSINE call; U = uncertain copy; '
         'TE:*, 5S, tRNA, CDS, INS = annotation / insertion nodes; an assembly gap is a path break, not a node. Structure index: +1 blocks, 0 no memory, < 0 alternation.', fontsize=8.8, color=INK2)
fig.subplots_adjust(left=0.07, right=0.99, top=0.84, bottom=0.06, wspace=0.12, hspace=0.45)
fig.savefig(PNG, dpi=A.dpi, facecolor='white'); fig.savefig(PNG.rsplit('.', 1)[0] + '.pdf', facecolor='white')   # 2026-11-01: >= 450 dpi + vector PDF
print('wrote', PNG, 'and', PNG.rsplit('.', 1)[0] + '.pdf')
