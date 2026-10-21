#!/usr/bin/env python3
"""One array, its window features and segmentation (step 53). usage: plot53_array.py <out folder> <test53 dir> <array_id> <png>"""
import sys, csv, gzip, collections
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT, T, AID, PNG = sys.argv[1:5]
csv.field_size_limit(1 << 30)
INK, INK2, GRID, SURF = '#1f1f1e', '#5f5e58', '#e4e2dc', '#ffffff'
CAT = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4']; OTHER = '#b4b2a9'
PROF = {'homogeneous_monomeric': '#2a78d6', 'variant_mixed': '#1baf7a', 'alternating': '#6250d6', 'interleaved_families': '#eb6834', 'unresolved': '#b4b2a9'}
NAMES = {}
import os
if os.environ.get('SATCHAIN_NAMES'):   # optional family -> display-name table (pipeline_family_id, manuscript_name)
    for r in csv.DictReader(open(os.environ['SATCHAIN_NAMES'], encoding='utf-8'), delimiter='\t'): NAMES[r['pipeline_family_id']] = r['manuscript_name']

win = [r for r in csv.DictReader(gzip.open(f'{T}/53_window_features.tsv.gz', 'rt'), delimiter='\t') if r['array_id'] == AID]
seg = [r for r in csv.DictReader(open(f'{T}/53_segments.tsv'), delimiter='\t') if r['array_id'] == AID]
cps = [r for r in csv.DictReader(open(f'{T}/53_change_points.tsv'), delimiter='\t') if r['array_id'] == AID]
cop = []
with open(f'{OUT}/03_monomer_copies.tsv', newline='') as fh:
    rd = csv.reader(fh, delimiter='\t'); h = next(rd); ix = {k: i for i, k in enumerate(h)}
    for row in rd:
        if row[ix['array_id']] == AID:
            try: idt = float(row[ix['best_alignment_identity']])
            except ValueError: idt = 0
            cop.append(((int(row[ix['genomic_start_0based']]) + int(row[ix['genomic_end']])) / 2e6, row[ix['monomer_class']] + row[ix['strand']], idt,
                        row[ix['complete_or_truncated']] == 'complete', row[ix['family_id']]))
chrom = win[0]['chrom']; period = float(seg[0]['period_02'])
x = [(int(r['start0']) + int(r['end'])) / 2e6 for r in win]
f = lambda c: [float(r[c]) if r[c] != 'NA' else float('nan') for r in win]
conf = [c for c in cop if c[3] and c[2] >= 80]
top = [k for k, _ in collections.Counter(c[1] for c in conf).most_common(5)]
fam = collections.Counter(c[4] for c in conf).most_common(1)[0][0]

plt.rcParams.update({'font.family': 'Aptos', 'font.size': 10, 'axes.edgecolor': GRID, 'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2})
panels = [('copies', 'monomer\nclasses'), ('local_period', 'local period (bp)'), ('period_strength', 'periodicity\nstrength'),
          ('GC', 'GC fraction'), ('kmer_pc1', '4-mer profile\nPC1'), ('mean_identity', 'copy identity (%)'), ('selfloop', 'class self-loop\nfraction')]
fig, axs = plt.subplots(len(panels) + 1, 1, figsize=(13, 11.5), sharex=True, gridspec_kw={'height_ratios': [0.45, 1.6] + [1] * (len(panels) - 1)})
fig.patch.set_facecolor(SURF)
# segment band
ax = axs[0]
for s in seg:
    a0, a1 = int(s['start0']) / 1e6, int(s['end']) / 1e6
    ax.add_patch(Rectangle((a0, 0), a1 - a0 - 0.0008, 1, color=PROF[s['profile']], lw=0))
    if a1 - a0 > 0.05 * (int(seg[-1]['end']) - int(seg[0]['start0'])) / 1e6: ax.text((a0 + a1) / 2, 0.5, f"seg {s['segment']}: {s['profile'].replace('_', ' ')}", ha='center', va='center', fontsize=9, color='white' if s['profile'] != 'unresolved' else INK)
ax.set_ylim(0, 1); ax.set_yticks([]); ax.set_ylabel('segments', rotation=0, ha='right', va='center')
for k, (col, lab) in enumerate(panels):
    ax = axs[k + 1]
    if col == 'copies':
        for i, cl in enumerate(top + ['other']):
            pts = [c for c in conf if (c[1] == cl) or (cl == 'other' and c[1] not in top)]
            ax.scatter([c[0] for c in pts], [i] * len(pts), s=4, marker='|', color=CAT[i] if cl != 'other' else OTHER, lw=0.8)
        unc = [c for c in cop if not (c[3] and c[2] >= 80)]
        ax.scatter([c[0] for c in unc], [len(top) + 1] * len(unc), s=4, marker='|', color=OTHER, lw=0.6)
        ax.set_yticks(range(len(top) + 2)); ax.set_yticklabels(top + ['other class', 'uncertain copy']); ax.invert_yaxis()
    else:
        y = f(col); ax.plot(x, y, color=INK, lw=1.4)
        if col == 'local_period': ax.axhline(period, color=CAT[1], lw=1, ls='--'); ax.text(x[0], period, f' array period {period:.0f} bp', color=INK2, va='bottom', fontsize=8.5)
    ax.set_ylabel(lab, rotation=0, ha='right', va='center'); ax.grid(axis='y', color=GRID, lw=0.6)
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
for ax in axs:
    for c in cps:
        ax.axvline(int(c['position']) / 1e6, color=INK if c['status'] == 'accepted' else INK2, lw=1.1 if c['status'] == 'accepted' else 0.7,
                   ls='-' if c['status'] == 'accepted' else ':', zorder=0)
axs[-1].set_xlabel(f'{chrom} position (Mb)')
fig.suptitle(f'{AID}, {chrom} centromere: {NAMES.get(fam, fam)} ({fam}), {len(cop):,} copies; window {int(win[0]["end"]) - int(win[0]["start0"]):,} bp, step half a window',
             x=0.01, ha='left', fontsize=13, color=INK, fontweight='bold')
fig.text(0.01, 0.945, 'solid line = accepted boundary (>= 2 feature blocks, >= 1 sequence-only block); dotted = change point in one block only', fontsize=9.5, color=INK2)
fig.subplots_adjust(left=0.15, right=0.99, top=0.93, bottom=0.05, hspace=0.18)
fig.savefig(PNG, dpi=150, facecolor=SURF)
print('wrote', PNG)
