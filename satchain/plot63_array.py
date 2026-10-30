#!/usr/bin/env python3
"""plot63_array.py: tracks of 63_expansion_tracks.py for one or more arrays, side by side.
Per array: unit haplotypes (broad, fine), exact-match units in a window of 10, then the self-identity matrix (MinHash Jaccard of bins).
usage: plot63_array.py <dir> <out.png> name1[:title] [name2[:title] ...]"""
import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
D, OUT, names = sys.argv[1], sys.argv[2], sys.argv[3:]
PAL = ['#4c72b0', '#dd8452', '#55a868', '#c44e52', '#8172b3', '#937860', '#da8bc3', '#8c8c8c']
fig = plt.figure(figsize=(7.5 * len(names), 10.5)); fig.patch.set_facecolor('white')
for c, nt in enumerate(names):
    name, title = (nt.split(':', 1) + [nt])[:2]
    U = pd.read_csv(os.path.join(D, f'63_units_{name}.tsv'), sep='\t'); M = np.load(os.path.join(D, f'63_dotplot_{name}.npy'))
    B = pd.read_csv(os.path.join(D, f'63_dotplot_bins_{name}.tsv'), sep='\t'); lo, hi = B.start.min(), B.end.max()
    hk = [x for x in U.columns if x.startswith('hap_')]
    x0 = 0.08 + 0.91 * c / len(names); w = 0.80 / len(names)
    rows = hk + ['exact_win']; th = 0.035
    for r, col in enumerate(rows):
        ax = fig.add_axes([x0, 0.88 - r * (th + 0.01), w, th])
        if col == 'exact_win':
            for u in U.itertuples(): ax.axvspan((u.start - lo) / 1e3, (u.end - lo) / 1e3, color=plt.cm.YlOrRd(u.exact_win / 10), lw=0)
            lab = 'exact-match units\nin window of 10'
        else:
            for u in U.itertuples(): ax.axvspan((u.start - lo) / 1e3, (u.end - lo) / 1e3, color=PAL[int(getattr(u, col)) % len(PAL)], lw=0)
            lab = f'unit haplotypes\nK = {col.split("_")[1]}'
        ax.set_xlim(0, (hi - lo) / 1e3); ax.set_yticks([]); ax.set_xticks([])
        if c == 0: ax.set_ylabel(lab, rotation=0, ha='right', va='center', fontsize=11)
        if r == 0: ax.set_title(title, fontsize=13, fontweight='bold', loc='left', pad=28)
        if r == 0: ax.text(0, 1.15, f'{len(U)} units, {U.exact.mean():.0%} with an exact copy elsewhere in the array', transform=ax.transAxes, fontsize=10.5)
    ax = fig.add_axes([x0, 0.08, w, 0.88 - len(rows) * (th + 0.01) - 0.10])
    ax.imshow(np.sqrt(M), cmap='Greys', origin='upper', extent=[0, (hi - lo) / 1e3, (hi - lo) / 1e3, 0], vmin=0, vmax=1, interpolation='nearest')
    ax.set_xlabel('kb along the array', fontsize=11); ax.tick_params(labelsize=10)
    if c == 0: ax.set_ylabel('kb along the array', fontsize=11)
cax = fig.add_axes([0.30, 0.025, 0.4, 0.012])
fig.colorbar(plt.cm.ScalarMappable(cmap='Greys'), cax=cax, orientation='horizontal').set_label('self-identity: sqrt(MinHash Jaccard, k = 21) between bins', fontsize=10.5)
fig.savefig(OUT, dpi=130); print('ok')
