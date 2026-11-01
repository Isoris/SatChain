#!/usr/bin/env python3
"""plot63_array.py: tracks of 63_expansion_tracks.py for one or more arrays, side by side.
Per array: unit haplotypes (broad, fine), exact-match units in a window (share of the units actually in the window: the last windows are
truncated, 63 writes win_units), then the self-identity matrix (MinHash Jaccard of bins).
2026-11-01 (D26): saved at 450 dpi plus a PDF next to it (project rule); haplotype colours are never recycled (K > 8 draws from tab20) and
every colour is in a key (haplotype legend per K, colour bar for the window share).
usage: plot63_array.py <dir> <out.png> name1[:title] [name2[:title] ...]"""
import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
D, OUT, names = sys.argv[1], sys.argv[2], sys.argv[3:]
PAL = ['#4c72b0', '#dd8452', '#55a868', '#c44e52', '#8172b3', '#937860', '#da8bc3', '#8c8c8c']
def palette(K): return PAL[:K] if K <= len(PAL) else [matplotlib.colors.to_hex(plt.cm.tab20(i / 20)) for i in range(K)] if K <= 20 else [matplotlib.colors.to_hex(plt.cm.gist_ncar(i / K)) for i in range(K)]
fig = plt.figure(figsize=(7.5 * len(names), 10.5)); fig.patch.set_facecolor('white')
for c, nt in enumerate(names):
    name, title = (nt.split(':', 1) + [nt])[:2]
    U = pd.read_csv(os.path.join(D, f'63_units_{name}.tsv'), sep='\t'); M = np.load(os.path.join(D, f'63_dotplot_{name}.npy'))
    B = pd.read_csv(os.path.join(D, f'63_dotplot_bins_{name}.tsv'), sep='\t'); lo, hi = B.start.min(), B.end.max()
    hk = [x for x in U.columns if x.startswith('hap_')]
    if 'win_units' not in U.columns: U['win_units'] = 10          # files from before 2026-11-01
    x0 = 0.08 + 0.91 * c / len(names); w = 0.80 / len(names)
    rows = hk + ['exact_win']; th = 0.035
    for r, col in enumerate(rows):
        ax = fig.add_axes([x0, 0.88 - r * (th + 0.01), w, th])
        if col == 'exact_win':
            for u in U.itertuples(): ax.axvspan((u.start - lo) / 1e3, (u.end - lo) / 1e3, color=plt.cm.YlOrRd(u.exact_win / max(u.win_units, 1)), lw=0)
            lab = f'exact-match units\nin window of {int(U.win_units.max())} (share)'
        else:
            K = int(col.split('_')[1]); pal = palette(max(K, int(U[col].max()) + 1))
            for u in U.itertuples(): ax.axvspan((u.start - lo) / 1e3, (u.end - lo) / 1e3, color=pal[int(getattr(u, col))], lw=0)
            lab = f'unit haplotypes\nK = {K}'
            ax.legend(handles=[Patch(color=pal[j], label=f'hap {j}') for j in sorted(U[col].unique())], loc='center left', bbox_to_anchor=(1.0, 0.5),
                      fontsize=7, frameon=False, ncol=1 if K <= 6 else 2, handlelength=1, borderaxespad=0.2)
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
cax2 = fig.add_axes([0.75, 0.025, 0.15, 0.012])
fig.colorbar(plt.cm.ScalarMappable(cmap='YlOrRd'), cax=cax2, orientation='horizontal').set_label('share of exact-match units in window', fontsize=9)
fig.savefig(OUT, dpi=450); fig.savefig(os.path.splitext(OUT)[0] + '.pdf'); print('ok', OUT, os.path.splitext(OUT)[0] + '.pdf')
