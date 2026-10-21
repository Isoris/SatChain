import sys, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d = pd.read_csv(sys.argv[1], sep='\t'); png = sys.argv[2]; title = sys.argv[3]
INK, INK2, GRID = '#1f1f1e', '#5f5e58', '#e4e2dc'
RAMP = {0.9: '#0d366b', 0.8: '#2a78d6', 0.7: '#6da7ec', 0.6: '#b7d3f6'}   # one-hue sequential, strict = dark
plt.rcParams.update({'font.family': 'Aptos', 'font.size': 10, 'axes.edgecolor': GRID, 'xtick.color': INK2, 'ytick.color': INK2})
comps = list(dict.fromkeys(d.comparison)); fig, axs = plt.subplots(2, len(comps), figsize=(15, 6.6), sharex=True, sharey='row')
for j, c in enumerate(comps):
    q = d[d.comparison == c]; ax = axs[0, j]
    for p, g in q.groupby('p_in'):
        ax.plot(g.k, g.diag_total, color=RAMP[p], lw=2, marker='o', ms=4, label=f'in >= {p:.0%}, out <= {1 - p:.0%}' if p != 0.6 else 'in >= 60%, out <= 40%')
    ax.plot(q[q.p_in == 0.8].k, q[q.p_in == 0.8].null_max_total, color=INK2, lw=1, ls='--', label='shuffled labels (max of 5)')
    ax.set_title(f"{c.replace(' of boundary', '')}\n(n = {q.n_a.iloc[0]:,} / {q.n_b.iloc[0]:,} copies)", fontsize=10.5, color=INK, loc='left')
    a2 = axs[1, j]; g = q[q.p_in == 0.8]
    a2.plot(g.k, g.best_kmer_balanced_acc, color='#eb6834', lw=2, marker='o', ms=4, label='observed')
    a2.plot(g.k, g.null_max_acc, color=INK2, lw=1, ls='--', label='shuffled labels')
    a2.set_ylim(0.45, 1.0); a2.set_xlabel('k-mer length')
    for ax_ in (ax, a2):
        ax_.grid(color=GRID, lw=0.6); [ax_.spines[s].set_visible(False) for s in ('top', 'right')]
axs[0, 0].set_ylabel('diagnostic k-mers\n(both sides)', color=INK2); axs[1, 0].set_ylabel('best single k-mer\nbalanced accuracy', color=INK2)
axs[0, 2].legend(frameon=False, fontsize=8.5, loc='upper left'); axs[1, 2].legend(frameon=False, fontsize=8.5, loc='center right')
fig.suptitle(title, x=0.01, ha='left', fontsize=13, fontweight='bold', color=INK)
fig.subplots_adjust(left=0.07, right=0.99, top=0.84, bottom=0.09, hspace=0.25, wspace=0.08)
fig.savefig(png, dpi=150, facecolor='white'); print('wrote', png)
