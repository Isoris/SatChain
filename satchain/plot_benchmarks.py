#!/usr/bin/env python3
"""plot_benchmarks.py (2026-11-01): one figure per decision of a SatChain run, drawn from the shipped tables, so every rule can be inspected.
Panels (PNG 450 dpi + PDF, <dest>/bench_<name>.*):
  k_choice        ARI across k, ARI across restarts, silhouette vs K (55_K_choice.tsv) with the rule lines; cluster sizes at each K (55_cluster_summary.tsv)
  gene_arrays     every family with a gene hit on the bp-ratio x unit-share plane, cut lines, named families (59b_gene_arrays_<G>.tsv)
  gates           runs passing all gates over the 27 sweep settings (59_gates_sweep_<G>.tsv)
  val56           val56 benchmark: runs kept / Mb added / fraction of v4.5 kept vs q (val56_benchmark_<G>.tsv)
  band_test       tandem share by lag-identity band vs the null (threshold_band_test_<G>_summary.tsv)
  interruptions   class Mb over the 59e sweep (59e_sweep_<G>.tsv)
  order2          order-2 hits vs chain length and q cut (53b_o2_sweep.tsv)
  te_layer        TE insertions inside runs: two carry rules x three intrinsic cuts (59c_te_<G>.tsv in cat59 and cat59_minrule)
usage: plot_benchmarks.py <run folder> <dest> [--tags Gar,Mac]
Numbered-story layout (2026-11-01, user: "it is the best for understanding, it's for noobs"): every figure carries a top banner "N  <one-line conclusion>" numbered 1-8 in
pipeline order (1 band test, 2 val56, 3 gates, 4 gene arrays, 5 TE layer, 6 interruptions, 7 K choice, 8 order 2); the conclusion is written from the tables at render
time; the per-panel titles stay the technical line.  Banner contract = figures/paper_figures/code/assemble_v2.py: BAN 7 mm, bold number + bold sentence (8 pt), 0.5-pt
navy hairline across the figure, house style satellite_style.py (Arial, navy #17252D, no top / right spines)."""
import os, sys, argparse, csv, collections, textwrap
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ap = argparse.ArgumentParser(); ap.add_argument('V'); ap.add_argument('DEST'); ap.add_argument('--tags', default='Gar,Mac'); a = ap.parse_args()
V, D, TAGS = a.V, a.DEST, a.tags.split(','); os.makedirs(D, exist_ok=True)
NAVY = '#17252D'; INK, INK2, GRID = NAVY, '#5f5e58', '#e4e2dc'; HAP = {'Gar': '#c0392b', 'Mac': '#1f5fa8'}; SP = {'Gar': 'C. gariepinus', 'Mac': 'C. macrocephalus'}
plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.edgecolor': GRID, 'xtick.color': INK2, 'ytick.color': INK2, 'axes.labelcolor': INK2, 'text.color': NAVY,
                     'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
MM = 1 / 25.4; BAN = 7.0   # banner height (mm) for one line; +3.5 mm per extra wrapped line


def tidy(ax, title=None):
    ax.grid(color=GRID, lw=0.5); [ax.spines[s].set_visible(False) for s in ('top', 'right')]
    if title: ax.set_title(title, loc='left', fontsize=9, color=INK)


def banner(fig, num, text):
    """top banner of the numbered-story layout: hairline at the top edge, bold number at 0.5 mm, bold sentence at 5 mm (assemble_v2.banner), wrapped to the figure width;
    the panels are packed below it with tight_layout(rect=...) so the banner never covers data"""
    W, H = fig.get_size_inches() / MM
    lines = textwrap.wrap(text, width=max(40, int((W - 8) / 1.42)))           # ~1.42 mm per character at 8 pt bold Arial (measured on the 2026-11-01 renders)
    ban = BAN + 3.5 * (len(lines) - 1)
    fig.add_artist(Line2D([0, 1], [1, 1], color=NAVY, lw=0.5, transform=fig.transFigure))
    fig.text(0.5 / W, 1 - 1.0 / H, num, fontsize=8, fontweight='bold', ha='left', va='top', color=NAVY)
    fig.text(5.0 / W, 1 - 1.0 / H, '\n'.join(lines), fontsize=8, fontweight='bold', ha='left', va='top', color=NAVY, linespacing=1.25)
    fig.add_artist(Line2D([0, 1], [1 - ban / H] * 2, color=NAVY, lw=0.5, transform=fig.transFigure))
    fig.tight_layout(rect=[0, 0, 1, 1 - (ban + 1.0) / H])


def save(fig, name):
    fig.savefig(f'{D}/bench_{name}.png', dpi=450, facecolor='white', bbox_inches='tight'); fig.savefig(f'{D}/bench_{name}.pdf', bbox_inches='tight'); plt.close(fig); print('wrote', name)


def pct(x): return f'{100 * x:.0f}'


def rng(vals, f=pct):
    lo, hi = f(min(vals)), f(max(vals)); return lo if lo == hi else f'{lo}-{hi}'


def rd(p): return pd.read_csv(p, sep='\t', comment='#') if os.path.exists(p) else None


# ---- k_choice
k = rd(f'{V}/lin55/55_K_choice.tsv'); cs = rd(f'{V}/lin55/55_cluster_summary.tsv')
if k is not None:
    fig, axs = plt.subplots(1, 3, figsize=(10, 3.2))
    ax = axs[0]; ax.plot(k.K, k.min_ARI_across_k, 'o-', color='#0d366b', label='min ARI across k = 5, 6, 7'); ax.plot(k.K, k.min_ARI_restarts, 's-', color='#2a78d6', label='min ARI between restart seeds')
    ax.axhline(0.8, color=INK2, ls='--', lw=0.8); ax.set_ylim(0, 1.05); ax.set_xlabel('K'); ax.set_ylabel('adjusted Rand index'); ax.legend(frameon=False, fontsize=7); tidy(ax, 'stability of the partition (rule: both >= 0.8)')
    ax = axs[1]; ax.plot(k.K, k.mean_silhouette, 'o-', color='#eb6834'); ax.axhline(0.9 * k.mean_silhouette.max(), color=INK2, ls='--', lw=0.8); ax.set_xlabel('K'); ax.set_ylabel('mean silhouette (k = 5-7)'); tidy(ax, 'silhouette (rule: >= 0.9 x best)')
    ch = k[k.chosen == 'yes'].K.iloc[0]; axs[0].axvline(ch, color='#eb6834', lw=1); axs[1].axvline(ch, color='#eb6834', lw=1); axs[0].text(ch + 0.1, 0.05, f'chosen K = {ch}', color='#eb6834', fontsize=7)
    ax = axs[2]
    if cs is not None:
        c6 = cs[cs.k == 6]; bottom = np.zeros(len(c6.K.unique())); Ks = sorted(c6.K.unique())
        for j in range(max(Ks)):
            vals = [c6[(c6.K == K_) & (c6.cluster == f'C{K_}.{j + 1}')].n.sum() for K_ in Ks]
            ax.bar(Ks, vals, bottom=bottom, color=plt.cm.tab20(j % 20), edgecolor='white', lw=0.3, label=f'cluster {j + 1}'); bottom += np.array(vals)
        ax.set_ylim(0, bottom.max() * 1.02)   # 2026-11-09: every colour in the key, bars not clipped
        ax.set_xlabel('K'); ax.set_ylabel('sampled copies'); ax.legend(frameon=False, fontsize=6, title='cluster (by size)', title_fontsize=6.5, loc='upper left', bbox_to_anchor=(1.0, 1.0)); tidy(ax, 'cluster sizes at each K (k = 6)')
    # 2026-11-09: banner text from the table (the old text said every other K failed the ARI rule; in the final run K = 4 is stable and fails on silhouette)
    rest = k[k.K != ch]; sil_cut = 0.9 * k.mean_silhouette.max()
    ari_ok = (rest.min_ARI_across_k >= 0.8) & (rest.min_ARI_restarts >= 0.8); fail_ari = rest[~ari_ok].K.tolist(); fail_sil = rest[ari_ok & (rest.mean_silhouette < sil_cut)].K.tolist()
    why = '; '.join(t for t in [f'K = {", ".join(map(str, fail_ari))} unstable (ARI < 0.8)' if fail_ari else '',
                                 f'K = {", ".join(map(str, fail_sil))} stable but silhouette below {sil_cut:.2f}' if fail_sil else ''] if t)
    comp = ''
    if cs is not None:
        cc = cs[(cs.k == 6) & (cs.K == ch)]
        comp = '; ' + ', '.join(f'{r.cluster} = {r.n:,} copies (top {r.top_families.split(";")[0].split(":")[0]})' for r in cc.itertuples())
    banner(fig, '7', f'Step 55, K choice: the k-mer clusters of the sampled monomers of candidate centromeric arrays pass the rule only at K = {ch} '
                     f'({k[k.chosen == "yes"].note.iloc[0] or "rule"}; ARI {k[k.K == ch].min_ARI_across_k.iloc[0]:.2f}{comp}); {why}; so the clusters are a grouping, not a history')
    save(fig, 'k_choice')
# ---- gene arrays
fig, axs = plt.subplots(1, len(TAGS), figsize=(5 * len(TAGS), 4), sharey=True); ga_txt = []
for ax, G in zip(np.atleast_1d(axs), TAGS):
    g = rd(f'{V}/cat59/59b_gene_arrays_{G}.tsv')
    if g is None: continue
    isga = g.element_class.astype(str).str.startswith('gene_array'); ga_txt.append(f'{isga.sum()} of {len(g)} ({SP[G]}, {(isga & g.element_class.str.contains("tRNA")).sum()} tRNA)')
    col = g.top_class.map({'5S_rRNA': '#e6a817', 'tRNA': '#2e86c1', 'snRNA_U1': '#7d3c98', 'snRNA_U2': '#9b59b6', 'snRNA_U4': '#9b59b6', 'snRNA_U5': '#9b59b6', 'snRNA_U6': '#b07cc6', 'miR-430': '#1e8449', 'rRNA_other': '#8e5a2b'}).fillna('#999999')
    ax.scatter(g.observed_over_expected.clip(upper=3.5), g.unit_share, s=np.sqrt(g.bp / 1e3) * 2 + 4, c=col, alpha=0.75, edgecolor='white', lw=0.3)
    ax.axvline(0.7, color=INK2, ls='--', lw=0.8); ax.axhline(0.5, color=INK2, ls='--', lw=0.8); ax.axvline(0.2, color=GRID, ls=':', lw=0.8)
    for r in g[(g.bp > 150000) | g.family.isin(('MacSat06', 'GarSat49', 'GarSat04'))].itertuples(): ax.annotate(f'{r.family} ({r.top_class})', (min(r.observed_over_expected, 3.5), r.unit_share), fontsize=6, xytext=(3, 3), textcoords='offset points', color=INK)
    ax.set_xlabel('observed / expected bp share (one gene per unit)'); ax.set_ylabel('share of units with an annotated gene'); ax.set_xlim(-0.05, 3.6); ax.set_ylim(-0.03, 1.05)
    tidy(ax, f'{SP[G]}: gene array = ratio >= 0.7 AND unit share >= 0.5; area = family bp')
    for lab, c in (('5S rRNA', '#e6a817'), ('tRNA', '#2e86c1'), ('snRNA U1/U2/U4/U5/U6', '#9b59b6'), ('miR-430', '#1e8449'), ('rRNA, other', '#8e5a2b'), ('other gene class', '#999999')): ax.scatter([], [], c=c, label=lab)
    ax.legend(frameon=False, fontsize=7, loc='lower right')
banner(fig, '4', 'Step 59b, gene arrays: tandem copies of a gene are not satellites; two criteria call ' + ' and '.join(ga_txt) + ' families with gene hits gene arrays (5S rDNA = GarSat49, MacSat06), '
                 'while GarSat04 stays a satellite (a 5S segment in the unit, but no annotated gene)')
save(fig, 'gene_arrays')
# ---- gates
fig, axs = plt.subplots(1, len(TAGS), figsize=(5 * len(TAGS), 3.4), sharey=True); gt_def, gt_rng = [], []
for ax, G in zip(np.atleast_1d(axs), TAGS):
    s = rd(f'{V}/cat59/59_gates_sweep_{G}.tsv')
    if s is None: continue
    s['pct'] = 100 * s.n_pass_all / s.n_runs; x = np.arange(len(s)); gt_rng += [s.pct.min(), s.pct.max()]
    col = s.G2_identity.map({80: '#b7d3f6', 85: '#2a78d6', 90: '#0d366b'})
    ax.bar(x, s.pct, color=col, width=0.8); ax.set_ylim(0, 100); ax.set_ylabel('% runs passing G1, G2 and G3'); ax.set_xticks(x); ax.set_xticklabels([f"{r.G1_gc_points}/{r.G2_identity}/{r.G3_band}" for r in s.itertuples()], rotation=90, fontsize=5.5)
    d = s[(s.G1_gc_points == 5) & (s.G2_identity == 85) & (s.G3_band == '0.35-0.65')]
    if len(d): ax.bar(d.index, d.pct, color='#eb6834', width=0.8); ax.text(d.index[0], d.pct.iloc[0] + 2, f'default {d.pct.iloc[0]:.0f} %', fontsize=7, ha='center', color='#eb6834'); gt_def.append(f'{d.pct.iloc[0]:.0f} % ({SP[G]})')
    for v, c in ((80, '#b7d3f6'), (85, '#2a78d6'), (90, '#0d366b')): ax.bar([0], [0], color=c, label=f'G2 identity >= {v}')
    ax.legend(frameon=False, fontsize=7, loc='lower left'); tidy(ax, f'{SP[G]}: G1 GC points / G2 identity / G3 band')
banner(fig, '3', f'Step 59, QC gates: the default setting (orange) keeps {" and ".join(gt_def)} of the runs; over the 27 settings {min(gt_rng):.0f}-{max(gt_rng):.0f} % pass, '
                 'and it is the identity gate G2 (80 / 85 / 90 %) that moves the count, not the GC or band gates')
save(fig, 'gates')
# ---- val56
fig, axs = plt.subplots(1, 3, figsize=(10, 3.2)); v_kept, v_old, v_add = [], [], []
for G in TAGS:
    b = rd(f'{V}/hmm56_{G}/val56_benchmark_{G}.tsv')
    if b is None: continue
    v_kept += list(b.runs_kept / b.runs_total); v_old += list(b.frac_36_bp_kept.clip(upper=1)); v_add += list(b.added_beyond_36_Mb)
    axs[0].plot(b.q, 100 * b.runs_kept / b.runs_total, 'o-', color=HAP[G], label=SP[G]); axs[1].plot(b.q, b.added_beyond_36_Mb, 'o-', color=HAP[G]); axs[2].plot(b.q, 100 * b.frac_36_bp_kept, 'o-', color=HAP[G])
axs[0].set_ylabel('% HMM runs kept'); axs[1].set_ylabel('Mb added beyond v4.5'); axs[2].set_ylabel('% of v4.5 bp kept')
for ax, t in zip(axs, ('runs kept', 'sequence added', 'old sequence kept (summed per run; > 100 % = overlapping runs)')): ax.set_xlabel('background quantile q'); tidy(ax, t)
axs[0].legend(frameon=False, fontsize=7)
banner(fig, '2', f'val56, background check: every run is re-tested against windows outside all candidate intervals; at q = 0.8-0.95 {rng(v_kept)} % of the runs survive, '
                 f'{rng(v_old)} % of the v4.5 satellite sequence is kept and {rng(v_add, lambda v: f"{v:.1f}")} Mb is new, so the check removes little and the quantile is not a sensitive choice')
save(fig, 'val56')
# ---- band test
fig, ax = plt.subplots(figsize=(6, 3.4)); bt_cliff, bt_null = [], []
for G in TAGS:
    t = rd(f'{V}/benchmarks/threshold_band_test_{G}_summary.tsv')
    if t is None: continue
    bands = t[t.band != 'random_outside_arrays']; x = np.arange(len(bands))
    bt_cliff += list(bands[bands.band.isin(('0.30-0.32', '0.32-0.34'))].tandem_above_null_p95); bt_null += list(t[t.band == 'random_outside_arrays'].tandem_above_null_p95); ax.plot(x + (0.08 if G == 'Mac' else -0.08), 100 * bands.tandem_above_null_p95, 'o-', color=HAP[G], label=f'{SP[G]} (above null p95)')
    nullv = 100 * t[t.band == 'random_outside_arrays'].tandem_above_null_p95.iloc[0]; ax.axhline(nullv, color=HAP[G], ls=':', lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels(bands.band)
ax.axvline(2.5, color=INK2, ls='--', lw=0.8); ax.text(2.55, 8, 'old cut 0.34', fontsize=7, color=INK2); ax.set_ylabel('% windows tandem on the sequence'); ax.set_xlabel('lag-identity band of the v4.5 scan'); ax.legend(frameon=False, fontsize=7); tidy(ax, 'dotted: random windows outside arrays, runs and TRF')
banner(fig, '1', f'Band test: the old lag-identity cut (0.34) lost real satellite; windows just below it, at 0.30-0.34, are tandem on the sequence in {rng(bt_cliff)} % of cases, '
                 f'random windows in {rng(bt_null)} %')
save(fig, 'band_test')
# ---- interruptions
fig, axs = plt.subplots(1, len(TAGS), figsize=(5.5 * len(TAGS), 3.6), sharey=True); it_txt = []
for ax, G in zip(np.atleast_1d(axs), TAGS):
    s = rd(f'{V}/cat59/59e_sweep_{G}.tsv')
    if s is None: continue
    tot = s.groupby(['class_cut', 'intrinsic_cut', 'max_gap']).Mb.first(); tsh = s[s['class'] == 'TE_insertion'].groupby(['class_cut', 'intrinsic_cut', 'max_gap']).Mb_class.sum() / tot
    it_txt.append(f'{100 * tsh.min():.0f}-{100 * tsh.max():.0f} % of the {tot.min():.1f}-{tot.max():.1f} Mb in {SP[G]}')
    s['setting'] = s.class_cut.astype(str) + '/' + s.intrinsic_cut.astype(str) + '/' + (s.max_gap // 1000).astype(str) + 'k'; piv = s.pivot_table(index='setting', columns='class', values='Mb_class', aggfunc='sum').fillna(0)
    piv = piv.loc[sorted(piv.index, key=lambda v: (float(v.split('/')[0]), float(v.split('/')[1]), int(v.split('/')[2][:-1])))]
    cols = {'TE_insertion': '#e6a817', 'missed_TE_54': '#b9770e', 'other:intron': '#5dade2', 'other:exon_CDS': '#1f5fa8', 'other:ncRNA': '#7d3c98', 'other:SSR': '#58d68d', 'unannotated': '#bdbdbd'}
    bottom = np.zeros(len(piv))
    for c in [c for c in cols if c in piv]: ax.bar(np.arange(len(piv)), piv[c], bottom=bottom, color=cols[c], label=c, width=0.85); bottom += piv[c].values
    ax.set_xticks(np.arange(len(piv))); ax.set_xticklabels(piv.index, rotation=90, fontsize=5.5); ax.set_ylabel('Mb of interruptions'); tidy(ax, f'{SP[G]}: class cut / intrinsic cut / max gap')
axs[0].legend(frameon=False, fontsize=6.5)
banner(fig, '6', 'Step 59e, interruptions: what sits between two consecutive runs of the same family is mostly inserted TEs (' + '; '.join(it_txt) + ' over the 27 settings); '
                 'the rest is intron or unannotated sequence, and the wider the allowed gap, the more non-TE sequence is swept in')
save(fig, 'interruptions')
# ---- order 2
o = rd(f'{V}/cat53/o2/53b_o2_sweep.tsv')
if o is not None:
    fig, axs = plt.subplots(1, len(TAGS), figsize=(5 * len(TAGS), 3.2), sharey=True)
    for ax, G in zip(np.atleast_1d(axs), TAGS):
        g = o[(o.genome == G) & (o.excess_floor_bits == 0.05)]
        for q, c in ((0.05, '#2a78d6'), (0.01, '#0d366b')):
            h = g[g.q_cut == q]; ax.plot(h.min_transitions_variant_only, 100 * h.n_hits / h.n_reported, 'o-', color=c, label=f'q < {q}')
            for r in h.itertuples(): ax.text(r.min_transitions_variant_only, 100 * r.n_hits / r.n_reported + 1.5, f'{r.n_hits}/{r.n_reported}', fontsize=6, ha='center', color=c)
        ax.set_xscale('log'); ax.margins(y=0.18); ax.set_xlabel('minimum variant-only transitions per run'); ax.set_ylabel('% of reported runs with order-2 structure'); ax.legend(frameon=False, fontsize=7); tidy(ax, f'{SP[G]} (excess >= 0.05 bits)')
    oo = o[(o.excess_floor_bits == 0.05) & o.genome.isin(TAGS)]; osh = oo.n_hits / oo.n_reported
    banner(fig, '8', f'Step 53b, order 2: a run has second-order structure when the next variant depends on the current one; this holds for a minority, {rng(osh)} % of the reported runs, '
                     'whichever minimum chain length or q cut is used (BH within each reported set; excess >= 0.05 bits)')
    save(fig, 'order2')
# ---- TE layer
fig, ax = plt.subplots(figsize=(6.5, 3.4)); x = 0; ticks = []; te_v = {'cat59_minrule': [], 'cat59': []}
for G in TAGS:
    for rule, folder, c in (('v4.6 rule (>= 50 % of min(copy, feature))', 'cat59_minrule', '#b7d3f6'), ('feature covers >= 50 % of the copy', 'cat59', '#0d366b')):
        te = rd(f'{V}/{folder}/59c_te_{G}.tsv')
        if te is None: continue
        tb = te.bp.sum(); vals = [100 * te[f'TE_insert_bp_{cut}'].sum() / tb for cut in ('0.2', '0.3', '0.5')]; te_v[folder] += vals
        ax.bar([x, x + 1, x + 2], vals, color=c, width=0.9, label=rule if G == TAGS[0] else None); ticks += [(x + 1, f'{SP[G]}\n{rule.split(" (")[0][:14]}')]; x += 4
ax.set_xticks([t[0] for t in ticks]); ax.set_xticklabels([t[1] for t in ticks], fontsize=6.5); ax.set_ylabel('TE insertions, % of run bp'); ax.legend(frameon=False, fontsize=7); tidy(ax, 'three bars per group = intrinsic cut 0.2 / 0.3 / 0.5')
f1 = lambda v: f'{v:.1f}'
banner(fig, '5', f'Step 59c, TE layer: how a TE label is carried into a run matters more than the intrinsic cut; TE insertions are {rng(te_v["cat59_minrule"], f1)} % of the run sequence under the v4.6 rule '
                 f'but {rng(te_v["cat59"], f1)} % when the feature must cover >= 50 % of the copy (both species, cuts 0.2-0.5)')
save(fig, 'te_layer')
