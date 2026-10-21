#!/usr/bin/env python3
"""Topology views of one array (step 53): (1) per-segment transition graphs before / after the class merge, major and minor
edges; (2) the transition matrix as a band along the array (one line per matrix entry, sliding copy windows);
(3) a filmstrip of local transition graphs, frame by frame along the array.
States = step-03 classes of confident copies (identity >= 80, complete), merged by 53_class_merges.tsv; other copies = U.
Path breaks where consecutive copies are separated by > 0.5 x period (no edge across a break).
usage: plot53_topology.py <03_monomer_copies.tsv> <53 run dir> <array_id> <segA,segB> <png> [frames]"""
import sys, csv, collections, math
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, Arc, Rectangle

COP, RUN, AID, SEGS, PNG = sys.argv[1:6]; NF = int(sys.argv[6]) if len(sys.argv) > 6 else 18
csv.field_size_limit(1 << 30)
INK, INK2, GRID = '#1f1f1e', '#5f5e58', '#e4e2dc'
segs = {r['segment']: r for r in csv.DictReader(open(f'{RUN}/53_segments.tsv'), delimiter='\t') if r['array_id'] == AID}
cps = [int(r['position']) for r in csv.DictReader(open(f'{RUN}/53_change_points.tsv'), delimiter='\t') if r['array_id'] == AID and r['status'] == 'accepted']
SEGCOL = {'homogeneous_monomeric': '#2a78d6', 'variant_mixed': '#1baf7a', 'alternating': '#6250d6', 'interleaved_families': '#eb6834', 'unresolved': '#b4b2a9'}
period = float(next(iter(segs.values()))['period_02'])
merge = {}
for r in csv.DictReader(open(f'{RUN}/53_class_merges.tsv'), delimiter='\t'):
    if r['array_id'] == AID and r['decision'] == 'merge': merge.setdefault(r['class_a'], set()).add(r['class_b'])
parent = {}
def find(x):
    parent.setdefault(x, x)
    while parent[x] != x: x = parent[x]
    return x
for a_, bs in merge.items():
    for b_ in bs: parent[find(b_)] = find(a_)
cop = []
with open(COP, newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if r['array_id'] != AID: continue
        try: idt = float(r['best_alignment_identity'])
        except ValueError: idt = 0
        conf = r['complete_or_truncated'] == 'complete' and idt >= 80
        cop.append((int(r['genomic_start_0based']), int(r['genomic_end']), r['monomer_class'] if conf else 'U'))
cop.sort()
allc = sorted({c for *_, c in cop if c != 'U'})
groups = collections.defaultdict(list)
for c in allc: groups[find(c)].append(c)
merged = {c: '+'.join(sorted(groups[find(c)])) for c in allc}; merged['U'] = 'U'
raw_states = allc + ['U']; mer_states = sorted(set(merged.values()) - {'U'}) + ['U']
PAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4']
col_raw = {c: PAL[i % 5] for i, c in enumerate(allc)}; col_raw['U'] = '#b4b2a9'
col_mer = {s: col_raw[s.split('+')[0]] if s != 'U' else '#b4b2a9' for s in mer_states}


def transitions(lo, hi, lab):
    t = collections.Counter(); n = collections.Counter(); prev = None; pe = None
    for s, e, c in cop:
        if s < lo or e > hi: continue
        st = lab[c] if c in lab else c
        n[st] += 1
        if prev is not None and s - pe <= 0.5 * period: t[(prev, st)] += 1
        prev, pe = st, e
    return t, n


def draw_graph(ax, t, n, states, cols, title, small=False):
    ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.6, 1.6 if small else 2.15); ax.set_aspect('equal'); ax.axis('off')
    K = len(states); pos = {s: (math.cos(math.pi / 2 + 2 * math.pi * i / K), math.sin(math.pi / 2 + 2 * math.pi * i / K)) if K > 1 else (0, 0) for i, s in enumerate(states)}
    out = collections.Counter()
    for (a, b), v in t.items(): out[a] += v
    tot = max(1, sum(n.values()))
    for (a, b), v in sorted(t.items(), key=lambda kv: kv[1]):
        p = v / max(1, out[a]); major = p >= 0.15; lw = (0.4 + 5 * p) * (0.6 if small else 1)
        ls = '-' if major else (0, (2, 2)); c = cols[a] if major else INK2
        if a == b:
            x, y = pos[a]; d = 0.42; ang = math.atan2(y, x) if (x or y) else math.pi / 2
            cx, cy = x + d * math.cos(ang), y + d * math.sin(ang)
            ax.add_patch(Arc((cx, cy), 0.55, 0.55, color=c, lw=lw, ls=ls if major else '--', zorder=1))
            if not small: ax.text(cx + 0.38 * math.cos(ang), cy + 0.38 * math.sin(ang), f'{p:.2f}', ha='center', va='center', fontsize=7.5, color=INK2)
        else:
            ax.add_patch(FancyArrowPatch(pos[a], pos[b], connectionstyle='arc3,rad=0.22', arrowstyle='-|>', mutation_scale=8 if small else 12,
                                         color=c, lw=lw, ls=ls, shrinkA=14 if small else 22, shrinkB=14 if small else 22, zorder=1))
            if not small and p >= 0.05:
                mx, my = (pos[a][0] + pos[b][0]) / 2, (pos[a][1] + pos[b][1]) / 2; dx, dy = pos[b][0] - pos[a][0], pos[b][1] - pos[a][1]
                ax.text(mx - 0.18 * dy, my + 0.18 * dx, f'{p:.2f}', fontsize=7.5, color=INK2, ha='center', va='center')
    for s in states:
        r = 0.12 + 0.33 * math.sqrt(n[s] / tot) if n[s] else 0.06
        ax.add_patch(Circle(pos[s], r * (0.8 if small else 1), facecolor=cols[s], edgecolor='white', lw=1.5, zorder=3))
        if not small:
            ang = math.atan2(pos[s][1], pos[s][0]) if K > 1 else -math.pi / 2; d = r + 0.32
            ax.text(pos[s][0] + d * math.cos(ang - 1.25), pos[s][1] + d * math.sin(ang - 1.25), f'{s}\n{n[s]:,} copies', ha='center', va='center',
                    fontsize=8.5, color=INK, bbox=dict(facecolor='white', edgecolor='none', pad=0.6, alpha=0.85), zorder=4)
    ax.set_title(title, fontsize=9 if small else 10.5, color=INK, loc='center')


plt.rcParams.update({'font.family': 'Aptos', 'font.size': 10, 'axes.edgecolor': GRID, 'xtick.color': INK2, 'ytick.color': INK2})
sA, sB = SEGS.split(',')
fig = plt.figure(figsize=(16, 14.5)); gs = fig.add_gridspec(5, NF // 2, height_ratios=[1.9, 0.12, 1.3, 1.15, 1.15], hspace=0.45, wspace=0.05)
# (1) segment graphs
for j, (sg, merged_view) in enumerate(((sA, False), (sA, True), (sB, False), (sB, True))):
    ax = fig.add_subplot(gs[0, j * (NF // 8):(j + 1) * (NF // 8)]); S = segs[sg]; lo, hi = int(S['start0']), int(S['end'])
    lab = merged if merged_view else {c: c for c in raw_states}
    t, n = transitions(lo, hi, lab)
    draw_graph(ax, t, n, mer_states if merged_view else raw_states, col_mer if merged_view else col_raw,
               f"segment {sg} ({S['profile'].replace('_', ' ')})\n{lo / 1e6:.3f}-{hi / 1e6:.3f} Mb, {'after' if merged_view else 'before'} class merge")
    ax.add_patch(Rectangle((-1.6, 2.07), 3.2, 0.08, color=SEGCOL[S['profile']], clip_on=False))
# (2) transition matrix band along the array
ax = fig.add_subplot(gs[2, :]); x0, x1 = cop[0][0], cop[-1][1]
step = max(1, len(cop) // 160); win = 40
keys = [(a, b) for a in mer_states if a != 'U' for b in mer_states if b != 'U']
series = {k: [] for k in keys}; xs = []; ufr = []
for i in range(0, len(cop) - win, step):
    lo, hi = cop[i][0], cop[i + win][1]; t, n = transitions(lo, hi, merged); out = collections.Counter()
    for (a, b), v in t.items(): out[a] += v
    xs.append((lo + hi) / 2e6); ufr.append(n['U'] / max(1, sum(n.values())))
    for k in keys: series[k].append(t[k] / out[k[0]] if out[k[0]] else float('nan'))
LS = {True: '-', False: '--'}
for k in keys:
    ax.plot(xs, series[k], lw=1.6, color=col_mer[k[0]], ls=LS[k[0] == k[1]], label=f'P({k[1]} | {k[0]})')
ax.plot(xs, ufr, lw=1, color='#b4b2a9', label='uncertain-copy fraction')
for c in cps: ax.axvline(c / 1e6, color=INK, lw=1.1)
ax.set_ylim(-0.02, 1.02); ax.set_xlim(x0 / 1e6, x1 / 1e6); ax.set_ylabel('transition\nprobability', color=INK2)
ax.set_title('transition matrix along the array (windows of 40 copies; solid = stay in the same state, dashed = switch)', fontsize=10.5, color=INK, loc='left')
ax.legend(frameon=False, fontsize=8.5, ncol=5, loc='upper center', bbox_to_anchor=(0.5, -0.2)); ax.grid(color=GRID, lw=0.6)
for sp_ in ('top', 'right'): ax.spines[sp_].set_visible(False)
for S in segs.values(): ax.add_patch(Rectangle((int(S['start0']) / 1e6, 1.03), (int(S['end']) - int(S['start0'])) / 1e6, 0.04, color=SEGCOL[S['profile']], clip_on=False))
# (3) filmstrip
edges = [x0 + (x1 - x0) * i / NF for i in range(NF + 1)]
for f in range(NF):
    ax = fig.add_subplot(gs[3 + f // (NF // 2), f % (NF // 2)]); lo, hi = edges[f], edges[f + 1]
    t, n = transitions(lo, hi, merged); mid = (lo + hi) / 2
    prof = next((S['profile'] for S in segs.values() if int(S['start0']) <= mid < int(S['end'])), 'unresolved')
    draw_graph(ax, t, n, mer_states, col_mer, f'{lo / 1e6:.2f}-{hi / 1e6:.2f}', small=True)
    ax.add_patch(Rectangle((-1.6, -1.6), 3.2, 3.2, fill=False, edgecolor=SEGCOL[prof], lw=2.5, clip_on=False))
fig.suptitle(f'{AID}: transition-graph topology, segment {sA} vs segment {sB}, with a frame-by-frame filmstrip along the array',
             x=0.01, ha='left', fontsize=13, fontweight='bold', color=INK)
fig.text(0.01, 0.955, f"states after merge: {', '.join(s for s in mer_states if s != 'U')}; U = uncertain copy (identity < 80 % or truncated). "
         'Solid coloured edge = major (>= 15 % of the source state\'s transitions), dashed grey = minor; width = probability; node area = copies. '
         'Frame border = segment colour.', fontsize=9, color=INK2, wrap=True)
fig.savefig(PNG, dpi=140, facecolor='white', bbox_inches='tight'); print('wrote', PNG)
