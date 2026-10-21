"""Summarise the 53b second-order test (H1 - H2 vs simulations of each array's own first-order chain) per array, sweeping the z cut."""
import sys, pandas as pd, numpy as np
S = sys.argv[1]; rows = []
for G in ('Gar', 'Mac'):
    cat = pd.read_csv(f'{S}/cat59/59_catalogue_{G}.tsv', sep='\t', usecols=['array_id', 'name_v46', 'family_v46', 'element_class', 'unit_bp'], low_memory=False)
    fam = cat.groupby('array_id').agg(family=('family_v46', 'first'), name=('name_v46', 'first'), cls=('element_class', 'first'), unit_bp=('unit_bp', 'median')).reset_index()
    for src in ('existing', 'new'):
        d = pd.read_csv(f'{S}/cat53/{G}_{src}/53b/53b_topologies.tsv', sep='\t', low_memory=False)
        a = d[d.unit == 'array'].copy(); a['G'] = G; a['src'] = src; a['bp'] = a.end - a.start0
        rows.append(a.merge(fam, on='array_id', how='left'))
A = pd.concat(rows); A['order2_z'] = pd.to_numeric(A.order2_z, errors='coerce'); A['order2_excess_bits'] = pd.to_numeric(A.order2_excess_bits, errors='coerce')
A['memory_z'] = pd.to_numeric(A.memory_z, errors='coerce')
A.to_csv(f'{S}/53b_o2_arrays.tsv', sep='\t', index=False)
for G, g in A.groupby('G'):
    t = g[g.order2_z.notna()]
    print(f'{G}: {len(g)} arrays, {len(t)} testable (order-2 z defined), testable bp {t.bp.sum()/1e6:.2f} of {g.bp.sum()/1e6:.2f} Mb')
    for zc in (2, 3, 5, 10):
        for ex in (0.0, 0.05, 0.1):
            h = t[(t.order2_z >= zc) & (t.order2_excess_bits >= ex)]
            print(f'   z>={zc:<2} excess>={ex:<4}: {len(h):5d} arrays ({len(h)/max(len(t),1)*100:4.1f}%), {h.bp.sum()/1e6:6.2f} Mb ({h.bp.sum()/max(t.bp.sum(),1)*100:4.1f}% of testable bp)')
    h = t[t.order2_z >= 3]; m = t[t.memory_z >= 3]
    print(f'   overlap with first-order memory (z>=3): order2 hits {len(h)}, of which with memory {len(h[h.memory_z >= 3])}; memory hits {len(m)}')
    print('   by size_flag:', t.groupby('size_flag').apply(lambda x: f'{(x.order2_z >= 3).sum()}/{len(x)}').to_dict())
    print('   top families by order-2 bp (z>=3):'); f = h.groupby('family').agg(n=('array_id', 'size'), bp=('bp', 'sum'), z=('order2_z', 'median'), ex=('order2_excess_bits', 'median')).sort_values('bp', ascending=False).head(8)
    print(f.to_string(float_format=lambda v: f'{v:.3g}'))
    print('   top contexts:', h.order2_top_context.value_counts().head(6).to_dict())
