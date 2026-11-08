"""62_hor_calls_<H>.tsv carry the v4.5 copy-table family label (retired CgaSat/CmaSat names). Rename that column to
legacy_name_do_not_use and add `family` (published id of the array, from arrays.tsv) as column 5. Run once; LF."""
import os, sys
D = r'F:\data\satellite_FISH\satchain_v4.7_2026-11-01\final\hor_narrow_2026-11-09'
fam = {}
for l in open(os.path.join(D, 'arrays.tsv'), encoding='utf-8').read().rstrip('\n').split('\n')[1:]:
    f = l.split('\t'); name = f[1].split('_', 3)[3]
    for r in f[7].split(','): fam[r] = name
for H in ('Gar', 'Mac'):
    p = os.path.join(D, 'o62', f'62_hor_calls_{H}.tsv'); L = open(p, encoding='utf-8', newline='').read().rstrip('\n').split('\n')
    h = L[0].split('\t')
    if 'legacy_name_do_not_use' in h: sys.exit(f'{p}: already fixed')
    i = h.index('family'); h[i] = 'legacy_name_do_not_use'; h.insert(i, 'family'); out = ['\t'.join(h)]
    for l in L[1:]:
        f = l.split('\t'); f.insert(i, fam[f[0]]); out.append('\t'.join(f))
    open(p, 'w', encoding='utf-8', newline='').write('\n'.join(out) + '\n'); print(p, len(out) - 1, 'rows')
