"""palette_v47 (2026-11-01): the v4.6 family palette with the lineage lookup replaced by the SatChain v4.7 sequence clusters.
Colour groups: cluster C2.1 (the centromeric families; GarSat03 / MacSat03 / GarSat09 ... ) keeps the red / orange range that 'lineage A'
had, cluster C2.2 (GarSat04 alone) the blue range of 'lineage B', families without a cluster the green / teal / magenta range. The named
family colours (FAMCOL) are unchanged so every earlier figure keeps its hues. Reads <run folder>/cat59/59_catalogue_<G>.tsv (seq_cluster)."""
import os, sys, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import palette_v46 as P46
from palette_v46 import *                                      # noqa: F401,F403  (FAMCOL, CLASS, TE, SSR, cls, dark, shade, ramp, RDNA, HATCHES ...)
RUN = os.environ.get('SATCHAIN_RUN', 'F:/data/satellite_FISH/satchain_v4.7_2026-11-01')
GROUP_OF = {'C2.1': 'A', 'C2.2': 'B'}                          # hue range only (A = reds / oranges, B = blues); the labels say the cluster, not a lineage
GROUP_LABEL = {'A': 'cluster C2.1 (centromeric families)', 'B': 'cluster C2.2 (GarSat04, 5S-derived)', 'none': 'no cluster assigned'}
LINEAGE_COL = P46.LINEAGE_COL
_CL = None


def cluster_of(f):
    global _CL
    if _CL is None:
        _CL = {}
        for sp in ('Gar', 'Mac'):
            p = f'{RUN}/cat59/59_catalogue_{sp}.tsv'
            if not os.path.exists(p): continue
            t = pd.read_csv(p, sep='\t', usecols=['family', 'family_v46', 'seq_cluster', 'element_class']).drop_duplicates('family_v46')
            for k, v, e in zip(t.family_v46, t.seq_cluster, t.element_class):
                if e == 'satellite' and (v.startswith('C') or k not in _CL): _CL[k] = v
    return _CL.get(str(f).split('_in')[0], 'none')


def lineage_of(f):                                              # keeps palette_v46._build_table() working: it groups by 'A' / 'B' / 'none'
    return GROUP_OF.get(cluster_of(f), 'none')


P46.lineage_of = lineage_of                                     # the table builder inside palette_v46 calls the module-level name
P46._TAB = None


def family(f): return P46.family(f)


def family_hatch(f): return P46.family_hatch(f)
