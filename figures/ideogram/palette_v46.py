"""One colour scheme for every satellite / FISH figure (2026-10-27). Import: sys.path.insert(0, ".../runs_of_record_v4.6_2026-10-27/code"); import palette_v46 as P
Family colours = the slide palette (slidekit.FAMCOL, unchanged): lineage A (centromeric) blues / teals, lineage B (5S-derived) violets / plum, 5S genes gold.
Any other satellite family = neutral grey-blue (same rule as the slides). Repeat classes (rlvl1) and positions have their own hues that no family uses:
  micro = green, mini = orange-brown, gene array = gold, unresolved = grey; position C / I / T = ink / mid grey / light grey.
Shade rule: sub-levels of one hue are mixes with white (shade(c, 0..1)), so a level-2 / level-3 colour is always recognisable as its parent."""
import pandas as pd, matplotlib.colors as mc
# Family colours = centromere_composition.png (user: "that is excellent ... more serious and professional"): lineage A = reds / oranges, lineage B = blues, same hue in both species.
FAMCOL = {"MacSat03": "#9e1b32", "GarSat03": "#c0392b", "GarSat09": "#e67e22", "MacSat162": "#d9822b", "GarSat91": "#8e5a2b",
          "GarSat04": "#3a7bbf", "MacSat06": "#2f5d9e", "MacSat1272": "#7a4fa0"}
LINEAGE_COL = {"A": "#e8a0a0", "B": "#9cbbe0"}                 # any other family of the lineage ("other lineage A / B" in the reference figure)
OTHER_SAT = "#8a9a5b"                                         # family without lineage
TE = {"LINE": "#e6a817", "LTR": "#2e86c1", "TIR": "#7d3c98", "Helitron": "#1e8449"}
CLASS = {"satellite": "#0b7285", "minisatellite": "#6baed6", "microsatellite": "#9ecae1", "gene array": "#e6a100", "unresolved": "#a8a29e"}   # = the thesis genome-composition colours (Satellite DNA teal, mini / micro blues, tRNA / 5S gold)
CLASS_ALIAS = {"mini": "minisatellite", "micro": "microsatellite", "gene_array": "gene array"}
# microsatellite MOTIFS each get their own colour (user 2026-10-27: "we still lack segregating colours for micro and mini sat"). The colours carry NO meaning about repeat length: six hues spread over the wheel
# (blue, cyan, green, lime, violet, magenta; none in the red / orange / gold range that lineage A families and gene arrays use), muted, chosen so that the smallest pairwise CIELAB distance is 31 (a hand-grouped
# blue-shade set had 3 near-identical blues). CLASS["microsatellite"] stays the class colour for places where the motifs are summed.
SSR = {"AT": "#3c60aa", "AC": "#55c3c3", "AG": "#328f59", "AAT": "#98bf4a", "AGAT": "#9460c7", "AAAT": "#b5408e"}


def ssr(motif):
    """colour of one microsatellite motif; accepts '(AC)n', 'AC', 'dinuc_AC', 'SSR (AC)n (kb per Mb)'"""
    import re
    m = re.search(r"\(([ACGT]+)\)n", str(motif)) or re.search(r"(?:^|_)([ACGT]{2,4})$", str(motif)) or re.search(r"^([ACGT]{2,4})$", str(motif))
    return SSR.get(m.group(1), CLASS["microsatellite"]) if m else CLASS["microsatellite"]


_NCLS = None


def class_of(f):
    """repeat class of a family by unit length (micro <= 6 bp, mini 7-60 bp, satellite > 60 bp; nomenclature table): 'satellite' / 'minisatellite' / 'microsatellite' / 'gene array' / 'unresolved'"""
    global _NCLS
    if _NCLS is None:
        t = pd.read_csv("F:/data/satellite_FISH/runs_of_record_v4.6_2026-10-27/satellite_family_nomenclature_v46.tsv", sep=chr(9))
        _NCLS = {k: CLASS_ALIAS.get(v, v) for k, v in zip(t.pipeline_family_id, t.unit_class)}
    return _NCLS.get(str(f).split("_in")[0], "unresolved")


def class_colour(f): return CLASS.get(class_of(f), CLASS["unresolved"])
POSITION = {"centromeric": "#1b1b1b", "interstitial": "#7a7a7a", "terminal": "#bdbdbd"}
_LIN = None


def lineage_of(f):
    global _LIN
    if _LIN is None:
        _LIN = {}
        for sp in ("Gar", "Mac"):
            t = pd.read_csv(f"F:/data/fClaHyb_func_annotation_v1.0/fClaHyb_{sp}/fClaHyb_{sp}_satelliteDNA/fClaHyb_{sp}.satDNA_v46.catalogue.tsv", sep="	", usecols=["family_v46", "lineage"]).drop_duplicates("family_v46")
            for k, v in zip(t.family_v46, t.lineage):
                if v in ("A", "B") or k not in _LIN: _LIN[k] = v
    return _LIN.get(f, "none")


_TAB = None


def _build_table():
    """Every satellite family of the FISH package gets its OWN colour (user: 'each family has the same colour: that is wrong'). Reference families keep the approved colours (centromere_composition scheme).
    The others get distinct hues inside the hue range of their lineage (lineage A: reds / oranges / pinks / browns; lineage B: blues / violets; no lineage: greens, teals, magentas, olive, slate),
    spaced evenly and alternating in lightness, so neighbours in the list never look alike; ordered by satellite Mb so the large families get the most separated colours.
    Where two colours still look close, figures add a hatch glyph (family_hatch) or a row number."""
    global _TAB
    import colorsys
    t = dict(FAMCOL); fam = []
    try:
        fs = pd.read_csv("F:/data/satellite_FISH/runs_of_record_v4.6_2026-10-27/FISH_package/family_summary.tsv", sep=chr(9)).sort_values("sat_Mb", ascending=False); fam = [f for f in fs.family if f not in t and "_in" not in f]
    except Exception: pass
    groups = {"A": [], "B": [], "none": []}
    for f in fam: groups[lineage_of(f) if lineage_of(f) in ("A", "B") else "none"].append(f)
    HR = {"A": [(0, 50), (335, 360)], "B": [(190, 275)], "none": [(75, 175), (280, 335), (50, 75)]}                  # hue ranges in degrees
    import numpy as np
    placed = [np.array(mc.to_rgb(c)) for c in t.values()]
    for g, fs_ in groups.items():                                                  # greedy farthest-point choice inside the lineage's hue range: every new colour is as far as possible from all colours already used
        if not fs_: continue
        cand = []
        for lo, hi in HR[g]:
            for h in range(lo, hi, 3):
                for light in (0.32, 0.40, 0.48, 0.56, 0.64):
                    for sat in (0.40, 0.55, 0.70):
                        r, g_, b = colorsys.hls_to_rgb(h / 360, light, sat); cand.append(np.array([r, g_, b]))
        cand = np.array(cand)
        for f in fs_:
            dmin = np.min(np.linalg.norm(cand[:, None, :] - np.array(placed)[None, :, :], axis=2), axis=1); k = int(np.argmax(dmin)); c = cand[k]
            t[f] = "#%02x%02x%02x" % tuple(int(round(x * 255)) for x in c); placed.append(c)
    _TAB = t; return t


HATCHES = ["", "///", "xxx", "...", "\\\\\\", "+++"]


def family(f):
    f = str(f).split("_in")[0]
    if f in FAMCOL: return FAMCOL[f]
    tab = _TAB or _build_table()
    if f in tab: return tab[f]
    h = sum(ord(c) * (i + 1) for i, c in enumerate(f)) % 360                       # family outside the package: stable hash colour
    import colorsys; r, g_, b = colorsys.hls_to_rgb(h / 360, 0.5, 0.6); return "#%02x%02x%02x" % (int(r * 255), int(g_ * 255), int(b * 255))


def family_hatch(f):
    """glyph that goes with the colour when several families are on one figure: reference families none; the others cycle through HATCHES by their rank"""
    f = str(f).split("_in")[0]
    tab = _TAB or _build_table()
    if f in FAMCOL or f not in tab: return ""
    order = [k for k in tab if k not in FAMCOL]
    return HATCHES[(order.index(f) // 12) % len(HATCHES)]
def cls(k): return CLASS.get(CLASS_ALIAS.get(k, k), CLASS["unresolved"])
def dark(c, k=0.7):
    r = mc.to_rgb(c); return tuple(k * x for x in r)
def shade(c, t):
    """t = 0 (the colour itself) .. 1 (white)"""
    r = mc.to_rgb(c); return tuple(x + (1 - x) * t for x in r)
def ramp(c, v):
    """0 = white, 1 = the colour (v clipped to 0..1)"""
    v = min(1.0, max(0.0, v)); return shade(c, 1 - v)
INK, INK2, GRID, CHROM_FILL, CHROM_EDGE = "#1b1b1b", "#5b5b5b", "#d9d6cf", "#eceae4", "#9a9890"
RDNA = {"5S rDNA clusters": "#e6a817", "45S rDNA": "#6c3483", "interstitial telomere motifs": "#1e8449", "telomere": "#1e8449"}
# centromere = pale warm grey block, never black (user 2026-10-27: black centromeres and black outlines hide the focal arrays); array outlines = dark(family colour), never black
CEN, CEN_EDGE, CHROM_LINE = "#c9c3b4", "#7a7467", "#dcd8ce"
