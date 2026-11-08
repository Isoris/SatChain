"""Shared house style of the satellite / RTE paper figures (2026-10-27; re-pointed to SatChain v4.7 on 2026-11-01).
SOURCE OF THE RULES: the pop-gen figures' house_style() in F:\\data\\fClaHyb_population_v1.0\\_heterozygosity_roh_ngsFHMM\\code\\build_ROH_per_sample_216.py (reference figure
Figure_ROH_per_sample_216_ancestry.png; collage F:\\data\\fClaHyb_population_v1.0\\_population_structure_pcangsd_ngsadmix\\figures_v2_2026-10-27\\Figure_popstructure_collage_216):
Arial 7 pt, navy #17252D text, 0.5 pt hairlines, no top / right spines, bold panel letters, ROH colours #F1E492 #EBDB87 #95BADD #71A29C #283B7C #83251D, plus terracotta #D08C6E, lavender #B8A9CC, grey #A7AEB3.
Every panel is a module  fig<N>_panel_<x>.py  with   draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None) -> dict(numbers)   drawing INTO a shared matplotlib figure at a position given in millimetres
from the top-left of the figure (this is how the pop-gen collage is built: one matplotlib figure, data drawn directly, not pasted pictures).  The assembly scripts call house_style() once and then every draw().
usage:  import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "code")); import satellite_style as ST; ST.house_style()

TERMINOLOGY (SatChain v4.7, 2026-11-01): the word "lineage" is dropped.  The k-mer clustering of centromeric monomers (step 55) is stable only at K = 2: C2.2 = GarSat04 alone,
C2.1 = every other centromeric family (GarSat03, GarSat09, GarSat10, MacSat03, MacSat162 ...).  The two colour groups of the figures are therefore
  CEN    = the GarSat03 / MacSat03 group, "the RTE-carrying centromeric families" (seq_cluster C2.1): reds
  FIVE_S = GarSat04 (5S-derived satellite, seq_cluster C2.2) and, where a panel needs one colour for the 5S-related arrays, the 5S rDNA gene arrays GarSat49 / MacSat06: blues.
MacSat06 is the C. macrocephalus 5S rDNA array (element_class gene_array:5S_rRNA), a gene array, never a satellite family group."""
import os, sys, pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

MM = 1 / 25.4
NAVY = "#17252D"; GRID = "#E5E7EB"; MUTED = "#8A949C"; RED = "#b91c1c"
ROH = ["#F1E492", "#EBDB87", "#95BADD", "#71A29C", "#283B7C", "#83251D"]
TERRACOTTA = "#D08C6E"; LAVENDER = "#B8A9CC"; GREY = "#A7AEB3"
CEN = "#83251D"; FIVE_S = "#283B7C"                       # the two groups (see the module docstring); same hues as palette_v46 / centromere_composition, tuned to the house palette
CEN_LIGHT = TERRACOTTA; FIVE_S_LIGHT = "#95BADD"
GROUP = {"cen": CEN, "5S": FIVE_S}; GROUP_LIGHT = {"cen": CEN_LIGHT, "5S": FIVE_S_LIGHT}
GROUP_NAME = {"cen": "GarSat03 / MacSat03 group", "5S": "GarSat04 (5S-derived)"}
SPECIES = {"Mac": "#4a9a5f", "Gar": "#9b5f9b"}          # species colours of the RTE copy tree (chosen away from the group reds / blues)
RTE = "#e6a817"                           # the RTE LINE (gold, same as palette_v46 TE LINE); the user wants it called RTE in every figure, caption and file label, never the old name
CENTROMERE = "#c9c3b4"; CEN_EDGE = "#7a7467"            # pale warm grey (never black), as palette_v46
NUC = {"A": "#2E8B57", "C": "#3B5BDB", "G": "#D08C6E", "T": "#C0392B"}      # nucleotide view: A / C / G / T (muted versions of the IGV scheme used on thesis slide G07c)
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "ideogram"))    # F:/data/SatChain/figures/ideogram/palette_v46.py (family colours; identical to runs_of_record_v4.6_2026-10-27/code/palette_v46.py)

# ---- SatChain v4.7 run of record (2026-11-01): the tables every catalogue-derived panel reads
V47 = pathlib.Path("F:/data/satellite_FISH/satchain_v4.7_2026-11-01/final")   # 2026-11-02: final run (clean candidate set, nothing carried)
RUNS = V47 / "runs_of_record"                                              # {Gar,Mac}/runs_of_record.tsv, family_totals.tsv, chrom_totals.tsv, centromere_check.tsv
NOM = RUNS / "satellite_family_nomenclature.tsv"                           # haplotype, family, name, element_class, unit_bp_median, bp, n_runs, seq_cluster
CAT59 = V47 / "cat59"                                                      # 59_catalogue_{Gar,Mac}.tsv (column lineage -> seq_cluster; family_v46, name_v46, element_class, tier)
RUNS_VERSION = "SatChain v4.7 runs of record (satchain_v4.7_2026-11-01)"
# chromosomes with NO row in the v4.7 runs of record / catalogue (chrom_totals: 0 runs); a data gap of the v4.7 run, reported in every caption that depends on it
MISSING_CHROMS = {"Mac": [], "Gar": []}                                    # 2026-11-01 20:20: the six chromosomes were restored (val56 track bug fixed); kept empty for the "no rows" fallback


def runs(sp, **kw):
    import pandas as pd
    return pd.read_csv(RUNS / sp / "runs_of_record.tsv", sep="\t", **kw)


def catalogue(sp, **kw):
    import pandas as pd
    return pd.read_csv(CAT59 / f"59_catalogue_{sp}.tsv", sep="\t", **kw)


def nomenclature():
    """the v4.7 nomenclature table (both haplotypes), with a column `is_gene_array`"""
    import pandas as pd
    global _NOM
    if "_NOM" not in globals():
        _NOM = pd.read_csv(NOM, sep="\t"); _NOM["is_gene_array"] = _NOM.element_class.astype(str).str.startswith("gene_array")
    return _NOM


def family_cluster(f):
    """seq_cluster of a family in v4.7 (C2.1 / C2.2 / none / na), '' if unknown"""
    n = nomenclature(); x = n[n.family == f]
    return str(x.seq_cluster.iloc[0]) if len(x) else ""


def family_group(f):
    """'cen' (seq_cluster C2.1), '5S' (GarSat04 = C2.2, or a 5S rDNA gene array), 'gene' (other gene array) or 'other'"""
    n = nomenclature(); x = n[n.family == f]
    if not len(x): return "other"
    ec = str(x.element_class.iloc[0]); c = str(x.seq_cluster.iloc[0])
    if ec.startswith("gene_array:5S"): return "5S"
    if ec.startswith("gene_array"): return "gene"
    return "cen" if c == "C2.1" else "5S" if c == "C2.2" else "other"


def is_gene_array(f):
    n = nomenclature(); x = n[n.family == f]
    return bool(x.is_gene_array.iloc[0]) if len(x) else str(f).startswith("gene_array")



MUTE = float(os.environ.get("SAT_MUTE", "0.3"))      # 2026-11-01 user choice: 30 % muted + bevel on every block panel        # 0 = house colours as they are; 0.25 = 25 % towards the ink-and-paper tone (less saturated, slightly darker)


def tone(fc):
    """family colour pulled MUTE of the way towards a muted print tone (saturation down, value slightly down); MUTE = 0 returns fc"""
    import colorsys
    from matplotlib.colors import to_rgb, to_hex
    if MUTE <= 0: return fc
    h, l, sat = colorsys.rgb_to_hls(*to_rgb(fc)); sat = sat * (1 - 0.8 * MUTE); l = l * (1 - 0.25 * MUTE)
    return to_hex(colorsys.hls_to_rgb(h, l, sat))


def block(ax, x, y, w, h, fc, zorder=2, bevel=0.16):
    """A satellite-run block with a subtle bevel (2026-11-01, user: flat blocks "look cartoon"): the flat colour, a translucent white highlight
    over the upper part, a translucent dark band along the lower edge and a hairline outline in a darker shade. Vector-safe (three patches)."""
    from matplotlib.patches import Rectangle
    from matplotlib.colors import to_rgb
    fc = tone(fc); r, g, b = to_rgb(fc); dark = (r * 0.62, g * 0.62, b * 0.62)   # the block panels (1G) are the only muted ones
    ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor=dark, lw=0.25, zorder=zorder))
    ax.add_patch(Rectangle((x, y + h * 0.55), w, h * 0.45, facecolor="white", edgecolor="none", lw=0, alpha=bevel, zorder=zorder + 0.01))
    ax.add_patch(Rectangle((x, y), w, h * 0.18, facecolor="black", edgecolor="none", lw=0, alpha=bevel * 0.7, zorder=zorder + 0.01))


def family_colour(f):
    """each satellite family keeps its OWN muted colour (palette_v46.family: reference families fixed, others unique)"""
    import palette_v46 as P
    return P.family(f)                                        # house colours everywhere (user 2026-11-01: "before was more colorful, it was better"); only block() mutes


def house_style():
    av = {f.name for f in fm.fontManager.ttflist}
    face = "Arial" if "Arial" in av else "TeX Gyre Heros" if "TeX Gyre Heros" in av else "Liberation Sans" if "Liberation Sans" in av else "DejaVu Sans"
    matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42, "font.family": face, "font.size": 7,
                                "text.color": NAVY, "axes.labelcolor": NAVY, "axes.edgecolor": NAVY, "xtick.color": NAVY, "ytick.color": NAVY,
                                "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2.4, "ytick.major.size": 2.4,
                                "axes.grid": False, "grid.color": GRID, "grid.linewidth": 0.4, "axes.spines.top": False, "axes.spines.right": False,
                                "legend.frameon": False, "legend.handlelength": 1.3, "figure.facecolor": "white", "savefig.facecolor": "white",
                                "axes.titlelocation": "left", "axes.titleweight": "bold",
                                "mathtext.fontset": "custom", "mathtext.rm": face, "mathtext.it": face + ":italic", "mathtext.bf": face + ":bold"})
    return face


def new_figure(w_mm=183.0, h_mm=200.0):
    return plt.figure(figsize=(w_mm * MM, h_mm * MM), dpi=500)


def ax_mm(fig, left, top, w, h, **kw):
    """axes at (left, top) with size (w, h), all in mm measured from the TOP-LEFT of the figure"""
    W, H = fig.get_size_inches() / MM
    return fig.add_axes([left / W, 1 - (top + h) / H, w / W, h / H], **kw)


def letter(fig, left, top, s):
    """bold CAPITAL panel letter, 8 pt, at (left, top) mm"""
    W, H = fig.get_size_inches() / MM
    fig.text(left / W, 1 - top / H, s.upper(), fontsize=8, fontweight="bold", ha="left", va="top", color=NAVY)


def hairline_axes(ax):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_linewidth(0.5)
    return ax


def save(fig, base, dpi=500):
    """PNG at 500 dpi + vector PDF + SVG (fonts as text)"""
    fig.savefig(base + ".png", dpi=dpi); fig.savefig(base + ".pdf"); fig.savefig(base + ".svg")


NOM_V46 = pathlib.Path("F:/data/satellite_FISH/runs_of_record_v4.6_2026-10-27/satellite_family_nomenclature_v46.tsv")   # ONLY to translate the v4.6 manuscript codes of tables built before v4.7


def pipeline_name(name):
    """ONE family name everywhere (decided 2026-10-28): the pipeline name (MacSat03, GarSat04 ...). Manuscript-style codes are mapped back through the v4.7 nomenclature (name -> family)
    and the v4.7 catalogues (name_v46 -> family).  The composition tables of centromere_composition_2026-10-24 carry the v4.6 manuscript codes, which v4.7 renumbered (e.g. CmaSat002-214 =
    MacSat06 became 'gene_array:5S_rRNA|MacSat06', CmaSat011-215 = MacSat162 became CmaSat012-215), so the v4.6 nomenclature is read as the LAST fallback for those codes only
    (no v4.6 number enters any figure).  Pipeline names and unknown names pass through unchanged."""
    import pandas as pd
    global _PN
    if "_PN" not in globals():
        _PN = {}
        if NOM_V46.exists():
            o = pd.read_csv(NOM_V46, sep="\t"); _PN.update(dict(zip(o.name_v46.astype(str), o.pipeline_family_id))); _PN.update(dict(zip(o.manuscript_name.astype(str), o.pipeline_family_id)))
        for sp in ("Mac", "Gar"):
            c = catalogue(sp, usecols=["family", "name_v46"]).drop_duplicates("name_v46")
            _PN.update(dict(zip(c.name_v46.astype(str), c.family)))
        n = nomenclature(); _PN.update(dict(zip(n.name.astype(str), n.family)))
        for k in list(_PN):
            if "|" in k: _PN[k.split("|")[-1]] = _PN[k]
    return _PN.get(str(name), str(name))
