"""satpaths: one resolver for the annotation folder layout (2026-11-01, audit CB16 / C17 / D08 / CA25).
The annotation folders were renamed in 2026-10 from `<FUNC>/_repeats_EDTA/` to `<FUNC>/fClaHyb_<G>_repeats_EDTA/`; several
scripts still hard-coded the old names and either crashed or silently calibrated on an unmasked null. Every script now calls
func_path(FUNC, G, sub, name): it tries the new layout first, then the old one, and raises (never returns a guess) when the file is
missing. overlay_path() does the same for the annotation_overlays/<G> folder. usage in a script:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from satpaths import func_path
    gff = func_path(FUNC, G, '_repeats_EDTA', f'fClaHyb_{G}.TEanno.sorted.gff3.gz')"""
import os


def func_path(func, g, sub, name, required=True):
    """<func>/fClaHyb_<g><sub>/<name> (current layout) or <func>/<sub>/<name> (pre-2026-10 layout)."""
    for d in (os.path.join(func, f'fClaHyb_{g}{sub}'), os.path.join(func, sub)):
        p = os.path.join(d, name)
        if os.path.exists(p): return p
    if required: raise FileNotFoundError(f'{sub}/{name} not found under {func} (tried fClaHyb_{g}{sub} and {sub})')
    return None


def func_dir(func, g, sub):
    for d in (os.path.join(func, f'fClaHyb_{g}{sub}'), os.path.join(func, sub)):
        if os.path.isdir(d): return d
    raise FileNotFoundError(f'folder {sub} not found under {func}')


def overlay_path(ovl, g, name, required=True):
    """<ovl>/<name>; <ovl> is annotation_overlays/<G> (or the parent, then <ovl>/<G>/<name>)."""
    for p in (os.path.join(ovl, name), os.path.join(ovl, g, name)):
        if os.path.exists(p): return p
    if required: raise FileNotFoundError(f'{name} not found in {ovl}')
    return None


def genome_fasta(func, g): return func_path(func, g, '_genome', f'fClaHyb_{g}_LG.fa')


def gaps_bed(func, g): return func_path(func, g, '_genome', f'fClaHyb_{g}_LG.gaps.detgaps.bed')
