# Paper figures: layout, panels, contracts (2026-10-27; re-pointed to SatChain v4.7 on 2026-11-01)

TERMINOLOGY (2026-11-01): the word "lineage" is dropped (SatChain v4.7: the k-mer clustering of centromeric monomers is stable only at K = 2, C2.2 = GarSat04 alone, C2.1 = everything else). The two groups of the figures are the GarSat03 / MacSat03 group (the RTE-carrying centromeric families, colour ST.CEN, reds) and GarSat04 (the 5S-derived satellite, ST.FIVE_S, blues) with the 5S rDNA gene arrays GarSat49 / MacSat06 (gene arrays, never a satellite family group). Catalogue-derived panels read F:/data/satellite_FISH/satchain_v4.7_2026-11-01/{runs_of_record,cat59} through code/satellite_style.py (ST.RUNS, ST.CAT59, ST.NOM); the v4.7 tables have NO rows for Mac LG14 / 17 / 18 and Gar LG23 / 27 / 28 (ST.MISSING_CHROMS), which every affected panel reports.


Style = the pop-gen house style (see code/satellite_style.py: Arial 7 pt, navy #17252D, 0.5 pt hairlines, no top/right spines, bold panel letters, house colours).
Every figure is ONE matplotlib figure 183 mm wide, data drawn directly (not pasted pictures). Each panel is a module `panels/figN_panel_X.py` exposing

    draw(fig, left_mm, top_mm, w_mm, h_mm, letter=None) -> dict   # draws INTO `fig`, mm from the top-left of the figure; returns the numbers it shows
    if __name__ == "__main__":  renders the panel alone at its final size (500 dpi PNG + PDF) into panels/_standalone/ for QA

Inside a panel: no titles bigger than 7 pt bold, no connector lines, no frames around the panel, panel letter = `ST.letter(fig, left_mm, top_mm, letter)` at the panel's top-left corner (the panel content starts ~4 mm below).
Evidence levels, parameters, sources and sweeps go to the CAPTION STUB `panels/figN_panel_X_caption.md` (not into the picture) plus the table `panels/figN_panel_X_numbers.tsv`.

## Figure 1: genomic organisation of the satellite DNA families in the F1 hybrid catfish (183 x 205 mm; assemble_v2 layout 183 x 282 mm with banners)
| panel | left, top, w x h (mm) | content |
|---|---|---|
| a | 0, 0, 92 x 56 | trees of the two groups (5S-related families: GarSat04 + the 5S rDNA arrays GarSat49 / MacSat06; centromeric families: GarSat03 / MacSat03 group), clades S1.. / R1.., tips coloured by family, species markers, 5S rRNA references, scale bar |
| b | 95, 0, 88 x 56 | nucleotide view of monomer alignments: centromeric (GarSat03, MacSat03: same stripes in both species) and 5S-related (GarSat04 5S-derived satellite, MacSat06 5S rDNA array: 5S-related segment bracketed); example run = most homogeneous of the 5 largest v4.7 runs; A/C/G/T muted colours |
| c | 0, 60, 88 x 84 | reference centromeres: composition of every centromere in both haplotypes (ideogram-like: satellite families, strand, RTE / TE insertions, homogenised core) |
| d | 91, 60, 92 x 84 | organisation of the major centromeres to scale in BOTH species (runs, interruptions, RTE fragments, gaps); gene arrays shown as gene arrays, not satellites |
| e | 0, 148, 66 x 52 | the MacSat03 unit carries the RTE: the composite unit (satellite block + following interruption) repeating along the centromeres |
| f | 69, 148, 56 x 52 | 5S rDNA clusters along the chromosomes of both haplotypes, against the 5S-related arrays (GarSat04, GarSat49, MacSat06) |
| g | 128, 148, 55 x 52 | the 5S gene mapped onto the 5S-related monomers (GarSat04, MacSat06, GarSat49; GarSat03 / MacSat03 = controls): where the 80-110 bp 5S-related segment sits and how conserved it is |

## Figure 2: the RTE LINE built into the MacSat03 centromeres (183 x 172 mm; assemble_v2 layout 183 x 198 mm with banners)
| panel | left, top, w x h (mm) | content |
|---|---|---|
| a | 0, 0, 88 x 52 | the element: RTE-clade LINE, 3,117 bp consensus, ORF, APE motif, RT domain, 3' end; where it sits in the unit |
| b | 91, 0, 92 x 52 | what each MacSat03 centromere keeps of the element, and the unit anatomy (type I, type II) |
| c | 0, 56, 56 x 62 | RTE core copies inside / outside the centromere by family group of the centromere (cen = GarSat03 / MacSat03 group, 5S = GarSat04, other; v4.7 intervals via code/rte_core_copies_v47.py) (+ window sweep) |
| d | 59, 56, 66 x 66 | tree of RTE copies, both species |
| e | 128, 56, 55 x 66 | positions of the type II insert copies along the chromosomes of both haplotypes |
| f | 0, 126, 100 x 46 | one interrupted array (Mac LG04) to scale: classes, RTE fragments, spacer types A, B, C |
| g | 103, 126, 38 x 46 | age / divergence of RTE copies (relative, sweep-stable) |
| h | 144, 126, 39 x 46 | homogenisation of the unit within centromeres (shuffle null) |

(Panel letters of Figure 2 are a to h; the old collage letters differ.)

## Rules every panel respects
Evidence level L1 (sequence) unless the picture says otherwise (LG04 spacer reading = L2); every threshold / window / k that a panel depends on is swept and the stability is reported in the caption stub; each satellite family keeps its own muted colour (palette_v46.family, repo copy figures/ideogram/palette_v46.py); micro / mini / satellite classes by the class palette; no EDTA calls or annotation ids in the picture (write RTE / RTE LINE); pale-grey centromeres, no black blocks; text must not overlap; ins = insertion, int = interruption; all numbers come from repo tables (named in the caption stub); catalogue data = SatChain v4.7 runs of record / cat59, non-catalogue inputs (treefiles, LTR_MacSat03 tables, centromere_composition tables, array_definition benchmarks, RTE hits) read in place. Machine: ~1.5 GB free RAM, 3 CPUs: keep every process small. E: read-only; no WSL; no HMMER; no git.
