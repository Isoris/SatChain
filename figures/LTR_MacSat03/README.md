# FISH mapping folder: satellite DNA of the F1 catfish hybrid (updated 2026-10-27)

This folder started as the analysis of one element (the RTE-type LINE that alternates with MacSat03 satellite arrays; its original notes are in `README_MacSat03_RTE_element.md`).
It is now also the **single place for FISH mapping**: every ideogram, every per-array file set, every per-chromosome and per-marker page, the sequences with flanks, and the scripts.
`FILEMAP.md` says where each file is; `FILEMAP.tsv` lists every file with its status.

> **Everything here is in silico (computed from the genome assemblies). No probe has been hybridised yet.** Treat every number as a prediction to test, not as a result.

## The 5-minute orientation

1. **The animals.** The F1 hybrid carries one haploid set from each parent: *Clarias gariepinus* (**Gar**, 28 chromosomes) and *C. macrocephalus* (**Mac**, 27 chromosomes). So the table has 55 chromosome columns: `G01`-`G28`, `M01`-`M27`.
2. **The markers** (the rows of the table) are the repeats we can target: satellite DNA families (53 with an array of 10 kb or more), the 5S and 45S ribosomal DNA, telomere motif blocks, and microsatellites (dinucleotides such as (AC)n, trinucleotides, tetranucleotides).
3. **The question per cell:** *does this marker give a signal on this chromosome, where on it, and could the probe also light up somewhere else?* The in silico matrix gives a prediction for each cell; you record what you expect or what you see.
4. **Cross-hybridisation.** Families whose consensus sequences are at least 80 % identical are treated as one group, because one primer pair would amplify all of them. The figures mark windows of the whole group elsewhere in the genome.

## How to work (suggested order)

1. Open the workbook `FISH_mapping/01_FILL_IN_TABLE/marker_x_chromosome_FILL_IN.xlsx` (being built; the pictures of the same matrices are in `figure/diagnostics/centromeres/matrices/`). Fill it little by little: status, date, initials and notes columns are provided.
2. For a **marker**: open its page in `figure/diagnostics/fish_views/markers/`. It shows every chromosome of both genomes, the marker's arrays with their size in kb, and the other windows where it (or a similar family) also occurs.
3. For a **chromosome**: open `figure/diagnostics/fish_views/chromosomes/<chromosome>_markers.png`. It shows the chromosome to scale with every marker on it and a table (kb, arrays, position, unique or not, cross-hybridisation).
4. For a **centromere array** (54 of them): open its file set in `figure/diagnostics/centromeres/{Mac,Gar}/` (see below). Use the flanks in the `.fa` / `.gb` for locus-specific primers.
5. **Primers:** the checked primer pairs are in `F:/data/satellite_FISH/primer_design_v46/`. Only **7 of 28** pairs are usable satellite pairs (3 fully clean: GarSat4282, GarSat278, MacSat108; 4 with caveats written in that folder's README). 14 pairs fall in gene arrays (tRNA, U1, U2, U6, 5S, miR-430), which are not satellite probes.

## One array = one file set

In `figure/diagnostics/centromeres/{Mac,Gar}/`, names start with the chromosome and the array id, e.g. `LG05_Mac_SAT00380`:

| file | what it shows |
|---|---|
| `.fa` | the sequence: 5 kb flank (lowercase) + the array (UPPERCASE) + 5 kb flank |
| `.gb` | GenBank of the same sequence with features: array, runs, flanks, monomers, RTE hits, assembly gaps |
| `_igv.bed` | track for IGV (genome coordinates) |
| `_rlvl1_context.png` | the array to scale, class composition per block, its place on the ideograms of both genomes, windows of the family and its similar families elsewhere |
| `_chromosome.png` | the chromosome with every marker on it |
| `_marker_matrix.png` / `_marker_matrix_empty.png` | the marker table of that chromosome, filled and blank (being added) |
| `_divmap.png` / `_divmap_vs_centromeres.png` | identity map of the array and its flanks; identity to all other centromeres (being added) |
| `_rlvl2_stack.png` / `_rlvl3_linear.png` / `_rlvl3_topology.png` | internal structure: blocks and domains, every monomer copy, transitions between monomer classes |

## Other folders of the hub (`figure/diagnostics/`)

- `ideograms/`: the whole genome at three levels (level 1: what repeat; level 2: how it is arranged; level 3: which exactly, every label with its own colour and a hatch glyph when colours repeat).
- `ideograms_v45/`: the earlier ideograms that were approved (every run in its IGV family colour), with the code that redraws them.
- `markers_rdna_telomere/`: 45S rDNA, 5S rDNA and telomere sites, one file set per site (being added).
- `markers_microsatellites/`: dinucleotides, trinucleotides, tetranucleotides and mini / micro arrays, in their own folder (being added).
- `genome_views/`, `figures/`, `docs/`, `package/`, `slides/`, `probes_primers/`, `reference_figures/`: see `FILEMAP.md`.

## Colours

Family colours follow the approved reference figure (`figure/diagnostics/reference_figures/centromere_composition.png`): lineage A in reds and oranges, lineage B in blues. Repeat classes follow the thesis: satellite teal, minisatellite mid blue, microsatellite light blue, gene arrays and 5S gold. Where two things would look alike, a small number or a hatch pattern is added. Details: `figure/diagnostics/docs/COLOUR_GUIDE.md` (being written).

## What is not finished yet (honest list)

- Per-chromosome and per-genome marker matrices in the array folders, the per-array divergence maps, the 45S / telomere / microsatellite sets, the fill-in workbook and the numbered `FISH_mapping/` bundle are being built right now. `FILEMAP.md` is regenerated (`python -P make_filemap.py`) when they land, and shows what exists.
- Microsatellite and minisatellite pages are first drafts.
- The thesis slides (G10, G11, G11d / G11e ...) are in `F:/data/thesis_defense_slides_tpl/`; release R34 is being rebuilt with the new colours.
- Not used for primer design: motif views (removed as stale), anything marked outdated, superseded or check in `FILEMAP.tsv`.

## Names

- Families have two names: the pipeline id (`GarSat04`, `MacSat06`) and the manuscript name (`CgaSat001-216`, `CmaSat002-214`; Cga = *C. gariepinus*, Cma = *C. macrocephalus*, the number after the dash = unit length in bp).
- Position classes: **C** centromeric, **I** interstitial, **T** terminal (within 5 % of the chromosome end).
- `rlvl1 / rlvl2 / rlvl3`: resolution level 1, 2, 3 (see above).
- The assembly: haplotype FASTA files are under `F:/data/fClaHyb_func_annotation_v1.0/fClaHyb_{Gar,Mac}/fClaHyb_{Gar,Mac}_genome/` (read-only).
