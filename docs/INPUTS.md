# Inputs expected by SatChain (current layout)

The scripts were written for the fClaHyb annotation folders. The haplotype tag `<G>` (e.g. `Gar`, `Mac`) is part of the file names:
`fClaHyb_<G>`. Rename your files to this pattern, or edit the paths in the scripts (they are near the top of each script).

Per haplotype, an annotation folder `<func>/` with:

| path | content |
|---|---|
| `_genome/fClaHyb_<G>_LG.fa` (+ `.fai`) | genome |
| `_genome/fClaHyb_<G>_LG.gaps.detgaps.bed` | assembly gaps |
| `_repeats_EDTA/fClaHyb_<G>.TEanno.sorted.gff3.gz`, `fClaHyb_<G>.TElib.fa` | EDTA annotation and TE library |
| `_rRNA_barrnap/fClaHyb_<G>.rRNA.sorted.gff3.gz` | barrnap rRNA (5S etc.) |
| `_tRNA_tRNAscan/fClaHyb_<G>.tRNA.sorted.gff3.gz` | tRNAscan-SE |
| `_miRNA/fClaHyb_<G>.miRNA.curated.sorted.gff3.gz` | miRNA |
| `_ncRNA_rfam/fClaHyb_<G>.rfam.all_classes.sorted.gff3.gz` | Rfam (snRNA ...) |
| `_genes_lifton/fClaHyb_<G>.lifton.sorted.gff3.gz` | gene models (exon / CDS) |
| `_TR/GMATA_SSR/fClaHyb_<G>.SSR.{di,tri,tetra,penta_hexa}.bed.gz` | microsatellites |
| `_TR/TRF/fClaHyb_<G>.TRF.bed.gz`, `.TRF.mini_7-60bp.bed.gz`, `.TRF.sat_gt60bp.bed.gz` | TRF calls (name `<motif>x<copies>` or `<n>bpx<copies>`) |

An overlay folder `<overlays>/` with `fClaHyb_<G>.SINE.AnnoSINE.merged.bed.gz` (and copies of the rRNA / tRNA / miRNA / Rfam files).

A pipeline output folder per haplotype (`<out>/fClaHyb_<G>_satellites/` for existing arrays, or any folder for new candidates) with:
`02_satellite_arrays.tsv` (array_id, chrom, start_0based, end, array_length_bp, period_mode_bp), `03_monomer_copies.tsv` and
`03_array_consensus.tsv` (written by `03_monomer_copies.py`), and, for putative centromeric arrays, `31_centromere_summary.tsv`
(chrom, cen_start_0based, cen_end).
