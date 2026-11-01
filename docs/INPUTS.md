# Inputs expected by SatChain (current layout, 2026-11-01)

The scripts were written for the fClaHyb annotation folders. The haplotype tag `<G>` (e.g. `Gar`, `Mac`) is part of the folder and file
names: `fClaHyb_<G>`. Every script resolves annotation files through `satchain/satpaths.py` (`func_path(FUNC, G, sub, name)`), which
tries the current layout `<func>/fClaHyb_<G><sub>/<name>` first and the pre-2026-10 layout `<func>/<sub>/<name>` second, and raises
when the file is missing (no script calibrates or classifies on a silently missing mask any more; audit CB16 / C17 / D08 / CA25).

Per haplotype, an annotation folder `<func>/` (e.g. `fClaHyb_func_annotation_v1.0/fClaHyb_Gar/`) with:

| path (`<sub>` / `<name>`) | content | used by |
|---|---|---|
| `fClaHyb_<G>_genome/fClaHyb_<G>_LG.fa` (+ `.fai`) | genome | 57, 58b, 03, 59, 61, 61b, 62-65 |
| `fClaHyb_<G>_genome/fClaHyb_<G>_LG.gaps.detgaps.bed` | assembly gaps (N runs and zero-length contig joins) | 56, 53b, 59, 59e, 61, 61b, 64 |
| `fClaHyb_<G>_repeats_EDTA/fClaHyb_<G>.TEanno.sorted.gff3.gz`, `fClaHyb_<G>.TElib.fa` | EDTA annotation and TE library | 58, 58b, 59, 59c, 59e, 53b, 60a |
| `fClaHyb_<G>_rRNA_barrnap/fClaHyb_<G>.rRNA.sorted.gff3.gz` | barrnap rRNA (5S etc.) | 59e, 60a (59b reads the overlay copy) |
| `fClaHyb_<G>_tRNA_tRNAscan/fClaHyb_<G>.tRNA.sorted.gff3.gz` | tRNAscan-SE | 59e |
| `fClaHyb_<G>_miRNA/fClaHyb_<G>.miRNA.curated.sorted.gff3.gz` | miRNA | (overlay copy: 59b) |
| `fClaHyb_<G>_ncRNA_rfam/fClaHyb_<G>.rfam.all_classes.sorted.gff3.gz` | Rfam (5S, tRNA, snRNA, ...; attribute `rfam_acc`) | 59e, 53b (overlay copy: 59b, 55b) |
| `fClaHyb_<G>_ncRNA_all/fClaHyb_<G>.ncRNA_all.sorted.gff3.gz` | merged ncRNA | 59e (optional) |
| `fClaHyb_<G>_genes_lifton/fClaHyb_<G>.lifton.sorted.gff3.gz` | gene models (gene / exon / CDS) | 59e, 60a, 53b |
| `fClaHyb_<G>_TR/GMATA_SSR/fClaHyb_<G>.SSR.{di,tri,tetra,penta_hexa}.bed.gz` | microsatellites | 58, 59e |
| `fClaHyb_<G>_TR/TRF/fClaHyb_<G>.TRF.bed.gz`, `.TRF.mini_7-60bp.bed.gz`, `.TRF.sat_gt60bp.bed.gz` | TRF calls (name `<motif>x<copies>` or `<n>bpx<copies>`) | 58, 58b, val58, 58c |

An overlay folder `<overlays>/<G>/` (`satellite_FISH/annotation_overlays/<G>` on fClaHyb) with `fClaHyb_<G>.SINE.AnnoSINE.merged.bed.gz`
and the copies `fClaHyb_<G>.rRNA.sorted.gff3.gz`, `.tRNA.sorted.gff3.gz`, `.miRNA.curated.sorted.gff3.gz`, `.rfam.all_classes.sorted.gff3.gz`
read by 59b (gene arrays) and 55b (gene overlap of the sampled monomers); `overlay_path()` resolves them.

A pipeline output folder per haplotype (`<out>/fClaHyb_<G>_satellites/` for the existing arrays; any folder for the new candidates, written
by 58c + 03) with: `02_satellite_arrays.tsv` (array_id, chrom, start_0based, end, array_length_bp, period_mode_bp), `03_monomer_copies.tsv`
and `03_array_consensus.tsv` (written by `03_monomer_copies.py`; run `03b_dedup_copies.py` on tables made before 2026-11-01 to remove the
StringDecomposer batch duplicates), `36_characterized_intervals.bed` (v4.5 runs, comparison only) and, for the candidate centromeric
intervals, `31_centromere_summary.tsv` (chrom, cen_start_0based, cen_end).

Step-54 insertion tables (`54b_ins_members.tsv`, `54c_ins_families_classified.tsv`; on fClaHyb in `satellite_FISH/docs/design/step54_insertions/`)
are read by 59c and 59e for the 'missed known-TE insertion' class.

The 1-kb tandem-signal bedgraph of step 48 (`<G>_tandem_signal_1kb.bedgraph`, `satellite_FISH/igv/`) is read by val56.
