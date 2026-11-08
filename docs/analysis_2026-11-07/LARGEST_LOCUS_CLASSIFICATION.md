# Why the largest repetitive locus of a chromosome is not always the candidate centromere (2026-11-07)

Written for: the author and anyone reading the v4.7 catalogue or the ideograms.
Data: `satchain_v4.7_2026-11-01/final/cat59/59_catalogue_{Gar,Mac}.tsv` (final run, 2026-11-02) and `runs_of_record/<G>/centromere_check.tsv`.
Table: `chromosome_loci_top3_v47final.tsv` (every run cluster <= 100 kb apart, any element class, three largest per chromosome).

## What the catalogue classifies today (step 59, 59b, 59d)
Per run, three independent labels:
1. `unit_class` by unit length: micro (<= 6 bp), mini (7-99 bp), satellite (>= 100 bp).
2. `element_class` per FAMILY (59b, before tiers and before any centromere flag): `gene_array:<class>` when annotated genes (5S rRNA, tRNA,
   snRNA, miRNA; barrnap + tRNAscan + Rfam merged) cover >= 0.7 of what one gene per unit would give AND >= 50 % of the units carry a hit;
   `satellite` otherwise. Families with ratio >= 0.2 but short of the gene-array rule keep `element_class = satellite` and carry
   `gene_flag = gene_related:<class>` (tRNA-derived minisatellites: GarSat4282, MacSat1272, MacSat1282, MacSat108).
3. Two centromere columns, both location-only (L1): `putative_centromeric_array_L1` = satellite run overlapping the step-31 interval
   (June TRF/synteny candidate intervals, external to SatChain); `cen_candidate_ror` = run belongs to the largest cluster (<= 100 kb apart)
   of runs with element_class = satellite AND unit_class = satellite on its chromosome.

So gene arrays ARE already excluded from `cen_candidate_ror`. What is missing is a step that looks at the chromosome, not the run:
nothing ranks the loci of a chromosome, labels each by class, and says whether the largest satellite locus and the synteny interval agree.

## What the data show (largest run cluster per chromosome, all classes)
| | Gar (28) | Mac (27) |
|---|---|---|
| largest locus = satellite array inside the synteny interval | 21 | 15 |
| largest locus = satellite array outside the synteny interval | 3 (LG09 GarSat03 20.5 Mb vs interval 9.1 Mb; LG13 GarSat09; LG17 GarSat09) | 3 (LG10 MacSat51; LG12 MacSat21; LG24 MacSat33) |
| largest locus = tandem gene array | 3 (LG04 5S GarSat49 404 kb; LG12 snRNA-U2 GarSat154 487 kb; LG28 tRNA GarSat62 356 kb) | 6 (LG01/LG22/LG25 5S MacSat06; LG03/LG06/LG21 tRNA arrays); LG22's 5S array sits INSIDE the synteny interval |
| largest locus = tRNA-derived gene_related family | 1 (LG05 GarSat4282 269 kb) | 3 (LG02 MacSat1282; LG20 MacSat1272+1271; LG26 MacSat108) |
| synteny interval with no satellite run at all | LG04, LG13, LG17 | LG10, LG22 |

Where the largest satellite-only cluster (`cen_candidate_ror`) and the synteny interval disagree, both columns are set and nothing says which
to trust: Gar LG04/05/09/12/13/17, Mac LG01/03/06/21/22/25 (synteny interval holds a smaller MacSat03 block while the largest locus is elsewhere)
and Mac LG02/10/12/20/24/26.

## Proposed step 59g (chromosome loci), to add after 59d
Input: 59 catalogue + centromere_check. Output: `59g_chromosome_loci_<G>.tsv`, one row per locus (run cluster <= 100 kb, any class), ranked by bp per chromosome.
1. `locus_class` from the bp-majority of its runs: tandem gene array (<class>), tRNA-derived satellite (gene_related), satellite array (unit >= 100 bp), minisatellite cluster.
2. `cen_status` per chromosome, one value:
   - `concordant`: largest satellite-array locus overlaps the synteny interval (Gar 21, Mac 15);
   - `synteny_only`: interval holds a satellite block but a larger satellite-array locus lies elsewhere (list both; e.g. Mac LG01, Gar LG09);
   - `satellite_only`: interval has no satellite run; largest satellite locus elsewhere (Gar LG04/13/17, Mac LG10/22);
   - `gene_array_in_interval`: the interval's largest content is a tandem gene array (Mac LG22, 5S).
3. The ideogram right-hand column and any "largest array per chromosome" figure read `locus_class` and `cen_status`, so a 5S or tRNA locus is
   labelled as such and never as the candidate centromere; gene arrays stay black.
4. Sweep: cluster gap 50 / 100 / 200 kb; minimum locus 50 / 100 kb.
This is a reporting step (no new biology): the candidate centromere remains a location call (L1), and the choice between the synteny
interval and the largest satellite locus where they disagree is left open (needs CENP-A / Hi-C / FISH).
