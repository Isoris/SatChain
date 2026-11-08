# Classification decision and measurement ledger (2026-11-07)

Written for: the author. Follows the review rule agreed today: one scientific decision at a time; every measurement has a declared
purpose (decides something, supports a result or figure, is a distinct QC, or is explicitly exploratory) or it is not part of the
production pipeline; an operation is kept, simplified, corrected or documented, not replaced by another operation by default.

## 1. The decision: what is inside the repeat unit of a family

One decision per family, four outcomes, from three axes that stay separate columns (`59g_unit_content_<G>.tsv`, step
`satchain/59g_unit_content.py`, outputs in `satellite_FISH\satchain_v4.7_2026-11-01\final\unit_content_2026-11-07\`):

| axis | column | what it measures | source |
|---|---|---|---|
| sequence identity | `sequence_axis` | the unit consensus aligned to the host's own annotated gene copies (5S, tRNA, snRNA U1-U6, miR-430, other miRNA, 5.8S/18S/28S, Histone3, UnaL2, other Rfam): intact_gene (identity >= 0.90, gene coverage >= 0.80) / gene_derived (>= 0.70, >= 0.50) / none | 59g (new, targeted check) |
| annotation | `element_class_59b`, `ratio_59b` | annotated-gene overlap of the family's runs vs one gene per unit; threshold-sensitive when the ratio is within 0.1 of the 0.7 cut | 59b (existing) |
| organization | `unit_organization` | gene_plus_spacer (unit <= 2.5 x gene length) / gene_in_long_unit / no_gene | 59g |

Outcome: **gene_array** = intact gene AND 59b calls it AND not threshold-sensitive; **gene_derived_satellite** = gene-derived unit and no
annotation call; **satellite** = no gene homology and no annotation call; **unresolved** = the axes disagree, or the 59b ratio is
threshold-sensitive. TE content, other-family content and composite flags are reported (`te_frac`, `othersat_frac`,
`composite_sat_pair`, `segments`) and do not enter the decision.

### Result on the final v4.7 catalogue

| | Gar families (Mb) | Mac families (Mb) |
|---|---|---|
| satellite | 1,599 (12.90) | 1,451 (13.00) |
| gene_array | 20 (2.33) | 15 (1.79) |
| gene_derived_satellite | 2 (5.85: GarSat04 5.80, GarSat4234) | 1 (0.02: MacSat6314) |
| unresolved | 5 (0.62) | 9 (3.18, of which MacSat06 2.55) |

What the organization axis adds: of the 35 gene_array families only GarSat49 (5S) and MacSat315 (U1) are gene + spacer units; every
other one is a long unit (590-1,604 bp) carrying an intact tRNA / snRNA / miR-430 / 5S gene together with TE or other-family
segments. In the text these are "tandem units carrying an intact gene", not "gene arrays" in the 5S-rDNA sense. GarSat04 (the largest
Gar centromeric family) is a gene + spacer unit whose 5S part is 78 % identical over 81 % of the gene with no annotation: the
sequence-level evidence for "5S-derived" that the pipeline did not have before.

### The 14 unresolved families and why

| family | bp | 59b | sequence axis | reason |
|---|---|---|---|---|
| MacSat06 | 2.55 Mb | gene_array:5S (ratio 0.726) | intact 5S, gene + spacer | ratio within 0.1 of the cut. The unit says gene array; the annotation says borderline. Decide by hand: with 100 % identity over the whole 5S gene and a 214-bp gene + spacer unit the sequence evidence is the stronger one. |
| GarSat62, GarSat847, MacSat319, MacSat161, MacSat6296 | 8-387 kb | gene_array:tRNA | gene_derived (0.77-0.88) | annotation calls intact tRNAs, the unit consensus holds a diverged one: tRNAscan hits sit on copies, the family consensus has drifted. Family-level label hides within-family variation (see 3). |
| MacSat170 (Mac), GarSat4256 | 15-154 kb | gene_array:tRNA | none (coverage 0.28-0.45) | the best tRNA hit covers less than half a gene: annotation overlap came from gene copies not represented in the unit consensus. |
| MacSat108, GarSat4243 | 11-117 kb | satellite / gene_related (ratio 0.64) | intact tRNA (1.00 / 0.96) | the unit carries an intact tRNA but the annotation ratio is below the cut: tRNAscan missed copies, or units vary. |
| MacSat1287, MacSat320 | 18-84 kb | gene_array:tRNA | gene_derived (identity 1.0, coverage 0.53-0.62) | half a tRNA at full identity: truncated gene in the unit. |
| MacSat6328 | 39 kb | gene_array:U4 (ratio 0.753) | intact U4 | ratio threshold-sensitive. |

## 2. Measurement ledger (keep / simplify / correct / document / exploratory)

| measurement | step | purpose today | verdict |
|---|---|---|---|
| annotated-gene overlap ratio and unit share | 59b | decision axis 2 | keep; document that it is coverage-vs-expectation, not "% of the array covered by genes" |
| unit vs gene sequence identity / coverage | 59g | decision axis 1 (targeted check, 2 min) | keep |
| unit organization (gene + spacer vs long unit) | 59g | decision axis 3; wording of results | keep |
| TE / other-family segments of the unit, composite flags | 59g | reported columns; no decision | document; use only when a figure needs them |
| PCA, k-means / Ward, silhouette, cluster x label tables, PCA figure | 59g | exploratory: shows the labelled groups are separable (K = 2, silhouette 0.63) | exploratory; not in the production chain |
| sequence clusters of monomer k-mer profiles | 55a-d, 59f | figure grouping (C2.1 / C2.2 in the ideogram legend); the "lineage" claim was dropped | keep only if a figure uses it; otherwise exploratory. It does not separate satellite types (K = 2 fallback, one cluster is GarSat04 alone). |
| topology, segments, second-order test | 53, 53b, 53b_o2 | descriptive result (L2) in the repetitive part (ED Fig R3 lineage) | keep as result; not a classifier |
| gates G1-G3 | 59 | QC of runs | keep |
| TE layer, interruptions | 59c, 59e | results (TE insertion share, interruption classes) | keep |
| centromere columns L1 and ror | 59, runs_of_record | location call; disagree on 18 chromosomes | document; a chromosome-level reporting step (loci ranked and labelled by content) would make the disagreement visible; no new biology |
| tiers | 59, 59d | figure legends | keep |
| val56 sequence check, band test | val56, benchmark_band_test | QC of runs vs null | keep |
| 61 unit boundaries, 61b, 62-65, 53c, 60a | drafts | no stored output | not in production; remove the quoted "first genome-wide result" of 61 from the docs |
| dead code found by the function trace | 03 read_fasta, LIB_EPS; 59 --gaps; runs_of_record merged_bp; 11 mismatch variants; 02 --threads | none | delete or wire up before release |
| hard-coded conclusions in benchmark banners | plot_benchmarks | figure text | correct: read from the tables |

## 3. Known design limitation, not fixed today

Classification is per family. A family can hold arrays whose units carry an intact gene and arrays whose units have drifted (the
five "gene_derived in the consensus, intact in the annotation" cases above). A per-run version of axis 1 (align each run's own
dominant-class consensus) would show this; it is a one-line change in 59g and about the same run time, but it is a new measurement and
is not run until a result needs it.

## 4. The evolutionary axis is not measured

Nothing in the pipeline supports "expansion" or "contraction" claims. What exists: mean copy identity per family (a homogenization
proxy, column `copy_identity`), the 53b topology, and the 63 expansion tracks (draft, never run). Default label for the third axis is
"unresolved history". The author's references (TRASH, Wlodzimierz 2023; SatXplor, Volaric 2025; Satellite DNA Genomics 2025) were not
checked by me and are not cited here.

## 5. Added later on 2026-11-07: the per-family satellitome table (59h)

| measurement | step | purpose | verdict |
|---|---|---|---|
| per-family metrics (unit, genome share, copies, arrays and their length classes, loci, clustered/dispersed, centromere share, divergence, GC/palindromy, gene/TE association, low-copy flag) | 59h | results table (manuscript, comparison with other satellitomes; the metrics of Garrido-Ramos et al. 2025 sections 3-7) | keep; aggregates of existing tables only |
| library view (shared / haplotype-specific, copy ratio) | 59h | result | keep; the identity cut is swept (0.70-0.90) because GarSat04 <-> MacSat06 sits at 0.77 and the shared bp share moves from 46 % to 3 % across the sweep; homologs that exist only as v4.5 consensuses without a run of record are flagged, not counted as shared |
| abundance rank curve + array-length classes figure | 59h | figure | keep |

Numbers (final v4.7): Gar 1,602 satellite families = 18.75 Mb = 1.95 % (GarSat04 0.60 %; top 5 = 65 %); Mac 1,453 = 13.14 Mb = 1.47 % (MacSat03 0.38 %; top 5 = 52 %). Below the reporting line (< 10 copies or < 1 kb): 550 / 446. Single-locus families 1,411 / 1,276; dispersed (>= 5 loci) 58 / 40, dominated by the 42-bp MacSat01/02 and GarSat01. Satellite bp in arrays >= 50 kb: 53 % (Gar) vs 14 % (Mac). Median divergence 6.8 / 6.7 %. Library at 0.80: Gar 184 shared (16 % of bp) / 24 homolog-without-run / 1,418 specific; Mac 196 (55 %) / 21 / 1,259. GarSat03 and GarSat91 have no Mac homolog at any cut; MacSat03's best Gar homolog is GarSat156 (0.86, 69 kb); GarSat09/10 match MacSat03 at 0.81-0.82.

## 6. 2026-11-08: how the library was matched, and the check of the shared pairs

How 59h matched: one alignment per family (the unit consensus of the family's largest run, doubled, vs every family consensus of both
haplotypes after a 15-mer screen; parasail Smith-Waterman); shared = best other-haplotype family at >= 0.80 identity, alignment covering
>= 50 % of BOTH units (the both-unit rule was added on 2026-11-08: the one-sided rule let a shared tRNA / U1 segment inside a long unit
pass, GarSat4282 vs MacSat1286, MacSat315 vs GarSat4239, MacSat1272 vs GarSat155), homolog with a run of record, cut swept 0.70-0.90.
After the correction: Gar 172 shared (14.8 % of satellite bp at 0.80; 44.6 % at 0.70), Mac 172 (51.1 %; 57.4 % at 0.70).

`library_pair_check.py` checks each pair at four levels against nulls (`final\library_pair_check_2026-11-08\pairs_check.tsv`, one page
per pair): (1) monomer: unit vs unit alignment, and 150 confident copies of each family aligned to both units and to 20 random units of
the other haplotype; (2) composition: 4-mer JSD vs 300 random cross-haplotype unit pairs; (3) organisation from 59h; (4) 53b topology
of the largest run (satellite-class, existing-source runs only). Verdict rule: same_family = consensus >= 0.80, copies within 10 points of
their own unit, above 99 % of the null, reciprocal; sister_family = consensus >= 0.70 and above 95 % of the null; else not_supported.

| pair | consensus | copies A->own / A->B ; B->own / B->A | null pct | reciprocal | verdict | reading |
|---|---|---|---|---|---|---|
| GarSat45 - MacSat22 | 1.000 | 1.00 / 1.00 ; 1.00 / 1.00 | 100 | yes | same_family | identical unit in both haplotypes |
| MacSat15 - GarSat158 | 0.888 | 0.91 / 0.87 ; 0.93 / 0.85 | 100 | yes | same_family | |
| MacSat06 - GarSat49 (5S) | 0.848 | 0.88 / 0.80 ; 0.98 / 0.85 | 100 | yes | sister_family | same origin, ~15 % apart; organisation differs (78 arrays / 15 loci vs 32 / 5) |
| GarSat04 - MacSat06 | 0.790 | 0.88 / 0.77 ; 0.88 / 0.76 | 100 / 95 | no | sister_family | below the 0.80 cut; the 5S-related trio GarSat04 / GarSat49 / MacSat06 |
| GarSat09 - MacSat03 | 0.819 | 0.80 / 0.81 ; 0.81 / 0.79 | 95 / 100 | no | sister_family | GarSat09 copies are as close to the MacSat03 unit as to their own |
| GarSat10 - MacSat03 | 0.781 | 0.83 / 0.78 ; 0.81 / 0.78 | 100 | no | sister_family | |
| MacSat03 - GarSat156 | 0.856 | 0.80 / 0.78 ; 0.93 / 0.84 | 90 / 95 | no | not_supported | MacSat03 is heterogeneous (own copies 0.80); GarSat156 is a 69-kb family |
| GarSat154 - MacSat169 (U2 arrays) | 0.828 | 0.99 / 0.82 ; 0.99 / 0.84 | 100 | no | sister_family | |
| MacSat01 - GarSat246, MacSat04 - GarSat27, GarSat02 - MacSat1021, GarSat07 - MacSat166 (42-54 bp units) | 0.82-0.90 | gaps 1-10 pts | 100 | no | sister_family | |
| MacSat05 - GarSat196 | 0.839 | 0.85 / 0.77 ; 0.93 / 0.85 | 90 / 100 | no | not_supported | |

What this says: the Mac major centromeric family MacSat03 corresponds to a GROUP of minor Gar families (GarSat09, GarSat10, GarSat156,
each ~0.80 to it) and not to GarSat03, which has no Mac homolog at any cut; "reciprocal" therefore fails by construction for a
many-to-one relation and should not decide alone. The null uses 20 random units (coarse: 95 = one unit above); it separates true
homology from chance but not sister from same family; that distinction rests on the copy gap and the user's judgement.

## 7. Hand decision 2026-11-09: MacSat06 = gene_array

Rule output unchanged (59g: `unresolved`, because the 59b ratio 0.726 is within 0.1 of the cut). Author's decision: **gene_array**, recorded
in `docs/decisions/hand_decisions.tsv` (the one place hand decisions live; exporters read it). Evidence: the unit is an intact 5S gene +
spacer (identity 1.0, coverage 1.0); per run (`59g_run_axis.py`, `final/run_axis_2026-11-09/`) 34 of 78 runs (65.8 % of bp) carry an
intact 5S gene in their consensus, on LG01, LG02, LG08, LG20 and LG25. Locus note kept with the decision: the runs on LG03 and LG22
(44 runs, 0.87 Mb) hold 5S-derived units without an intact gene in any sampled copy, closer to MacSat06 units (median identity 0.82) than
to GarSat04 (0.74); pseudogene array or young 5S-derived satellite cannot be told apart from sequence. Controls: GarSat49 98 % intact
copies, GarSat04 0 %.
