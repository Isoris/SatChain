# DRAFT steps, narrow lanes and library release — decisions of 2026-11-09

Context: evaluation of `SatChain_20_Worklanes_Specs_v4_2026-11-08.zip` (do not dispatch 20 lanes; run SC-10 and SC-01 narrow; decide the
DRAFT steps first). Graduation freeze (2026-11-06) and the measurement-purpose rule apply: every step below is either run with a stated
purpose or labelled exploratory.

## Decision ledger

| step | decision | reason (purpose test) | changed |
|---|---|---|---|
| 53c_library_pool | EXPLORATORY | no result, figure or thesis claim uses a pooled topology library; per-haplotype 53b libraries are used | register, README |
| 60a_tx_kmers | EXPLORATORY | the transcription test was done 2026-10-23 outside the chain (poly(A) data cannot rule out non-polyA transcripts) | register, README |
| 61_unit_boundaries | EXPLORATORY, removed from run_satchain.sh | no result needs base-pair edges (composite units ~13.6 kb; the 56 edge is 9-154 bp inside the satellite) | workflow (commented), register, README, STEP_METHOD |
| 61b_jsd_boundaries | EXPLORATORY, removed from run_satchain.sh | same as 61 | same |
| 62_hor_detect | RUN NARROW | re-test the one HOR claim of the thesis (Ch3 Table 3.9: MacSat1272 on Mac LG20), which rests on v4.5 self-similarity matrices | driver `run_hor_narrow.py`, results below |
| 63_expansion_tracks | EXPLORATORY | no claim depends on it | register, README |
| 64_object_profiles | EXPLORATORY | no claim depends on it | register, README |
| 65_hor_mining | RUN NARROW (cross-check of 62) | second, independent HOR signal on the same arrays | same driver |
| 59g_run_axis (new, SC-10 narrow) | RUN | evidence for the MacSat06 hand decision (59g: unresolved) | `satchain/59g_run_axis.py`, results below |
| export_library_v1 (new, SC-01 narrow) | RUN | deposit-ready family reference set, no new ids, no HMMs | `satchain/export_library_v1.py` |

Code of the exploratory steps and their synthetic tests stay in the repository.

## SC-10 narrow: 5S sequence axis per run (`final/run_axis_2026-11-09/`)

Same references, k-mer screen, Smith-Waterman scoring and cut-offs as 59g (intact = identity >= 0.90 and reference coverage >= 0.80),
applied to every run consensus and to up to 30 confident copies per run (doubled; a copy can start inside the gene).

| family | runs | runs with intact 5S (bp share) | copies with intact 5S | median copy identity |
|---|---|---|---|---|
| GarSat49 (control: gene array) | 32 | 31 (0.999) | 0.981 | 1.00 |
| GarSat04 (control: 5S-derived satellite) | 58 | 0 (0.0) | 0.0 | 0.815 |
| MacSat06 | 78 | 34 (0.658) | 0.355 | 0.867 |

MacSat06 splits by locus, not along a continuum: LG01, LG02, LG08, LG20, LG25 = intact 5S gene + spacer runs (34 runs, 1.68 Mb, median
94 % of copies intact); LG03 and LG22 = 5S-derived units only (44 runs, 0.87 Mb, 0 % intact copies, consensus 207-215 bp). The derived
units are closer to MacSat06's own intact units (median identity 0.82) than to GarSat04 (0.74): a diverged 5S-unit lineage inside Mac,
not a Mac copy of GarSat04. Sequence alone cannot separate a pseudogene array from a young 5S-derived satellite.
59g labels are unchanged (no rule was added). Hand decision 2026-11-09 (author): MacSat06 = gene_array, with the LG03/LG22 derived loci noted (`docs/decisions/hand_decisions.tsv`).

## SC-01 narrow: library release v1.0.0 (`final/library_release/SatChain_satellite_library_v1.0.0/`)

2,062 families in `SatChain_consensus.fa` (satellite 2,023, gene_array 25, unresolved 11, gene_derived_satellite 3; MacSat06 by hand decision) + 987 low-copy
candidates; 4,343 aliases; every model `not_attempted`. Deviations from the proposed format are listed in its release notes.

## 62 / 65 narrow HOR re-test (`final/hor_narrow_2026-11-09/`)

Arrays: runs of the dominant family overlapping the candidate centromere interval of the thesis Table 3.9 chromosomes (9 HOR-like,
4 plateau controls; `arrays.tsv`). 62 = lag-identity profile against a shuffled-order null (BH q <= 0.05), swept over 36 settings
(min identity 60/65/80, pairs 75/150, permutations 200/1000, min excess 0.02/0.03/0.05); 65 = HiCAT-like monomer-pattern mining
with a shuffled null per identity threshold. Copy tables (`copies_*.tsv`) are kept locally only.

| array (role) | 62 call (monomers, ~bp) | 62 lag identity (lag 1 -> best) | 65 top pattern: length, copies, coverage (null p95) |
|---|---|---|---|
| Mac LG20 MacSat1272 (HOR-like) | E03313: 5, ~1,060 bp; E03314 not called (lag-5 id 0.98 but fails the 2n/3n check) | 0.33 -> 0.98 (lag 5), both runs | 5-mer, 62 copies, 0.47 (0.02); both runs |
| Gar LG25 GarSat03 (HOR-like) | E02761: 3, ~678 bp | small excess | 3-mer, 46 copies, 0.15 (0.02) |
| Gar LG10 GarSat09 (HOR-like) | none (k-mer 1,571 bp interval in E01264/E01265) | small excess | 7-mer, 13 copies, 0.17 (0.08) |
| Mac MacSat03, 7 HOR-like arrays | LG05 E00828: 2 (~450 bp); others none | 0.80-0.93 -> 0.88-1.0 | 2- or 3-mers on 6 arrays; LG23 12-mer (4 copies) |
| Gar LG11 GarSat03 (control) | none | — | 2-mer, 407 copies, 0.36 (0.22) |
| Gar LG09 GarSat04 (control) | none | — | 2-mer, 35 copies, 0.03 (0.08) |
| Mac LG13 MacSat03 (control) | E02243: 3 (~678 bp) | 0.80 -> 0.93 | 3-mer, 11 copies, 0.04 (0.01) |

Call counts are flat across the 62 sweep (Gar 1-2, Mac 3 at every setting).

Reading:
- **MacSat1272 on Mac LG20 is confirmed as a HOR.** Both methods agree on a 5-monomer period. The lag-1 identity is 0.33,
  because neighbouring monomers differ, while the lag-5 identity is 0.98. That contrast is the HOR signature, and it is the
  strongest signal in either genome. The thesis claim stands with the v4.7 runs.
- The 2- and 3-monomer patterns in MacSat03, GarSat03 and GarSat09 are weak. 62 also calls one on the plateau control
  Mac LG13, and 65 reports a 2-mer covering 36 % of the plateau control Gar LG11. Short patterns therefore do not separate
  HOR-like arrays from controls and are not evidence of HOR structure. No new HOR claim is made.
- 62 and 65 stay exploratory outside this re-test. 65 cannot serve as an independent check for short periods.

`o62/62_hor_calls_<H>.tsv` carried the v4.5 copy-table label (retired CgaSat/CmaSat names); `fix_hor_family_names.py` renamed that
column to `legacy_name_do_not_use` and added the published `family`.

## Dead code removed (ledger row of 2026-11-07)

`patch_dead_code_2026-11-09.py`: removed the unused `read_fasta` and `LIB_EPS` from 03_monomer_copies.py, the unused `--gaps`
option from 59_final_catalogue.py (no caller passes it; gap segments come from the 56 BED names) and the unused `merged_bp` from
runs_of_record.py. Outputs are unchanged (nothing read these names). `02 --threads` and the `11` mismatch variants belong to the
frozen v4.5 pipeline (`satellite_FISH/pipeline_v4_5_fixed`, deposited code) and stay as they are.

## README numbers and benchmark banners (final run)

- `benchmarks/regenerate_numbers.py` on `final/` reproduces `NUMBERS_2026-11-02.tsv` except one row: the bp-windowed 53 count
  included unsegmentable runs (1,480 / 1,883 -> 359 / 450; fixed). Current table: `benchmarks/NUMBERS_2026-11-09.tsv`.
- README: every number moved from the 2026-11-01 run to the final run; retired names replaced by published ids (legacy names in
  brackets); status section updated (no carried Mac runs, exploratory steps, narrow HOR result, provisional ids in 59/59b/59c).
- `plot_benchmarks.py` banner 7 claimed every K >= 3 fails the ARI rule; in the final run K = 4 and 6 are stable and fail on silhouette.
  The banner now names the failing criterion from `55_K_choice.tsv`; the cluster-size panel labels every colour. Re-rendered to
  `final/benchmarks/figures_2026-11-09/` (older figures kept).
- Also in `plot_benchmarks.py`: the gene-array key now covers snRNA U4/U5 and the "other" grey (U4 was drawn grey without a key entry);
  the val56 "old sequence kept" panel states that its > 100 % values are a per-run sum over overlapping runs, not a union.
  Overlapping point labels in the gene-array panel are left as they are.
