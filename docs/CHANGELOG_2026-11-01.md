# SatChain v4.7, 2026-11-01: code audit and fixes

Written for: anyone reading the repository history. The two code audits (2026-10-30, 2026-11-01; 90 findings, ids CA/CB/HB/C/D/F in
`F:\data\satellite_FISH\audit_code_2026-11-01\AUDIT_SUMMARY_2026-11-01.md` and the per-stage `FINDINGS.tsv`) were worked through on
branch `fix-audit-2026-11-01`. Every behaviour change is dated 2026-11-01 in the script docstring or a comment. Numbers quoted in
README / METHODS before this date are withdrawn; `benchmarks/regenerate_numbers.py` recomputes them from the tables of a run folder.

## Terminology
- **lineage → sequence cluster.** Step 55 is a k-mer composition clustering of monomers (k-means), not a phylogeny. Column `lineage` is
  now `seq_cluster`, clusters are named `C<K>.<rank>` and described by content, K is chosen by a stated stability rule (`55_K_choice.tsv`),
  a family is assigned only with a ≥ 70 % majority of its sampled copies, and origin claims (5S-derived, RTE-interleaved) are per-family
  homology statements. The old scripts are in `docs/superseded_2026-11-01/`.

## Critical fixes
- **CB01** StringDecomposer 5-kb batch duplicates: `03` `ordered_valid_calls()` drops exact duplicates and resolves overlapping tiles;
  `03b_dedup_copies.py` applies the same rule to older tables (Gar 880,498 → 849,224 rows; Mac and the Gar new candidates likewise);
  `56` and `59` skip exact duplicates on the fly as well.
- **C01 / C02 / F01** MacSat06 is the *C. macrocephalus* 5S rDNA array: `59b` merges the Rfam hits of every class into the barrnap /
  tRNAscan / miRNA tracks, derives the gene length from the data, and calls a gene array only when the bp ratio ≥ 0.7 AND ≥ 50 % of the
  units carry a hit (both criteria swept). `59b` runs inside `59` before tiers, the centromere flag and the clustering (C04).
- **C08 / F02** Overlapping runs resolved in `59` (same family merged, different families trimmed; column `overlap_note`).
- **C03 / C18** Lineage K = 3 hard-coded and unstable → `55c` sweeps K = 2–8 with cross-k, restart and per-haplotype ARI; the fClaHyb run
  falls back to K = 2 (C2.2 = GarSat04 alone).

## Major fixes
- **CB16 / C17 / D08 / CA25** one path resolver `satpaths.py` for the renamed annotation folders; missing masks raise.
- **HB05 / HB06 / HB13** `56` fits one model per haplotype on every table given (`--out` takes several), samples training arrays by id,
  reports mean_post over the trimmed run. **HB07** `gaps.py` `between(..., closed=True)`: a contig join on a tile edge stops a chain.
- **HB08 / CA34 / CA47 / HB19** `val56_runfilter.py` rewritten: any-lag k-mer excess from the 57 scan, bp-weighted mean, background
  outside all candidate intervals, one interface for existing and new runs, benchmark table shipped.
- **C06** G2 is a real cut (identity ≥ 85) and the three gates are swept (`59_gates_sweep_*.tsv`). **C11** EDTA overlap on the clipped
  union; unit class by the 100-bp rule. **C07** `cen_candidate_ror` from the runs of record beside the step-31 flag. **C22** dropped runs logged.
- **C10 / C14 / C15** `59c` / `59d` join on `array_id_v46`, use the dominant class computed by `59`, recompute tiers after the NewSat grouping.
- **C16** `59e` chains consecutive runs by family across candidate intervals, splits `other_annotation` into exon/CDS, intron, ncRNA and
  SSR, sweeps class cut × intrinsic cut × max gap. **C17** `59c` carry rule = feature covers ≥ 50 % of the copy (v4.6 rule kept via
  `SATCHAIN_59C_CARRY=min`; both runs shipped).
- **CB02 / CB12** `58c_new_arrays_02format.py` restored (rule new_array OR (7–39 bp AND TRF cover ≥ 0.5), tier column). **CB03** TRF
  support = coverage of compatible calls. **CB04** null quantile over lag ≥ 7 windows. **CB05** region unit = lag of the best window.
  **CB06** scales united before 58b. **CB09** EDTA merged in 58b. **CB11** tier cut-off sweep. **CB10 / CB14 / CB15** 57 docstring, non-N
  shuffle, right-aligned last window. **CB07 / CB08 / CB19** 03 status columns (library_rebuilt, converged, period_bound_hit, `--seed-multiples`),
  SD reliability flag.
- **HB01** `benchmark_band_test.py` with a null outside arrays, runs and TRF. **HB02 / CA19 / CA43** `53` run mode by default, copy-indexed
  windows for short runs, last window to the run end; sweeps on the runs of record. **HB03** composition + 4-mer = one block. **HB04 / CA26 /
  CA28 / HB14** `53b_o2_summary.py` carries the p/q/BH procedure (BH within each reported set, variant-only transitions, units = runs).
  **HB09 / CA27 / CA29** library cuts from a split-half null; `sweep53b` never overwrites production; `53c` takes arguments. **HB10 / CA21**
  one state definition (class) in 53 and 53b, minority-strand share per run. **CA20 / CA23** per-run seeded RNG, midpoint membership.
  **HB11 / HB12** plots on runs of record, 450 dpi + PDF, no "centromere" in titles.
- **D01 / D19** `61b` scan clipped to the spacer, ramp direction by junction type, w = 1 unit primary, null 1000. **D02 / D23 / D12**
  `chains.py` shared by 61 / 61b / 64 (family_v46, 0 ≤ gap, gap refusal, `~g` ids). **D03 / D16 / D07** `62` requires id(2n) and id(3n),
  per-lag null and max statistic, BH over runs, strand flag. **D04 / D18** `65` shuffled-string null, min_copies 3, plateau choice.
  **D05 / D17** `64` block permutation, squared JSD, BH, `--gaps`. **D06 / D11** `61` junction homology null, satellite_bp / interruption_bp.
  **D10 / D26** `63` adjacency blocks, exact-match null, `plot63` 450 dpi. **D09 / D13** `--sweep` tables per step; STEP_METHOD section 8
  rewritten. **D20** tests as test_* functions (gaps, chains, 61, 61b, 62, 63, 64, 65); environment pinned.

## Not changed / still open
- The *C. macrocephalus* new-candidate 03 table (`new03_Mac`) was lost before the dedup; its 2,076 runs are carried from the v4.6 draft
  (`carried_from` column) until 58c → 03 is rerun.
- 57 / 58 were rerun with the new null for *C. gariepinus* only (`cal58_Gar_rerun/`); the catalogue still rests on the v4.6 candidate set.
- Steps 61–65 have no genome-wide output on disk; 60a specificity against the other haplotype is not implemented.
- The 2026-10 segmentation / topology sweeps on 30 hand-picked intervals are kept for the record only.
- Run folder of record: `F:\data\satellite_FISH\satchain_v4.7_2026-11-01\` (dedup tables, 56, val56, cat59, cat59_minrule, lin55, cat53,
  benchmarks, logs); recovered v4.6 inputs in `runs_of_record_v4.6_2026-10-27\inputs_recovered_2026-11-01\`.
