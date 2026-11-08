# The element interleaved with MacSat03 arrays (EDTA Mac TE_00003409_INT), 2026-10-22

Question: what is the ~7.5-kb "LTR/unknown" element that alternates with MacSat03 (lineage A) satellite every ~12.5 kb, and where does it come from?

## Results so far
> **2026-10-23 correction: results 2–4 below are withdrawn.**
> - Why: the consensus contains about two lineage A satellite monomers (2,929–3,388; MacSat03 80%) and several internal repeats, and the first copy-chaining was not collinear. Satellite hits and fragments were therefore merged into "copies" (spans up to 486 kb).
> - What stands: result 1 (RTE-type RT, protein-level hits).
> - What replaces them: copies are now defined on the satellite-free, RT-bearing core 5,560–7,400. See `core/summary.txt` and `tree/tree_stats.txt` (code: 05_mask.sh, 06_core.sh).

1. **It is an RTE-type LINE, not an LTR element** (`domains/`).
   - The library entry has only an internal region (8,164 bp) and no LTR.
   - Its one clear domain is a reverse transcriptase whose best match is the RTE-clade LINE RT: 38% identity over 194 aa, E = 3e-34, reverse strand at 6,219–6,746 bp. Other LINE clades match at E ≤ 1e-19.
   - There is no Ty3/Gypsy or Copia RT, protease, integrase or RNase H.
   - Method: tblastn of the reference domain consensuses (CONS columns of the TEsorter REXdb / GyDB / Kapitonov profiles; no HMMER) against the element.
   - Only about 920 bp (11%) matches satellite families, and none of it is MacSat03, so the consensus is not a satellite chimera.
   - Limitations: the consensus has only one ORF ≥ 150 aa; the reference set has no LINE endonuclease profile.
2. **It is in both genomes** (`copies/summary.txt`; dc-megablast of the consensus, chained hits).

   | | Mac | Gar |
   |---|---|---|
   | copies of any length | 5,682 (11.1 Mb) | 11,038 (15.2 Mb) |
   | near-full copies (≥ 70% of the consensus) | 166 | 246 |
   | median identity of near-full copies to the Mac-built consensus | 77% | 73% |

   It predates the Gar–Mac split.
3. **It is centromeric only in lineage A centromeres** (`copies/by_centromere_lineage.tsv`).

   | Centromeres | Near-full copies inside the centromere ± 1 Mb |
   |---|---|
   | Gar, lineage A | 85 of 90 (94%) |
   | Gar, lineage B | 3 of 129 (2%) |
   | Mac, lineage A | 157 of 160 (98%) |

   On Gar lineage-B chromosomes the copies lie about 10 Mb from the centromere, in clusters, and 71 of 129 sit next to interstitial lineage A arrays (GarSat10/09/03).

   The element follows lineage A satellite wherever it is (L1).
4. **Homologous chromosome pairs** (`copies/homologous_chromosomes.tsv`): Mac chromosome → Gar partner (GENESPACE dominant), with copy counts and the satellite family of each centromere.

## Interpretation (to test)
- **L2:** an old RTE-like element that was associated with lineage A satellite in the common ancestor of the two species and is kept, or repeatedly amplified, inside lineage A arrays. Lineage A centromeres became interleaved arrays; lineage B centromeres did not.

## Next
- Ages: pairwise divergence of the near-full copies within and between species; a tree of copies (Gar and Mac copies interleaved = ancestral; species-specific clades = independent expansions).
- Per-locus orthology: anchor-bounded comparison at the Mac LG25/08/16/15/09 centromeres and their Gar partners.
- Outgroups: C. fuscus, C. dussumieri, H. fossilis.
- Endonuclease domain: needs a LINE APE/EN reference set.
- Full-length structure: 3' tail, target-site duplications.

## Corrected results (2026-10-23, RT-bearing core 5,560–7,400)
- **Core copies** (`core/summary.txt`):

  | | Mac | Gar |
  |---|---|---|
  | core copies | 316 | 1,244 |
  | inside a centromere (±1 Mb) | 212 (67%) | 110 (9%) |
  | next to lineage A satellite | 57% | 5% |

  All 212 Mac centromeric copies are in lineage A (MacSat03) centromeres.
- **Tree** (`tree/`; all 316 Mac + 400 Gar stratified, IQ-TREE -fast, SH-aLRT 1000; the full 1,560-copy UFBoot tree crashed for lack of WSL memory):
  - Only 3.6% of copies have their nearest relative in the other species (random labels: 49%).
  - Large single-species clades: Gar 266 copies, Mac 241 copies (SH-aLRT ≥ 90).
  - Median patristic distance: Gar–Gar 0.10, Mac–Mac 0.17, Gar–Mac 0.18.
- **Locality within species:**
  - Mac centromeric copies: 95% have a centromeric nearest relative (chance 56%), 65% on the same chromosome.
  - Gar: no such structure (centromeric 6 of 34, about chance; same chromosome 11%).
- **Interpretation (L2):**
  - The element was present in the common ancestor, but the copies seen today were made after the split, separately in each species.
  - Gar: a dispersed, younger retrotransposition burst (homogeneous, mostly on chromosome arms).
  - Mac: copies inside lineage A centromeres grew locally, chromosome by chromosome. This fits co-amplification with MacSat03 as part of a ~12.5-kb satellite + element unit (array-level duplication and homogenization), not fresh insertions.
- **To strengthen:**
  - the full tree when memory allows;
  - ages per clade;
  - the ~12.5-kb unit as a higher-order repeat (element + satellite period);
  - expression (k-mers, pending memory);
  - outgroups.
- **Composite repeat unit in Mac centromeres (L1, 2026-10-23):**
  - Consecutive centromeric core copies are a median 13.6 kb apart (IQR 12.5–15.3 kb; 188 intervals).
  - Per-chromosome period: 12.4 kb (LG07) to 15.3 kb (LG02, LG16).
  - Tree distance between Mac centromeric copies, by separation:

    | Separation | Median tree distance |
    |---|---|
    | same chromosome, ≤ 25 kb | 0.001 |
    | same chromosome, 25–100 kb | 0.002 |
    | same chromosome, 100 kb–1 Mb | 0.004 |
    | same chromosome, > 1 Mb | 0.173 |
    | different chromosomes | 0.156 |

  - Mac lineage A centromeres are tandem arrays of a ~13–15-kb element + MacSat03 unit, homogenized within each centromere (concerted evolution). Each chromosome has its own unit variant.
  - In Gar the element is dispersed and not part of the centromere unit.
- **Expression, k-mer test** (`expression/kmer_depth.tsv`; 2026-10-23):
  - Method:
    - Element-specific 31-mers: in ≥ 3 core copies, with every genome occurrence inside a core copy (Gar 8,791 k-mers, Mac 5,854).
    - Scale: single-copy EF1a k-mers (gapdh/rpl13a not matched by name; to add).
    - Reads: first 2 M spots of each public run.
  - Gar (liver ×4, spleen ×4) and most Mac samples (3–13, 16, 24 dpf whole body; 1-month liver):
    - element k-mer depth 0–0.06 per million reads, 3,000–80,000× below EF1a;
    - only 0–13% of element k-mers seen;
    - essentially silent.
  - Mac 14 and 17 dpf, both replicates:
    - 0.11–0.36 per million, 27–49% of element k-mers seen, 49–470× below EF1a;
    - EF1a itself is low in these libraries (17–50 vs ~100 per million).
  - Caveats:
    - The elevated samples need a genomic-DNA contamination control (single-copy intergenic/intronic k-mers in the same runs) before calling this expression.
    - Libraries are unstranded, so sense vs antisense is not tested.

## Why the period differs between chromosomes: decomposing the unit (2026-10-23; `unit/`, code 09, 10)
- **Units:** 177 consecutive centromeric units (one element core to the next, 8–25 kb) on 8 Mac chromosomes. Each unit split into satellite (validated runs), element (consensus hits, excluding its built-in satellite) and other.
- **Two unit types:**

  | Type | Chromosomes | Element in the unit | Satellite | Other |
  |---|---|---|---|---|
  | I | LG07, 04, 01, 23, 02 | complete, 5.6–7.4 kb (73–76 of 82 bins) | 4.4–6.4 kb | ≤ 2.9 kb |
  | II | LG08, 25, 16 | truncated to the RT core, 1.7–2.5 kb | 3.9–5.4 kb | ~7.7 kb |

  - Within type I, the period tracks satellite block length (LG07 4.8 kb satellite, 12.3-kb period; LG23 6.4 kb, 13.9 kb).
- **Calibrated nodes (100-bp element bins):** 16 present on every chromosome, all in the RT core; 63 chromosome-variable (the type II truncation). See `unit/node_presence.tsv`.
- **The type II insert is one sequence shared by three non-homologous centromeres:**
  - ~5.2–6.1 kb segments from LG08, LG25 and LG16 are 98–99% identical within a centromere and 92–94% between centromeres.
  - Not identified: best hit is EDTA Unknown TE_00000238 (1.2–1.8 kb), plus fragments of DNA transposons (hAT, Harbinger, Helitron).
  - Interpretation (L2): the insert entered the unit once and spread between centromeres, then was homogenized within each.

## How the unit spread and why it is kept (2026-10-23; `insert/`, `figure/composite_unit.png`, code 11, 12)
- **Insert copies** (`insert/`):
  - Mac: 69 loci with ≥ 3 kb of the insert, all centromeric, only on LG08 (24), LG25 (23) and LG16 (22). Median identity 91%.
  - Gar: no locus ≥ 3 kb.
  - Short fragments occur genome-wide in both species, most outside centromeres.
- **Insert domains:** no protein-domain hits (tblastn, e ≤ 1e-3); 4 ORFs ≥ 100 aa. Origin not identified.
- **Satellite variant follows the unit type** (MacSat03 chromosome-consensus identity):

  | Pair | Identity |
  |---|---|
  | type II vs type II | 95.1% (LG08–LG25 96.9, LG08–LG16 95.1, LG16–LG25 94.2) |
  | type I vs type I | 90.2% |
  | type I vs type II | 82.7% |
  | type II vs other centromeres | 82.1% |

  Whole units, satellite included, moved between the three type II centromeres.
- **Interpretation (L2):**
  1. In the Mac lineage, the RTE-like element became part of the MacSat03 repeat unit.
  2. Units spread between non-homologous centromeres (shared variants across chromosomes).
  3. One unit variant gained a ~6-kb insert replacing most of the element. This Mac-specific type II variant spread to LG08, LG25 and LG16.
  4. Within each centromere, units are homogenized (copies ≤ 1 Mb apart differ by 0.001–0.004 tree distance).
  - Maintenance is by homogenization of the whole repeat, not by retrotransposition: the element is silent in most RNA-seq samples, and its centromeric copies are local, same-chromosome relatives.
- **Figure:** `figure/composite_unit.png` / `.pdf`.
  - A: one unit per centromere, to scale.
  - B: element nodes per chromosome.
  - C: unit length vs satellite.
  - D: MacSat03 variant by type.
  - E: homogenization vs distance.
  - F: summary.
- **Monomer level** (`figure/composite_unit_depth.png`, code 14; 5,782 MacSat03 monomers inside units, 20 subfamilies at 90%):
  - All three type II centromeres carry the same monomer subfamily, Macsf065: 88% (LG08), 100% (LG25), 86% (LG16).
  - Each type I centromere has its own dominant subfamily: LG07 and LG01 Macsf182 (81%, 100%), LG04 Macsf102 (65%), LG02 Macsf103 (82%), LG23 mixed (Macsf536, Macsf8563).
  - So the type II unit carries its own satellite variant: it spread as a whole unit (satellite + truncated element + insert).
- **Element map:**
  - The element has 13 short ORFs (≥ 100 aa) and no long intact reading frame. It is degraded, consistent with no autonomous activity.
  - Type II units keep only about 1.0–1.7, 2.0–2.3, 5.5–6.9 and 7.2–7.5 kb of the element.
- Unit boundaries (start of each element core copy) are drawn as thin black lines in `composite_unit_context.png` panel A and `composite_unit_depth.png` panel B.

## k-mer divergence and a proposed history (2026-10-23; `kmer_tree/`, `figure/kmer_divergence.png`, code 15)
Method:
- Canonical 15-mers, Mash distance, neighbour joining.
- Unit tree: 177 whole units, midpoint root.
- Element-core diversity: median pairwise distance within each group, scaled by the median Mac-vs-Gar core distance.

Results:
- **Unit tree:**
  - Every centromere forms its own clade (within-centromere distance 0.002–0.013).
  - The deepest split separates type II (LG08, LG25, LG16) from type I.
  - Type I branching: LG04 first, then LG02, then LG23, then LG07 + LG01 (sisters; both carry the Macsf182 monomer subfamily).
  - Type II branching: LG16 first, then LG08 + LG25.
- **Between centromeres of the same type:**
  - type II–type II about 0.04–0.05;
  - type I–type I about 0.12;
  - type II vs type I about 0.20.
  - Type II centromeres share a much more recent common unit than the type I centromeres do.
- **Element-core diversity** (1 = Mac-vs-Gar):

  | Group | Relative diversity |
  |---|---|
  | Mac arm copies | 1.27 |
  | Gar copies (sample) | 0.58 |
  | Mac centromeric copies, per centromere | 0.008–0.055 |

Proposed order of events (L3, relative only; no clock):
1. An old RTE-like LINE is dispersed in the ancestor. The Mac arm copies are relicts, older than the Gar–Mac divergence of copies.
2. Gar: a later, separate dispersed burst (Gar copies are less diverse than Mac arm copies).
3. Mac: the element is captured in the MacSat03 unit; unit lineages diverge into type I and type II.
4. Type I spreads to LG04, LG02, LG23 and later LG07/LG01. Type II (with the insert) spreads more recently: first LG16, then LG08 and LG25.
5. In every centromere the units are recently homogenized (concerted evolution). That keeps the arrays young-looking and erases their individual age.

Caveat: within-centromere homogeneity measures the time since the last homogenization, not the age of the centromere.

## Normalized homogenization (2026-10-23; code 16; replaces panels B/C of kmer_divergence.png)
Normalization:
- Every centromere rarefied to 13 units (the smallest centromere, LG02).
- Every unit reduced to a fixed 2,000-k-mer MinHash sketch, which removes the unit-length effect (type II units are longer).
- 200 replicates. Pairwise distances use only the two units compared, so there is no shared denominator.

Homogenization index H = 1 − within / between (`figure/homogenization.png`):
- LG07 0.90, LG02 0.93, LG23 0.94, LG04 0.96, LG25 0.96, LG08 0.97, LG16 0.97, LG01 0.98.
- All centromeres are strongly homogenized. LG07 is the least (more unit variation), LG01 the most.
- Without normalization the values differ by at most 0.01.

Element copy diversity (`figure/element_age.png`; 60 copies per group, 500-k-mer sketches):

| Group | Diversity |
|---|---|
| Mac, outside centromeres | 0.150 |
| Mac, in centromeres (pooled; includes between-centromere differences) | 0.097 |
| Gar, all copies | 0.065 |

- The Mac dispersed copies are the oldest.
- The Gar copies come from a younger burst.

## Divergence part by part (2026-10-23; `figure/divergence_by_part.png`, `kmer_tree/divergence_by_part.tsv`, code 17)
Method: 11-mer Mash distance per part (satellite block, element in 500-bp windows, insert); up to 300 random pairs per comparison.

- **Within a centromere:** every part is near-identical, 0.000–0.006 along the element, 0.016 for the satellite, 0.004 for the insert. Homogenization covers the whole unit uniformly.
- **Between centromeres, same type:**
  - most conserved: the RT core (0.07–0.10 at 6–7.5 kb);
  - most diverged: the built-in satellite stretch (0.23 at 3–3.5 kb);
  - rest of the element 0.10–0.17;
  - MacSat03 block 0.155.
- **Type I vs type II:**
  - strongest difference at 1–2.5 kb (0.24–0.31), the region type II nearly lost;
  - RT core 0.08–0.15;
  - satellite 0.23.
- **Dispersed copies outside centromeres:** RT core 0.15–0.19. This is more diverged than centromeric copies even from different centromeres (0.07–0.10). Centromeric copies are younger and/or kept alike; dispersed copies are older relicts.
- **Insert:** 0.047 between the three type II centromeres, the least diverged part between centromeres. It is the most recently shared component.

Ranking of divergence between centromeres (high to low): built-in satellite > element body > MacSat03 block > RT core > insert.

## Anchored unit boundaries (2026-10-23; supersedes the core-start units of steps 09–19; code 20, `unit_anchored/`)
- **Why:** units cut at the start of each core hit did not start at the same homologous base. Partial hits and minus-strand copies shifted the start by a median 645 bp (> 1 kb for 12% of units). Downstream divergence compared non-homologous segments. Steps 09–19 are kept for the record but are superseded; the sensitivity and bootstrap runs (18, 19) were stopped and not used.
- **Anchor:** element consensus position 6,500, inside the RT core and kept by both unit types. Its genome position is interpolated inside every element hit covering it, strand-aware. A unit runs from one anchor to the next on the same strand (8–25 kb), oriented so the element reads forward.
- **Checks:**
  - Anchor 6,000 gives the same units (median length difference 0 bp, 90th pct 15 bp; n = 272).
  - Anchor 7,000 falls in a region most type I units have lost, so it is not usable (n = 18).
  - Strands are consistent along arrays.
- **Result:** 293 units on 12 centromeres, versus 177 on 8 before.
  - Four type I centromeres were missed before: LG15, LG19, LG21, LG26. Their copies cover the anchor but fall below the 70%-of-core threshold.
  - Type II is still exactly LG08, LG16, LG25: every unit there carries the insert.
  - Median unit length per centromere: 12.0–16.5 kb, sd 1.6–3.2 kb; ≤ 2 units per centromere are > 1.6× the median.

## 2026-10-23 correction: satellite block definition (steps 21-22)
- Bug: satblock.py excluded MacSat03 monomers overlapping any element BLAST hit. The element consensus contains its own MacSat03 stretch
  (query 2,928-3,343), so that stretch hits every satellite monomer and the rule deleted almost the whole block
  (median 1 monomer per unit kept vs ~18 present).
- Consequence: the first step-21 rerun satellite numbers (span ~0.01 between centromeres, divergence 0.16-0.21) are VOID.
  Element and insert numbers were not affected.
- Fix: hits whose query interval lies only inside 2,880-3,400 are ignored when masking. Block now median 18 monomers / ~4.2 kb per unit
  (LG23 lower: 7 monomers). Steps 21-22 rerun with the fix.
- Figure status after this fix: figure/homologous_divergence.png (22:35) was drawn from the buggy run -> STALE, will be redrawn by the rerun.
  sensitivity.png / bootstrap.png (steps 18-19) were never final (old boundaries). Other figures do not use the satellite block.

## 2026-10-24 correction: element windows (steps 21-22)
- Bug: element windows concatenated every BLAST hit overlapping the window. Satellite monomers hit the element's built-in MacSat03
  stretch (~20 hits per unit) and internal repeats overlap in query space, so windows near 2,900-3,400 became multi-kb chimeras
  (order set by a sort on sequence). Step 22 pools windows per unit pair, so it gave lag-1 element K2P 0.035 (step 21 per-window: 0.002)
  and placed Gar copies (0.119) closer than other Mac centromeres (0.154). Both numbers are VOID.
- Fix: satblock.window_hits(): drop hits confined to the built-in stretch (element consensus only), keep non-overlapping hits per
  window, best bit score first. Steps 21-22 rerun.
- Satellite-block fix from 2026-10-23 holds: first corrected run gave lag-1 satellite edges K2P 0.010-0.014 at span ~1.0,
  between centromeres 0.18-0.20 at span ~0.83, type I vs II 0.23-0.24 at span 0.85 (to be confirmed by this rerun).

## Step 22b (2026-10-24): matched windows, other Mac chromosome vs Gar copy
- Same windows for unit i vs another-chromosome unit j and vs a Gar copy g (median 1 window = 500 bp shared by all three; thin).
- Typical pair: Mac other chromosome 0.121, Gar 0.116 (per unit, then over units); Gar is closer or equal in 9 of 11 centromeres
  (LG15 has no window shared with the Gar sample). Nearest relative (min of 20): Mac 0.075 vs Gar 0.096; Gar is nearest for 58 of 250 units.
- Reading (L2 at most): each centromere's nearest relatives are other Mac centromeres, but a typical pair of Mac centromeres is about as
  divergent as Mac vs dispersed Gar copies. Not yet an age statement: the needed control is Mac DISPERSED copies on the same windows
  (if they sit at the same distance, this is the age of the element family, not of the species split). Output: homologous/matched_cd.tsv.

## Step 32 (2026-10-24): the source element (code 32*.py/.sh)
- **The transposon is a ~3.1-kb RTE LINE, not the 8,164-bp EDTA entry.**
  - Dispersed copies in both species share sharp ends at EDTA positions ~5,550 and ~7,820 (Gar arm: 891 of 1,134 start at 5,550-5,600)
    (`core/*.copy_ends.tsv`). The rest of the EDTA consensus (0-5,550, 7,820-8,164) exists only in Mac centromeric units:
    EDTA TE_00003409_INT is a consensus of the centromeric UNIT (element remnant + satellite + unit sequence).
  - Source consensus: 37 Gar full-length dispersed copies (MAFFT L-INS-i), 3,117 bp (`seq/Gar_full_consensus.fa`; Mac 13 copies: 3,014 bp,
    `seq/Mac_full_consensus.fa`). One ORF of 1,027 aa on the minus strand (`seq/Gar_full_orfs_minus.faa`).
  - RT at ORF aa 548-954, best match RTE-clade LINE RT (blastp E = 7e-88, 39%; next RTEX 2e-58). N-terminal 547 aa (where RTE carries APE EN)
    intact but untested (no EN reference). 3' end in LINE sense ...AAAGAAGAAGAAGA ((AAG)n tail). TSDs absent (2 of 37), as usual for RTE.
- **Lineage (core tree, rooted on the Gar 266-copy clade):** all Mac centromeric copies descend from one Mac-only lineage (241 tips; SH 94)
  whose base is a grade of 10 single Mac arm copies (`tree/spine_arm_copies.tsv`, `tree/spine_arm_context.tsv`). No arm copy is full length
  (0.9-4.2 kb). A cluster of ~20 near-identical 2.3-kb Mac arm copies is sister to type II. Source subfamily = Mac dispersed RTE lineage;
  Gar and Mac source consensuses are about equally close to the units (Mac +0.3-0.5% identity).
- **Truncation state per centromere** (all 293 units, `homologous/source_state_per_unit.tsv`, `source_bins_per_chrom.tsv`), constant within
  each centromere:
  - LG04: near-complete source (6-3,058) = ancestral state.
  - Type II (LG08/16/25): 3' 2.4-2.7 kb, 5'-truncated at ~2,740 (typical LINE 5' truncation).
  - LG02/21/23 and LG01/07: shared internal deletion ~1,700-2,700 (same left breakpoint ~1,699) = one shared derived event; LG01 = LG07 profile.
  - LG19: own deletion ~800-2,500. LG15/LG26: RT core only (~255-1,500).
- **CORRECTION** to "Two unit types" above: type II is not "truncated to the RT core". Measured on the source it keeps ~2.7 kb, more than
  type I LG01/02/07/21/23 (~1.5 kb). The EDTA consensus lacks source 1,700-2,680, so that part of type II was counted as "other".
  The insert itself contains no source sequence (blastn, no hit).
- First-pass junction check (`homologous/tail_junctions.tsv`): in LG02/04/19/21 the element 3' end lies within ~0-50 bp of a MacSat03
  monomer, with one identical flanking sequence across these centromeres; type II centromeres share a different one. Orientation of the
  junction extraction not yet verified - L1 lead only.

## Step 33 (2026-10-24): split arrays into validated satellite and TE spacers (`split/`, code 33*)
Method: 19 array regions (anchored-unit spans +-15 kb). Source RTE (3,117 bp) and type II insert mapped by dc-megablast. Every pipeline
MacSat03 monomer call aligned (local, rotation-aware) to the MacSat03 consensus; valid = score fraction >= 0.30. Valid monomers merged
into blocks (gaps < 300 bp, >= 3 monomers); spacer = interval between blocks. Array orientation = majority sequence strand of valid monomers.
- **Monomer layer was noisy:** score fractions are bimodal (4,493 calls < 0.1; 172 in 0.1-0.4; 6,207 >= 0.4). 4,600 of 13,100 MacSat03
  calls in these arrays (35%; LG04 59%, LG19 57%, LG02 52%) are not satellite: stringdecomposer forced "monomers" over spacer sequence
  (mixed strands, no similarity to neighbours). Threshold 0.30 is in the empty part of the distribution (0.1-0.4 changes <= 172 calls).
  -> Everything that used unvalidated monomers (satblock blocks in steps 21-22, monomer subfamily panels, earlier "tight spacer" calls) needs
  rerunning on `split/monomer_validation.tsv` (valid == True).
- **One architecture in all 12 centromeres** (`split/spacers.tsv`, `spacer_summary.tsv`): satellite block | spacer | satellite block, one TE
  fragment per spacer, always in the same orientation (RTE 3' tail toward the preceding block).
  - Type I spacers 7.8-9.8 kb: left flank 1.6-4.7 kb | RTE fragment 0.8-3.0 kb | right flank 3.6-5.5 kb.
  - Type II spacers 10.2-11.1 kb: insert ~5.8-6.1 kb | RTE fragment 2.3-2.6 kb | ~0.5 kb flank with MacSat03-like sequence.
  - RTE end coordinates are constant within a centromere (81-100% of spacers within 10 bp), see step 32 states.
  - Satellite phase (MacSat03 consensus position) at the block end before and the block start after the spacer is constant within a centromere
    (85-100% for most) and differs between centromeres (e.g. LG01 55/56, LG02 75/76, LG04 172/175, LG07 185/186, LG23 125/128).
    LG01, LG02, LG07, LG23: the spacer sits between two consecutive monomer positions = one clean insertion point, inherited by every unit.
- **Flanks** (`split/flank_composition.tsv`): not satellite and not the source RTE. Covered 70-96% by the EDTA unit consensus and by
  high-numbered satellite "families" (MacSat62xx) that were themselves built from these units (circular, not an identification).
  Type I right flanks: ~1 kb of short hits (76-82%) to EDTA TE_00004794_INT (5,918 bp, "LTR/unknown"), which has a Rex1-clade LINE RT
  (tblastn E = 3e-70) and BEL/Pao-like RT hits -> candidate second, older TE remnant in the unit (L1 lead, not verified).
  Type II left flank = the insert (60% covered by DNA/DTA hAT library hits, as before).

## Step 34 (2026-10-24): switch to SatChain v4.6 as the satellite layer (code 34a, 34b)
- v4.6 (copy-level HMM, runs trimmed to confident copies, per-run sequence check) contains 98% of the monomers validated in step 33 and 0% of
  the 4,584 false calls; all 276 step-33 TE spacers match a v4.6 interruption (median overlap 95%; v4.6 class TE_insertion 274/277).
  From here on: satellite = v4.6 runs (MacSat03 = CmaSat001-225), spacers = v4.6 interruptions_between_runs.
- RTE fragment per interruption by iterative local alignment to the source (segments >= 100 bp, score >= 150) (`split/v46_fragments.tsv`,
  `v46_fragment_summary.tsv`):
  - 15 MacSat03 centromeres carry the RTE in (nearly) every spacer: LG01, 02, 03, 04, 07, 14, 15, 17, 19, 21, 23, 26 (type I) and
    LG08, 16, 25 (type II). New vs the anchored-unit set: LG03, LG14, LG17.
  - 9 MacSat03 centromeres have regularly interrupted arrays WITHOUT the RTE: LG05 (51 spacers), 06 (25), 09 (16), 12 (55), 13 (29),
    18 (23), 22 (16), 24 (8), 27 (15) -> a second spacer type to identify.
  - Fragment identity to the source 85-93%; end coordinates constant within a centromere (mostly 90-100% of spacers within 10 bp).
  - Shared deletion 1,698-2,686 (identical breakpoints): LG01, LG07, LG02 (52%), LG03 (with two more). LG17 853-2,797; LG19 760-1,114 + 1,283-2,486.
  - v4.6 reverse_fraction orientation is flipped for some chromosomes relative to step 33 sequence orientation (te_strand '-' rows:
    LG02, 04, 14, 15, 08, 16); after flipping, every centromere has the same TE orientation. Use sequence orientation for layouts.

## Steps 35-36 (2026-10-24): homogenized core figure and stability score
- Figure: `F:/data/satellite_FISH/centromere_composition_2026-10-24/centromere_composition.png/.pdf` (code 35b build, 35c plot; data tables there).
  Assembly only (no reads, no ChIP; the HiFi/ONT FASTQs carry no MM/ML methylation tags). "Homogenized core", never "active".
  Core stability (W x k x tau grid, cores.tsv): Mac 24-27 of 27 centromeres have a core at W >= 20 kb for every k and tau; W = 10 kb fails
  (shorter than the ~13-kb unit). Gar 4-24 of 27 depending on tau -> Gar cores not robust with this measure (h favours long-unit arrays).
- Stability score (`split/stability_scores.tsv`, edge-free check `stability_scores_noedge.tsv`; code 36b/36c), 24 MacSat03 centromeres:
  - copy confidence 95-100% confident (SatChain copies are reliable);
  - boundary stability: one junction usually exact, the other drifts (e.g. LG25 right 1.00 / left 0.05 at +-10 bp; LG01 0.85/0.91);
  - positional identity: leave-one-out position accuracy 0.32-1.00 with edge copies, 0.21-0.94 without edge copies and with only
    full-length (180-260 bp) copies (chance 0.14); permutation p <= 0.025 for 47 of 48 array ends (LG08 left p = 0.11).
    -> copies have a positional identity inside the unit (HOR-like); arrays resolve at copy level, not only unit level.

## Step 37 (2026-10-24): RTE structure and its 5' end (`figure/rte_structure.png/.pdf`, code 37)
- Map (LINE sense): one ORF 1,027 aa open to the 5' end of the consensus; first ATG at LINE 91 (double ATG, same site in the Mac consensus)
  -> ORF from Met = 997 aa (C. elegans RTE-1: 1,024 aa). RT aa 548-954 (RTE clade, E = 7e-88). APE-like motif 'CVQET' at aa 122 in the
  N-terminal region (sequence motif only; endonuclease domain not tested). 3' end ...AAGAAGAAGA ((AAG)n tail). No TSDs (2/37).
- 5' end: 19/37 Gar copies start within +25 +- 5 bp of the alignment boundary, 10 start 15-25 bp upstream; the Mac consensus extends 12 bp
  further with the same sequence -> true 5' terminus ~12-25 bp upstream of Gar consensus position 1. Upstream = AT-rich, non-conserved insertion
  sites (majority base share ~0.4). 5' UTR ~90-115 bp; no CCAAT / GC box, the only TATA lies inside the ORF; GC of the first 100 bp = element mean.
  Promoter not identified (L1; RTE promoters are generally internal or unmapped).
- Centromeric copies lose 5' material: type II lacks LINE 1-~380 (start codon and the start of the N-terminal region); the shared type I
  deletion removes LINE ~432-1,420 (N-terminal / expected endonuclease region, just 3' of the APE-like motif); LG04 keeps LINE 60 onward.

## Step 38 (2026-10-24): LG04 composite array dissected (figure/diagnostics/, split/spacer_decomp/)
- Scale: Mac_SAT00318, 371 kb; 28 satellite blocks = 599 copies (132 kb), 27 interruptions (239 kb); unit median 13.0 kb.
- Monomer classes checked against SatChain's own alphabet (step 65, Louvain on copy identity, copies inside runs, identity >= 80):
  ARI vs step-03 classes 0.79 (t 0.90, 8 communities), 0.84 (t 0.94, 14), 0.50 (t 0.97, 28: C01 splits into ~6 sub-variants).
  Fixed start robust at all thresholds (position 1 = one community in 96 % of runs; position 2 93-96 %); fixed end weaker (57-64 %).
  -> use the swept step-65 alphabet for topology; step-03 classes only as cross-check. Step-53 merges: all tested pairs separable
  (max_diag 34-365 vs threshold 5); thresholds not swept.
- Unit order (graph + position rasters): INT -> C06 (197 bp, half monomer + junction sequence) -> U (truncated C02) -> C03 -> C04 ->
  C02/C01 core (C02 -> C01 0.91-0.95) -> C08 -> C09 -> INT. Only the core length varies (6 to 50+ copies).
- Interruptions decomposed (medoid int15, 8,532 bp, 23 segments at feature boundaries + shared breakpoints): all 23 segments present
  in 24/27 interruptions at 96-100 %. Variants: int04 (22.8 kb) = segments S8-S22 (RTE + downstream flank, ~4.3 kb) in 4 tandem copies;
  a 245-bp insertion (265-bp version) just 5' of the RTE in int03, 05-07, 19-20, 24-27 (not satellite, not RTE; unidentified);
  int24 deletion S15-S19; int25/26 shared small S7 deletion; int01 = 1.3-kb truncated spacer.
- Spacer content: RTE 36 %, independent TE pieces 15 % (same positions in every unit -> present in the ancestral unit),
  sequence matched only by EDTA's TE_00003409 (built from this unit) 47 %, unannotated 1 %.
- RTE junctions (rte_junctions_LG04.tsv): same element (source 6-3058, 3,022 bp) at the same breakpoints in 25/26 interruptions;
  both 40-bp flanks identical in 25-26/26; a 15-bp direct repeat ATATAATTACAAGCA flanks the element (right copy adjacent to the 3' end,
  left copy 17 bp before the aligned 5' end) = candidate target-site duplication of ONE insertion, copied with the unit (L2).
  RTE orientation in LG04 (genomic): 5' -> 3' towards higher coordinates (the earlier spacer figure had inverted arrows; fixed).
- Across centromeres (exact matches only): the 15-bp repeat in LG04 72, LG02 40, LG19 27, LG23 23, LG21 20; absent from type II and LG01/07/15/26.
  Needs fuzzy junction alignment across all centromeres before a shared-insertion claim (and outgroup polarity before any origin claim).
- Step 38b, k-mer resolution of the interruptions (spacer_kmer_hires.py; figure spacer_kmer_hires_Mac_SAT00318_LG04.png): presence of every
  reference k-mer in each interruption; informative variant sites (absent in >= 2 and present in >= 2) = 89 / 88 / 84 / 76 / 46 at
  k = 15 / 21 / 31 / 51 / 101. Average-linkage spacer haplotypes: K = 2 and K = 3 partitions identical at every k (ARI 1.0); K = 4 ARI 0.95 at k = 101.
  H1 (with the 245/265-bp insertion) = interruptions 03-07, 19-20, 24-27; H0 = 02, 10-12, 14, 16-18, 22-23; H2 = 08, 09, 13, 15, 21.
  H1 mostly at both array ends, H0/H2 in the middle: resembles layered expansion (older edges, younger centre) - L2 hypothesis, one array,
  to test against satellite sub-variants (step-65 alphabet at 0.97) and in other centromeres. Candidate SatChain step: spacer variants + haplotypes.
- Step 38c, corrections to 38b: the spacer groups were relabelled A (insertion carriers, H1), B (H0), C (H2). Checked by diagnostic sites (k = 21):
  A vs rest 12 fully diagnostic sites (16 at 80/20); B vs C 9 fully diagnostic (14 at 80/20). C contains the reference int15 and matches it:
  B carries 9 extra shared differences (38b had it reversed) -> part of C is reference-likeness (reference bias; swap test pending).
  A is internally heterogeneous (mean pairwise k-mer difference 0.070 vs B 0.019). ~85 informative sites, 263-522 variant sites in total.
- Step 38d, in-register consistency (inregister_consistency.py; *_inregister_consistency_*.png/.tsv): unit = satellite block + interruption;
  satellite distance from in-register edge copies (first/last E copies), spacer distance = reference-free k-mer Jaccard; 25 units.
  Block BEFORE the interruption: Mantel r 0.42 / 0.47 / 0.52 (E 4 / 6 / 8, p <= 0.002); partial r (array distance removed) 0.28 / 0.40 / 0.46 (p <= 0.008).
  Block AFTER: r 0.15-0.30, partial 0.10-0.29. Same at k = 21 and 31.
  -> satellite and spacer of a unit share a history beyond proximity; the copied unit is "block then spacer" (phase supported by data).
  -> spacer group A is NOT a satellite-side group (within-A satellite distance 0.092-0.100 > A vs others 0.070-0.080): the insertion and its
     variants behave like spacer-specific changes, not a separate lineage of whole units. The layered reading is not supported for A. (L2)
