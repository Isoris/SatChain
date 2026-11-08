# Figure 1a caption stub (evidence level L1: sequence) - SatChain v4.7 wording, 2026-11-01

**Trees of the two satellite groups, with four levels of resolution.** Maximum-likelihood trees of the v4.5 monomer copies (IQ-TREE 2.3.6, GTR+G, 1,000 ultrafast bootstraps; treefiles
`share_LITE_v4.5_2026-10-15/phylogeny/{lineageB_5S_derived,lineageA}/*.treefile`; the file stems are the v4.5 input names, the word "lineage" is no longer used):
the 5S-related tree, 256 copies (GarSat04 = the 5S-derived satellite, SatChain v4.7 sequence cluster C2.2; the 5S rDNA gene arrays GarSat49 and MacSat06, which are gene arrays
and not satellite families; + 16 5S rRNA gene references); the centromeric tree, 624 copies (GarSat03, GarSat09, GarSat10, MacSat03, MacSat162 = the GarSat03 / MacSat03 group,
the RTE-carrying centromeric families, sequence cluster C2.1); up to 12 sampled copies per chromosome array. Unrooted trees, equal-daylight layout, one scale bar each (substitutions/site);
the two groups are not joined (no joint alignment; between-group consensus identity 53-58 %). Tips: fill = family, circle = *C. gariepinus*, square = *C. macrocephalus*, diamond = 5S rRNA gene.
(1) **Sub-clades** R1-R19 (centromeric tree) and S1-S6 (5S-related tree), numbered by size: maximal sub-trees with patristic diameter <= 0.75 (S) / 0.7 (R)
substitutions/site, UFBoot >= 95 and >= 5 copies (the 5S reference group is never a clade); 218/256 (S) and 619/624 (R) copies fall in a clade; label = name + copies,
or name + chromosome for single-array clades (all sampled copies of that array).
(2) **Chromosome origin**: ring around each clade, one segment per chromosome array (segment width = copies); hue = species (purple *C. gariepinus*, green *C. macrocephalus*); dark = array flagged
centromeric in the v4.7 runs of record (`centromeric_L1`), pale = other arrays, open = 5S rRNA gene loci (not arrays); two hues in one ring = clade with both species. (3) **Support**: UFBoot on the
backbone nodes (nodes joining >= 2 clades): filled >= 95, open 70-94, none < 70; the gold ring marks where the 5S rRNA gene group joins the 5S-related clades (UFBoot 96).

## Clade rule and sweep (parameter stability; table: `fig1_panel_a_sweep.tsv`)
Cut swept 0.10-1.20 in steps of 0.05 (23 values), clades per cut, S: 0.1:4, 0.15:8, 0.2:9, 0.25:11, 0.3:11, 0.35:12, 0.4:13, 0.45:6, 0.5:6, 0.55:6, 0.6:6, 0.65:6, 0.7:6, 0.75:6, 0.8:6, 0.85:6, 0.9:6, 0.95:6, 1:6, 1.05:6, 1.1:6, 1.15:6, 1.2:6. R: 0.1:30, 0.15:38, 0.2:37, 0.25:39, 0.3:42, 0.35:38, 0.4:39, 0.45:37, 0.5:26, 0.55:26, 0.6:26, 0.65:19, 0.7:19, 0.75:19, 0.8:11, 0.85:9, 0.9:9, 0.95:1, 1:1, 1.05:1, 1.1:1, 1.15:1, 1.2:1. The drawn cut is the middle of the longest run of IDENTICAL clade sets (non-trivial partitions only: >= 4 clades,
none holding > 60 % of assigned copies; ties: more copies assigned, then coarser): S 0.75 (plateau 0.45-1.1, 14 consecutive cuts, one clade set), R 0.7
(plateau 0.65-0.75, 3 cuts). Stability = share of assigned copies that keep their clade partners (best Jaccard >= 0.8 with a clade of the other cut):
S 1.00-1.00 over cuts 0.45-1.10; R 1.00-1.00 at the plateau, 0.63 at cut 0.60 and 0.40 at 0.80: the 5S-related clade set is
robust, whereas the 19-clade level of the centromeric tree is a narrow plateau (cut 0.60 gives 26 clades, 0.80 gives 11).
Clades present (Jaccard >= 0.8) at >= half of the cuts 0.30-1.00: S 6/6, R 15/19 (column `persistence_cuts_030_100` of the numbers table); the cut-dependent clades of the centromeric tree are R1, R2, R6, R8 (large multi-chromosome clades that merge / split with the cut).
UFBoot minimum swept at the drawn cut (clades): S >=70: 7, >=80: 7, >=90: 6, >=95: 6, >=99: 6; R >=70: 12, >=80: 12, >=90: 15, >=95: 19, >=99: 22. The clade set of the centromeric tree depends on the support threshold (95 = conventional); the trees carry UFBoot only (IQ-TREE was run with -B 1000, no SH-aLRT).

## What the levels show
* Species mixing (share of copies whose nearest patristic neighbour on ANOTHER chromosome array belongs to the other species; `nn_cross_species`): over all satellite copies 0.000 (S) and 0.000 (R), against 0.479 / 0.511 expected if species were unrelated to tree position.
  1 clade(s) contain both species (S3; share inside the clade 0.0, expected 0.667): the two species form separate sub-groups, they do not interleave at copy level (species-specific expansion rather than a shared pool of monomers).
* Chromosome origin: the centromeric-tree clades hold 34 centromeric and 14 other arrays; 4 array(s) without a covering v4.7 run (the v4.7 tables have no rows for Mac ; Gar , so copies on those chromosomes are drawn pale = unknown); clades made only of non-centromeric arrays: R13, R15, R16, R18, R19.
* 5S: the 16 5S rRNA genes (8 Gar on LG04, 8 Mac on LG01 + LG20) form one clade (UFBoot 100, monophyletic = True) on a long branch (0.30 substitutions/site); the lowest node that joins them to a satellite clade
  (UFBoot 96, 31 copies: the 5S genes + MacSat06 copies incl. S6) lies inside the MacSat06 part of the tree (MacSat06 = the *C. macrocephalus* 5S rDNA array, 72 % of its units carry an annotated 5S gene in v4.7);
  the nodes on the way there have UFBoot 91/73/72/46/89/96, so the 5S-satellite relationship is supported at the last node only; the direction of the derivation is not resolved because the trees are unrooted.
* Copies of the centromeric tree that are not in a v4.7 run of their own family take the centromeric flag of the run covering the copy, where one exists.

Sources: treefiles above; `SatChain v4.7 runs of record (satchain_v4.7_2026-11-01)`: `runs_of_record/{Gar,Mac}/runs_of_record.tsv` (centromeric_L1), `satellite_family_nomenclature.tsv` (element_class, seq_cluster);
`fig1_panel_a_numbers.tsv` (family / clade / backbone / 5S rows), `fig1_panel_a_sweep.tsv`; code `fig1_panel_a.py` (`fig1_panel_a_combined.py` = same content on one common scale).
