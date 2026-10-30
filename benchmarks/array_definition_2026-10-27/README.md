# What is an array? Step-02 interval (A) vs SatChain run chain (B) — 2026-10-27
Code: compare_array_methods.py; tables: summary.tsv, centromeres.tsv, compare_<sp>_<G>.tsv.
- A (inherited from v4.5 step 02): 10-kb windows / 5-kb step, lag identity >= 0.34 at lag 40-1000 bp, positive windows merged when <= 20 kb apart,
  edges refined to the TRF same-period union. Defined BEFORE any copy is decomposed; STEP_METHOD calls it "search space only".
- B (SatChain): neighbouring runs of record (56 + val56) of the same family, interruption <= G, never across an assembly gap. G swept 10/20/50 kb.
Results (multi-run A arrays; G = 20 kb): Mac 738 -> 518 identical, 207 split by B, 8 merged, 5 complex; Gar 570 -> 450 / 107 / 7 / 6.
Why B splits A: family change inside A (Mac 293 breaks, 156 A arrays mixed-family; Gar 109 / 83), assembly gap (Mac 43 breaks, 53 A arrays contain a
gap; Gar 6 / 39), interruption > G (Mac 234 / 67 / 6 at 10 / 20 / 50 kb; Gar 258 / 21 / 3).
Major centromere identical under A and B: Mac 24/27, Gar 22/27. Differences: Mac LG04 (int04 = 22.8-kb tandem spacer > 20 kb), Mac LG05 (2 assembly
gaps), Mac LG20 (MacSat1271/1272 alternate; 316- and 579-kb interruptions inside one step-02 interval), Gar LG12 / LG23 (assembly gap), Gar LG19
(GarSat285 | GarSat09 | GarSat285), Gar LG26 (GarSat115 -> GarSat03), Gar LG18 (A's runs interleaved with another array's run).
Lesson: a fixed interruption length is a weak rule (LG04's 22.8-kb interruption is 4 spacer copies, i.e. part of the array); interruptions should
be judged by content (spacer-like or not), and composition changes inside runs need the domain layer (domain_seg.py).

## PNG audit (2026-10-27; png_audit.py -> png_audit.tsv, png_audit_unreviewed_sheet.png)
Rules = what replaced what and why, plus "older than the v4.6 catalogue". Summary: current 180, duplicate 108 (centromere *_linear / *_topology =
same as *_rlvl3_*), outdated 63, superseded 22, needs review 7. Outdated: gallery pages (copy barcodes), composite_unit* / anchored_units (EDTA-unit
based), centromere_composition and residuals (EDTA / AnnoSINE calls), spacers / spacer_decomposition (EDTA), G11f/g primer slides (old names and
verdicts), 16 defense-slide figures older than v4.6, motif_views_v2 Gar LG28, rlvl1 overviews (strand calls corrected). Needs review: the 10-23 RTE /
unit divergence figures (unit boundaries pre run-based split). Strand table recomputed 2026-10-27 by k-mer sharing with the family representative
array: 5 calls changed (Mac LG24; Gar LG03, LG04, LG18, LG28).
