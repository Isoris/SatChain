"""Audit of figure PNGs: which are current, superseded or outdated (2026-10-27).
Scans the satellite figure folders, SatChain benchmarks, motif views and slide figures. Status from explicit rules (what replaced what and why)
plus a date rule: a PNG older than the v4.6 catalogue (fClaHyb_Mac.satDNA_v46.catalogue.tsv) predates the current runs and array definition.
Output: png_audit.tsv (path, modified, status, reason), summary printed.
usage: python -P png_audit.py"""
import os, re, glob, datetime, pandas as pd
ROOTS = ["F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/figure", "F:/data/SatChain/benchmarks", "F:/data/satellite_FISH/motif_views", "F:/data/satellite_FISH/motif_views_v2",
         "F:/data/thesis_defense_slides_tpl/figures", "F:/data/satellite_FISH/centromere_composition_2026-10-24"]
CAT = os.path.getmtime("F:/data/fClaHyb_func_annotation_v1.0/fClaHyb_Mac/fClaHyb_Mac_satelliteDNA/fClaHyb_Mac.satDNA_v46.catalogue.tsv")
STRAND = os.path.getmtime("F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/figure/diagnostics/centromeres/array_strand.tsv")
CHANGED_STRAND = ("LG24_Mac_SAT01699", "LG03_Gar_SAT00231", "LG04_Gar_SAT00280", "LG18_Gar_SAT01075", "LG28_Gar_SAT01588")
RULES = [   # (regex on path, status, reason) - first match wins
 (r"ideogram_(Mac|Gar)_rlvl[123]\.png", "current", "nested levels (k-mer recurrence), legend column"),
 (r"kmer_homogeneity\.png", "current, caveat", "runs only: degraded edge copies / flanks trimmed by the runs are not counted"),
 (r"inregister_consistency_v2", "current", "alignment-based spacer distance, both pairings"),
 (r"topology_order_profiles_Mac_SAT00318_LG04", "current", "LG04 topology v3 grouped by the following spacer (A / B / C)"),
 (r"topology_clean_", "superseded", "topology v1 (single INT node) -> topology v3 / rlvl3_topology"),
 (r"diagnostics/diag_Mac_SAT", "superseded", "step-03 class chain diagnostics -> rlvl2 stack + rlvl3 linear / topology"),
 (r"figure/(composite_unit|composite_unit_context|composite_unit_depth|anchored_units)\.png", "outdated", "composite unit drawn from the EDTA unit entry before the run-based split (user: unit not correct)"),
 (r"figure/(divergence_by_part|element_age|homogenization|homologous_divergence|kmer_divergence|per_unit_divergence|per_unit_kmer_grid)\.png", "needs review",
  "RTE / unit divergence on 10-23 element and unit copies (pre run-based split, pre strand fix); numbers may hold, unit boundaries do not"),
 (r"figure/rte_structure\.png", "current", "RTE source element structure (sequence-based)"),
 (r"centromere_composition\.png", "outdated", "family composition with EDTA TE insertions; replace with rlvl1 overview + residual families"),
 (r"inregister_consistency_Mac_SAT00318", "superseded", "v1 in-register test (k-mer Jaccard dominated by ins / RTE deletion / GT runs; blank rows) -> inregister_consistency_v2"),
 (r"diagnostics/topology_(Mac_SAT00318_LG04|clean|Mac_SAT00318_LG04_runs_INT)\.png", "superseded", "topology v1/v2 (single INT node, arbitrary halves, false transitions) -> topology v3 / centromeres rlvl3_topology"),
 (r"spacer_kmer_hires", "superseded", "reference-biased H0/H1/H2 spacer groups -> spacer_types_msa / spacer_int_view (reference-free)"),
 (r"diagnostics/spacers_Mac_SAT00318", "outdated", "draws EDTA calls (inconsistent between identical spacers)"),
 (r"spacer_decomposition", "outdated", "EDTA-cut segments S1-S23 -> sequence segments P1-P15 (spacer_content)"),
 (r"spacer_content_Mac_SAT00318", "needs update", "P-segment profile current; add curated module / microsatellite calls, the 245-bp ins and structural variants"),
 (r"spacer_types_msa_LG04", "current (supplement)", "reference-free variant matrix"),
 (r"spacer_int_view_LG04", "current (diagnostic)", "row view with flanking monomers; publication redesign pending (per-type architecture)"),
 (r"centromeres/.+_rlvl3_(linear|topology)\.png$", "current", "rlvl3 copy-level views"),
 (r"centromeres/.+_rlvl2_stack\.png$", "current", "rlvl2 stack with JSD domains"),
 (r"centromeres/.+_(linear|topology)\.png$", "duplicate", "same content as _rlvl3_linear / _rlvl3_topology"),
 (r"centromeres/.+_stack\.png$", "check", "prototype stack (Gar LG19/LG20) superseded by _rlvl2_stack"),
 (r"ideogram_types_(Mac|Gar)(_rlvl[123])?\.png", "superseded", "first type scheme (catalogue unit mixed micro / mini / satellite) -> ideogram_<sp>_rlvl1-3 (nested, k-mer recurrence)"),
 (r"ideogram_intervals_", "current", "step-02 interval vs run-chain comparison (G = 20 kb)"),
 (r"gallery_(Mac|Gar)_p", "outdated", "copy barcodes empty for NEW arrays and meaningless for micro / mini; replace with sequence-based (divergence) thumbnails"),
 (r"divmap_Gar_LG20", "outdated", "drawn with the v1 block rule (>= 0.8 blocks); redraw with the insulation-plateau call"),
 (r"divmap_", "current", "insulation-plateau array call (v2)"),
 (r"residuals\.png", "outdated", "relies on EDTA / AnnoSINE calls (judged unreliable); replace with residual families + manual curation"),
 (r"cross_Mac_LG04_x_LG24", "current", "strand-corrected"),
 (r"curate_spacers_", "current", "microsatellites masked"),
 (r"motif_views/", "superseded", "motif views v1 -> motif_views_v2"),
 (r"motif_views_v2/Gar_LG28", "outdated", "only 33 % of the plotted copies lie inside v4.6 runs (pre-extent-fix interval)"),
 (r"motif_views_v2/", "current, caveat", "haplotype strips untested (reference-based clustering)"),
 (r"slide_G11[fg]_primers", "outdated", "old family names, 14 gene arrays shown as satellites, v1 one-orientation in-silico PCR verdicts (verify_2026-10-24)"),
 (r"linear_states_Mac_SAT00318_LG04", "current", "A / B / C types, sequence features, no EDTA"),
 (r"genome_divmap", "current, caveat", "0.9 identity / 5-kb windows not yet swept; satellite and spacer windows mixed"),
]
rows = []
for root in ROOTS:
    for p in glob.glob(os.path.join(root, "**", "*.png"), recursive=True):
        q = p.replace("\\", "/"); mt = os.path.getmtime(p); st, why = "", ""
        for rx, s, w in RULES:
            if re.search(rx, q): st, why = s, w; break
        if not st:
            if mt < CAT: st, why = "outdated", "predates the v4.6 catalogue (pre-runs-of-record arrays / copies)"
            else: st, why = "unreviewed", "no rule; newer than the v4.6 catalogue"
        # per-array stacks / linear / topology orient copies from the copy table, not from array_strand.tsv: unaffected by the 2026-10-27 strand fix
        if q.endswith("centromeres_rlvl1_overview_Mac.png") or q.endswith("centromeres_rlvl1_overview_Gar.png"):
            if mt < STRAND: st, why = "outdated", "prints strand calls; 5 corrected 2026-10-27"
        rows.append(dict(path=q, modified=datetime.datetime.fromtimestamp(mt).strftime("%Y-%m-%d %H:%M"), status=st, reason=why))
A = pd.DataFrame(rows).sort_values(["status", "path"]); A.to_csv("F:/data/SatChain/benchmarks/array_definition_2026-10-27/png_audit.tsv", sep="\t", index=False)
A["folder"] = A.path.str.replace("F:/data/", "").str.rsplit("/", n=1).str[0]
print(A.status.value_counts().to_string()); print()
pd.set_option("display.width", 260); pd.set_option("display.max_colwidth", 120)
print(A.groupby(["folder", "status"]).size().unstack(fill_value=0).to_string()); print()
for s_ in ("outdated", "needs update", "check", "unreviewed"):
    x = A[A.status == s_].groupby("reason").path.apply(lambda v: f"{len(v)}: " + ", ".join(sorted(os.path.basename(t) for t in v)[:6]) + (" ..." if len(v) > 6 else ""))
    if len(x): print(f"== {s_}"); [print("  -", r, "|", v) for r, v in x.items()]; print()