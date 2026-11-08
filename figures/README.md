# figures/ (added 2026-11-01)

The figure code of the fClaHyb satellite work, moved into the repository because every picture is drawn from SatChain tables:

| folder | origin (read-only copies stay there) | what it draws |
|---|---|---|
| `paper_figures/` | `satellite_FISH/paper_figures_2026-10-27` | Figure 1 (satellite organisation) and Figure 2 (the RTE LINE in the centromeric arrays): one module per panel (`panels/figN_panel_X.py`, `draw()` contract in `LAYOUT.md`), `code/assemble_v2.py` builds the 183-mm figures at 500 dpi; captions and per-panel number tables beside the modules. Upgraded to v4.7 on 2026-11-01 (runs of record `satchain_v4.7_2026-11-01/runs_of_record`, "lineage" wording replaced). |
| `ideogram/` | `runs_of_record_v4.6_2026-10-27/code` | chromosome ideograms with satellite families and candidate centromeric arrays (`ideogram_share_v46.py`, palette) |
| `centromere_composition/` | `satellite_FISH/centromere_composition_2026-10-24` | composition of every candidate centromeric array (families, strand, RTE / TE insertions) |
| `LTR_MacSat03/` | `satellite_FISH/LTR_MacSat03_2026-10-22/code` | the MacSat03 composite unit: RTE LINE structure, spacer types, k-mer trees, topology views |
| `../docs/img/` | `satchain/plot53_*.py`, `plot53b_units.py`, `kmer_separation_sweep.py` | README showcase figures (array topology, profiles, k-mer separation) |

Inputs are the runs of record and catalogue of a SatChain run folder plus the annotation folders (`docs/INPUTS.md`); each script names its
inputs at the top. Project rules: 450-500 dpi PNG + vector PDF, every colour labelled, no "centromere" without "candidate".
