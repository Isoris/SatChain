"""Colours for the LTR / MacSat03 figures (2026-10-27 colour audit). Everything is taken from palette_v46 (runs_of_record_v4.6_2026-10-27/code) and class_colours.
Meaning (same in every figure of this folder):
  SAT      = MacSat03 satellite (the family colour of the reference figure, crimson)
  ELEM     = the RTE-like element (a LINE: TE LINE amber);  ELEM_CORE = its RT core (dark amber, = RTE_EDGE of class_colours)
  INS      = the ~6 kb type II insert (slate, no other meaning anywhere)
  TYPE1/2  = unit type I (complete element) = SAT crimson, type II (element cut + insert) = INS slate
  REL      = ordinal ramp for comparison levels (closest relatives darkest -> most distant lightest); never a family / class meaning
  CHROM    = one muted colour per chromosome (class_colours.CLASS_SEQ order), only where chromosomes must be told apart
  CENTROMERE / OTHER / unit-boundary = ink / light grey / ink"""
import sys
sys.path.insert(0, "F:/data/satellite_FISH/runs_of_record_v4.6_2026-10-27/code"); sys.path.insert(0, "F:/data/satellite_FISH/LTR_MacSat03_2026-10-22/figure/diagnostics")
import palette_v46 as P
SAT = P.FAMCOL["MacSat03"]; ELEM = P.TE["LINE"]; ELEM_CORE = "#7a5200"; INS = "#3d4f66"; OTHER = "#d9d9d9"; INK = P.INK
TYPE1, TYPE2 = SAT, INS
REL = ["#5b5b5b", "#7a7a7a", "#9a9a9a", "#bdbdbd", "#d9d9d9"]   # grey ramp, closest = darkest but never black
CHROM_SEQ = ["#1f4e9c", "#1b9e8a", "#7b4fa3", "#c0392b", "#7f9f2f", "#4a90d9", "#b5468a", "#2e8b3d"]   # = class_colours.CLASS_SEQ[:8] (copied: class_colours imports Bio and is owned by another agent)
def type_col(is_type1): return TYPE1 if is_type1 else TYPE2
