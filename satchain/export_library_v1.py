#!/usr/bin/env python3
"""export_library_v1.py  (SatChain v4.7, 2026-11-09; satellite reference library in the 'library release format v1' shape)

Purpose: a portable family reference set for the deposit (consensus FASTA + catalogue + aliases + model status + checksums), built ONLY
from existing v4.7 outputs. No new measurement and no new identifiers (decision 2026-11-09).
Declared deviations from the format proposal (SatChain_20_Worklanes v4, LIBRARY_RELEASE_FORMAT_v1.md):
  - family_id = the deposited GarSatNN / MacSatNN / NewSat<H>nnnnn ids, not SAT000001 accessions;
  - no profile HMMs: SatChain_profiles.hmm.gz is not written, every family has model_status = not_attempted;
  - consensus = the 59g unit exactly as step 03 wrote it (dominant-class consensus of the family's largest run); no rotation or reverse
    complement is applied (normalization_method = none).
Representative haplotype: the genome of the id series (GarSat / NewSatGar -> Gar, MacSat / NewSatMac -> Mac), else the genome that has it.
Curated FASTA: families whose representative row has a unit consensus and 59h low_copy_flag = no. Others -> SatChain_candidates.tsv and
optional/candidates_consensus.fa. Nothing is manually curated: curation_status = computational_draft for all.
usage: export_library_v1.py <final dir> <annotation root> <dest dir> [--version 1.0.0]"""
import os, sys, csv, json, gzip, hashlib, argparse, datetime, collections
csv.field_size_limit(1 << 30)
ap = argparse.ArgumentParser(); ap.add_argument('final'); ap.add_argument('ann'); ap.add_argument('dest'); ap.add_argument('--version', default='1.0.0')
a = ap.parse_args()
REL = os.path.join(a.dest, f'SatChain_satellite_library_v{a.version}')
if os.path.exists(REL): sys.exit(f'{REL} exists: releases are immutable, issue a new version')
os.makedirs(os.path.join(REL, 'optional'))
def tsv(p): return list(csv.DictReader(open(p, encoding='utf-8'), delimiter='\t'))
def sha(b): return hashlib.sha256(b).hexdigest()
def fsha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for blk in iter(lambda: f.read(1 << 24), b''): h.update(blk)
    return h.hexdigest()
TAXON = {'Gar': 'Clarias gariepinus', 'Mac': 'Clarias macrocephalus'}
ASM = {'Gar': 'GCA_048129225.1 (fClaHyb_Gar, fClaHyb v1.0)', 'Mac': 'GCA_048129245.1 (fClaHyb_Mac, fClaHyb v1.0)'}
U = {G: {r['family']: r for r in tsv(os.path.join(a.final, 'unit_content_2026-11-07', f'59g_unit_content_{G}.tsv'))} for G in ('Gar', 'Mac')}
M = {G: {r['family']: r for r in tsv(os.path.join(a.final, 'family_metrics_2026-11-07', f'59h_family_metrics_{G}.tsv'))} for G in ('Gar', 'Mac')}
CAT = {G: tsv(os.path.join(a.final, 'cat59', f'59_catalogue_{G}.tsv')) for G in ('Gar', 'Mac')}
NOM = tsv(os.path.join(a.final, 'runs_of_record', 'satellite_family_nomenclature.tsv'))
HAND = {r['family']: r for r in tsv(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'docs', 'decisions', 'hand_decisions.tsv'))}
RUNAX = {r['family']: r for r in tsv(os.path.join(a.final, 'run_axis_2026-11-09', 'run_axis_family.tsv'))} if os.path.exists(os.path.join(a.final, 'run_axis_2026-11-09', 'run_axis_family.tsv')) else {}

def rep_genome(f):
    pref = 'Gar' if f.startswith(('GarSat', 'NewSatGar')) else 'Mac' if f.startswith(('MacSat', 'NewSatMac')) else None
    if pref and f in U[pref]: return pref
    return 'Gar' if f in U['Gar'] else 'Mac'
big = {}
for G in ('Gar', 'Mac'):
    for r in CAT[G]:
        f = r['family_v46']
        if (G, f) not in big or int(r['bp']) > int(big[(G, f)]['bp']): big[(G, f)] = r
aliases = collections.defaultdict(list)   # family -> [(alias, status, rationale, source)]
for r in NOM:
    if r['name'] and r['name'] != r['family']:
        aliases[r['family']].append((r['name'], 'retired_do_not_use', f'v4.6 / pipeline label in {r["haplotype"]}; retired 2026-10-29 naming audit', f'satellite_family_nomenclature.tsv:{r["haplotype"]}'))
for G in ('Gar', 'Mac'):
    seen = set()
    for r in CAT[G]:
        if r['family'] != r['family_v46'] and (r['family'], r['family_v46']) not in seen:
            seen.add((r['family'], r['family_v46']))
            aliases[r['family_v46']].append((r['family'], 'internal_provisional', f'pre-naming id of step 59 in {G} (column family_provisional_59 of the deposited catalogue)', f'59_catalogue_{G}.tsv'))

def rrna(u):
    cl, ax = u.get('gene_best_class', ''), u.get('sequence_axis', '')
    if cl in ('5S', 'rRNA_other') and ax in ('intact_gene', 'gene_derived'): return f'{cl}_{"intact" if ax == "intact_gene" else "fragment"}'
    return 'none_detected'

fams = sorted(set(U['Gar']) | set(U['Mac']))
cat_rows, cand_rows, fa, fa_c = [], [], [], []
for f in fams:
    G = rep_genome(f); u = U[G][f]; m = M[G].get(f, {}); seq = u['unit_seq'].upper(); b = big[(G, f)]
    other = 'Mac' if G == 'Gar' else 'Gar'
    notes = [f"59g decision {u['unit_content_rule']} (sequence_axis {u['sequence_axis']}, unit_organization {u['unit_organization']}, 59b ratio {u['ratio_59b'] or 'NA'})"]
    if u.get('gene_best_class') and u['gene_best_class'] not in ('5S', 'rRNA_other') and u['sequence_axis'] in ('intact_gene', 'gene_derived'):
        notes.append(f"gene context {u['gene_best_class']} {u['sequence_axis']}")
    if f in U[other]: notes.append(f"also in {other}: own unit {len(U[other][f]['unit_seq'])} bp, {U[other][f]['n_runs']} runs")
    if m.get('homolog_other_haplotype'): notes.append(f"homolog in other genome {m['homolog_other_haplotype']} (identity {m.get('homolog_identity')}, {m.get('homolog_status')})")
    klass, cur = u['unit_content_rule'], 'computational_draft'
    if f in HAND:   # docs/decisions/hand_decisions.tsv: the author's call overrides the rule label; both are kept
        h = HAND[f]; klass, cur = h['hand_decision'], 'hand_decision'
        notes.append(f"hand decision {h['decided_on']}: {h['hand_decision']} (rule label {h['rule_label_59g']}); loci: {h['locus_note']}")
    row = dict(family_id=f, preferred_name=f, aliases=json.dumps(sorted({x[0] for x in aliases[f]})), representative_taxon=TAXON[G], source_assembly_id=ASM[G],
               source_array_ids=json.dumps([b['run_id'], b['array_id_v46'] or b['array_id']]), consensus_length_bp=len(seq), consensus_sha256=sha(seq.encode()),
               repeat_class=klass, rrna_context=rrna(u), model_id='NA', model_status='not_attempted', curation_status=cur,
               library_version=a.version, provenance_uri=f"fClaHyb_v1.0_{G}_satelliteDNA.tar.gz: fClaHyb_{G}.satDNA_v47.unit_content_59g.tsv (unit_seq)", notes='; '.join(notes))
    if len(seq) >= 7 and m.get('low_copy_flag') == 'no': cat_rows.append(row); fa.append((f, f'{G} {TAXON[G]} {len(seq)} bp', seq))
    else: row['notes'] += '; candidate: ' + ('low copy (59h: < 10 copies or < 1 kb)' if m.get('low_copy_flag') != 'no' else 'no unit'); cand_rows.append(row); fa_c.append((f, f'{G} candidate', seq))

COLS = 'family_id preferred_name aliases representative_taxon source_assembly_id source_array_ids consensus_length_bp consensus_sha256 repeat_class rrna_context model_id model_status curation_status library_version provenance_uri notes'.split()
def wtsv(p, cols, rows):
    """plain TSV, no CSV quoting (JSON cells stay raw JSON); a tab or newline inside a cell is an error"""
    with open(p, 'w', encoding='utf-8', newline='') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            vals = [str(r[c]) for c in cols]
            if any(('\t' in v or '\n' in v or '\r' in v) for v in vals): sys.exit(f'tab/newline inside a cell: {p} {vals[0]}')
            fh.write('\t'.join(vals) + '\n')
def wfa(p, recs):
    with open(p, 'w', encoding='utf-8', newline='') as fh:
        for n, d, s in recs: fh.write(f'>{n} {d}\n' + '\n'.join(s[i:i + 60] for i in range(0, len(s), 60)) + '\n')
wfa(os.path.join(REL, 'SatChain_consensus.fa'), fa); wtsv(os.path.join(REL, 'SatChain_catalog.tsv'), COLS, cat_rows)
wtsv(os.path.join(REL, 'SatChain_candidates.tsv'), COLS, cand_rows); wfa(os.path.join(REL, 'optional', 'candidates_consensus.fa'), fa_c)
AL = []
for f in fams:   # one row per (alias, family); a label seen in both genomes keeps both sources
    grp = collections.defaultdict(list)
    for al, st, ra, src in sorted(set(aliases[f])): grp[(al, st)].append((ra, src))
    for (al, st), v in grp.items():
        AL.append(dict(alias=al, family_id=f, alias_status=st, first_release=a.version, last_release='', rationale=' | '.join(x[0] for x in v), source_record_id=' | '.join(x[1] for x in v)))
wtsv(os.path.join(REL, 'SatChain_aliases.tsv'), 'alias family_id alias_status first_release last_release rationale source_record_id'.split(), AL)
MOD = [dict(family_id=r['family_id'], model_id='NA', model_status='not_attempted', training_copy_count='NA', training_array_count='NA', training_assembly_ids='NA',
            alignment_sha256='NA', model_sha256='NA', model_length_bp='NA', calibration_status='NA', threshold_bits='NA', calibration_set_id='NA',
            notes='v1.0.0 is consensus-only by decision (2026-11-09); HMMs need aligned real copies and held-out calibration') for r in cat_rows + cand_rows]
wtsv(os.path.join(REL, 'SatChain_models.tsv'), 'family_id model_id model_status training_copy_count training_array_count training_assembly_ids alignment_sha256 model_sha256 model_length_bp calibration_status threshold_bits calibration_set_id notes'.split(), MOD)

# ---- checks (format section 'Input and final reproducibility checks', items 1, 3, 6)
ids = [r['family_id'] for r in cat_rows]
assert len(ids) == len(set(ids)) and [n for n, _, _ in fa] == ids, 'FASTA ids must equal catalogue ids'
assert all(sha(s.encode()) == r['consensus_sha256'] for (_, _, s), r in zip(fa, cat_rows))
known = set(ids) | {r['family_id'] for r in cand_rows}
assert all(x['family_id'] in known for x in AL), 'every alias resolves to a family'
ret = collections.defaultdict(set)
for x in AL:
    if x['alias_status'] == 'retired_do_not_use': ret[x['alias']].add(x['family_id'])
assert all(len(v) == 1 for v in ret.values()), f'retired alias maps to >1 family: {[k for k, v in ret.items() if len(v) > 1][:5]}'

counts = collections.Counter(r['repeat_class'] for r in cat_rows)
notes_md = f"""# SatChain satellite library v{a.version} (fClaHyb v1.0, SatChain v4.7)

Computational draft, not manually curated. Built {datetime.date.today()} by SatChain/satchain/export_library_v1.py from the v4.7 final outputs
(59g unit content 2026-11-07, 59h family metrics 2026-11-07/08, 59_catalogue, family nomenclature, per-run 5S axis 2026-11-09).

- `SatChain_consensus.fa`: {len(fa)} families (one consensus each): {', '.join(f'{k} {v}' for k, v in counts.most_common())}.
- `SatChain_candidates.tsv` + `optional/candidates_consensus.fa`: {len(cand_rows)} low-copy families (59h: < 10 copies or < 1 kb), same columns.
- `SatChain_aliases.tsv`: {len(AL)} aliases; retired labels (CgaSat/CmaSat names, v4.6 gene-array labels) are `retired_do_not_use`; step-59 pre-naming ids are `internal_provisional`.
- `SatChain_models.tsv`: every family `not_attempted` (no HMMs in this version).

Deviations from the proposed format (LIBRARY_RELEASE_FORMAT_v1): family ids are the deposited GarSatNN / MacSatNN / NewSat<H>nnnnn ids, not SAT accessions;
no `SatChain_profiles.hmm.gz`; consensus sequences are not phase- or orientation-normalised (step-03 consensus as written).

Classes: `repeat_class` = the 59g decision (satellite / gene_array / gene_derived_satellite / unresolved). `rrna_context` = 5S or other rRNA
homology of the consensus (intact = identity >= 0.90 and reference coverage >= 0.80; fragment = >= 0.70 and >= 0.50), independent of the class.
A family found in both genomes has one row; the consensus comes from the genome of its id series and the notes give the other genome's unit.
Hand decisions (docs/decisions/hand_decisions.tsv) override the rule label and set curation_status = hand_decision: {', '.join(f"{k} -> {v['hand_decision']}" for k, v in HAND.items()) or 'none'}.
Loci can differ inside a family (MacSat06: intact 5S gene arrays on LG01, LG02, LG08, LG20, LG25; 5S-derived units only on LG03 and LG22).

Not validated: no held-out benchmark, no transfer test to other genomes, no thresholds. Consensus sequences are references, not evidence of
shared origin or function.
"""
open(os.path.join(REL, 'SatChain_release_notes.md'), 'w', encoding='utf-8', newline='\n').write(notes_md)
asm = {G: dict(accession=ASM[G], fasta=f'fClaHyb_{G}_LG.fa', sha256=fsha(os.path.join(a.ann, f'fClaHyb_{G}', f'fClaHyb_{G}_genome', f'fClaHyb_{G}_LG.fa'))) for G in ('Gar', 'Mac')}
files = sorted(os.path.relpath(os.path.join(r, f), REL).replace('\\', '/') for r, _, fs in os.walk(REL) for f in fs)
man = dict(library_name='SatChain satellite library (fClaHyb v1.0)', library_version=a.version, schema_version='LIBRARY_RELEASE_FORMAT_v1 (with declared deviations)',
           created_at=datetime.datetime.now().isoformat(timespec='seconds'), source_worklane_runs=['SatChain v4.7 final run 2026-11-02', '59g 2026-11-07', '59h 2026-11-07/08', '59g_run_axis 2026-11-09'],
           source_assembly_hashes=asm, total_families=len(cat_rows), total_candidates=len(cand_rows), total_profiles=0,
           files=[dict(path=p, sha256=fsha(os.path.join(REL, p)), bytes=os.path.getsize(os.path.join(REL, p))) for p in files],
           normalization_method='none (step-03 consensus as written)', builder_version='export_library_v1.py 2026-11-09', hmm_tool_versions=None,
           calibration_policy='no thresholds; no held-out evaluation', known_limitations=['computational draft', 'no HMMs', 'no phase/orientation normalisation', 'family-level class provisional where loci differ'],
           curation_approval='none (computational draft)')
json.dump(man, open(os.path.join(REL, 'SatChain_manifest.json'), 'w', encoding='utf-8'), indent=1)
files = sorted(os.path.relpath(os.path.join(r, f), REL).replace('\\', '/') for r, _, fs in os.walk(REL) for f in fs)
open(os.path.join(REL, 'SHA256SUMS'), 'w', newline='\n').write(''.join(f'{fsha(os.path.join(REL, p))}  {p}\n' for p in files))
chk = {l.split('  ')[1].strip(): l.split('  ')[0] for l in open(os.path.join(REL, 'SHA256SUMS'))}
assert all(fsha(os.path.join(REL, p)) == h for p, h in chk.items())
print(f'library v{a.version}: {len(cat_rows)} families, {len(cand_rows)} candidates, {len(AL)} aliases; classes {dict(counts)}; checks passed -> {REL}')
