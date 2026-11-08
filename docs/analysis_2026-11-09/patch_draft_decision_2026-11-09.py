"""Apply the 2026-11-09 DRAFT-step decision to run_satchain.sh, STEP_METHOD.md, README.md and the 2026-11-07 script register.
Each replacement must match exactly once; otherwise the script stops before writing anything."""
import os, sys, csv, io
R = r'F:\data\SatChain'
def sub1(text, old, new, label):
    if text.count(old) != 1: sys.exit(f'{label}: expected 1 match, found {text.count(old)}')
    return text.replace(old, new)

W = os.path.join(R, 'workflow', 'run_satchain.sh'); w = open(W, encoding='utf-8', newline='').read()
old61 = w[w.index('# 61 / 61b: unit boundaries'):w.index('# 62-65 (exploratory')]
new61 = ('# 61 / 61b: EXPLORATORY since 2026-11-09 (docs/analysis_2026-11-09/DRAFT_STEPS_DECISION_2026-11-09.md): no result needs base-pair\n'
         '# edges (composite units ~13.6 kb; the 56 edge is 9-154 bp inside the satellite). Not part of the v4.7 chain; code and synthetic tests kept.\n'
         + ''.join('# ' + l + '\n' if l.strip() else l + '\n' for l in old61.rstrip('\n').split('\n')[1:]))
w = w.replace(old61, new61)
w = sub1(w, '# 62-65 (exploratory until their sweeps are on disk): see docs/STEP_METHOD.md section 8 and each script\'s usage line',
         '# 62 / 65: run NARROW only (2026-11-09) on the HOR-like centromeric arrays of thesis Table 3.9, to re-test the MacSat1272 HOR claim;\n'
         '#   driver docs/analysis_2026-11-09/run_hor_narrow.py. 53c, 60a, 63, 64: EXPLORATORY (same decision note). See docs/STEP_METHOD.md section 8.', 'sh 62-65')
open(W, 'w', encoding='utf-8', newline='').write(w)

S = os.path.join(R, 'docs', 'STEP_METHOD.md'); s = open(S, encoding='utf-8', newline='').read()
s = sub1(s, '## 9. Topology is descriptive today',
         'Decision 2026-11-09 (docs/analysis_2026-11-09/DRAFT_STEPS_DECISION_2026-11-09.md): 53c, 60a, 61, 61b, 63 and 64 are EXPLORATORY, not\n'
         'part of the v4.7 results (no result, figure or thesis claim depends on them; 61/61b removed from run_satchain.sh). 62 and 65 are run\n'
         'narrow on the HOR-like centromeric arrays only, with their sweeps, to re-test the one HOR claim of the thesis (MacSat1272, Mac LG20).\n\n'
         '## 9. Topology is descriptive today', 'STEP_METHOD')
open(S, 'w', encoding='utf-8', newline='').write(s)

M = os.path.join(R, 'README.md'); m = open(M, encoding='utf-8', newline='').read()
m = sub1(m, 'family-specific k-mers for transcription tests in RNA-seq (in development) |',
         'family-specific k-mers for transcription tests in RNA-seq; **exploratory** (2026-11-09) |', 'README 60a')
m = sub1(m, 'expansion tracks, object profiles; exploratory until their `--sweep` tables are on disk (`docs/STEP_METHOD.md` section 8) |',
         'expansion tracks, object profiles; **2026-11-09:** 61, 61b, 63, 64 exploratory (not in the v4.7 chain), 62 / 65 run narrow on the HOR-like centromeric arrays only (`docs/analysis_2026-11-09/`) |', 'README 61-65')
m = sub1(m, 'windows, change points, topology (memory, structure index, core/accessory), topology library |',
         'windows, change points, topology (memory, structure index, core/accessory), topology library (53c: **exploratory**, 2026-11-09) |', 'README 53c')
open(M, 'w', encoding='utf-8', newline='').write(m)

G = os.path.join(R, 'docs', 'analysis_2026-11-07', 'SCRIPT_REGISTER_2026-11-07.tsv')
rows = list(csv.reader(open(G, encoding='utf-8', newline=''), delimiter='\t')); h = rows[0]
iS, iN, iSt = h.index('script'), h.index('fate_or_note'), h.index('status')
NEW = {'53c_library_pool.py': 'EXPLORATORY', '60a_tx_kmers.py': 'EXPLORATORY', '61_unit_boundaries.py': 'EXPLORATORY',
       '61b_jsd_boundaries.py': 'EXPLORATORY', '63_expansion_tracks.py': 'EXPLORATORY', '64_object_profiles.py': 'EXPLORATORY',
       '62_hor_detect.py': 'NARROW_RUN', '65_hor_mining.py': 'NARROW_RUN'}
hit = 0
for r in rows[1:]:
    if r[0] == 'SatChain_repo' and r[iS] in NEW:
        r[iSt] = NEW[r[iS]]; r[iN] += f' | 2026-11-09: {NEW[r[iS]]} (DRAFT_STEPS_DECISION_2026-11-09.md)'; hit += 1
if hit != len(NEW): sys.exit(f'register: {hit} of {len(NEW)} rows found')
out = io.StringIO(); csv.writer(out, delimiter='\t', lineterminator='\r\n').writerows(rows)   # the register is CRLF
open(G, 'w', encoding='utf-8', newline='').write(out.getvalue())
print('patched: run_satchain.sh, STEP_METHOD.md, README.md, register (8 rows)')
