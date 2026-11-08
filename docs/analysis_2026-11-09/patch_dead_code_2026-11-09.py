"""Remove dead code found by the 2026-11-07 function trace (CLASSIFICATION_DECISION ledger row 'dead code'):
03 read_fasta + LIB_EPS, 59 --gaps option, runs_of_record merged_bp. 02 --threads and 11 mismatch variants are in the frozen
v4.5 pipeline (satellite_FISH/pipeline_v4_5_fixed, deposited code) and are left unchanged. Run once; keeps each file's newlines."""
import sys
R = r'F:\data\SatChain\satchain'
EDITS = {
    '03_monomer_copies.py': [
        ("LIB_EPS       = 0.98    # library consensus similarity, reported\n", ""),
        ("# COMPARISON RECOMMENDATION [DESIGN_ONLY]: Keep only if required; compare against indexed sequence extraction with source checksums.\n"
         "def read_fasta(path):\n"
         "    n, s = None, []\n"
         "    for line in open(path):\n"
         "        line = line.rstrip()\n"
         "        if line.startswith('>'):\n"
         "            if n: yield n, ''.join(s).upper()\n"
         "            n = line[1:].split()[0]; s = []\n"
         "        else: s.append(line)\n"
         "    if n: yield n, ''.join(s).upper()\n\n\n", ""),
    ],
    '59_final_catalogue.py': [
        (" [--validated-name 56_runs_validated_q0.9.bed] [--gaps <bed>]", " [--validated-name 56_runs_validated_q0.9.bed]"),
        ("; ap.add_argument('--gaps', default='')", ""),
    ],
    'runs_of_record.py': [
        ("def merged_bp(iv):\n"
         "    iv = sorted(iv); t = 0; cur = None\n"
         "    for s, e in iv:\n"
         "        if cur and s <= cur[1]: cur[1] = max(cur[1], e)\n"
         "        else:\n"
         "            if cur: t += cur[1] - cur[0]\n"
         "            cur = [s, e]\n"
         "    return t + (cur[1] - cur[0] if cur else 0)\n\n\n", ""),
    ],
}
for f, eds in EDITS.items():
    p = f'{R}\\{f}'; raw = open(p, 'rb').read(); crlf = b'\r\n' in raw
    s = raw.decode('utf-8').replace('\r\n', '\n')
    for old, new in eds:
        if s.count(old) != 1: sys.exit(f'{f}: expected 1 match, found {s.count(old)}: {old[:60]!r}')
        s = s.replace(old, new)
    if crlf: s = s.replace('\n', '\r\n')
    open(p, 'wb').write(s.encode('utf-8')); print(f, 'patched', 'CRLF' if crlf else 'LF')
