#!/usr/bin/env python3
"""03b_dedup_copies.py  (2026-11-01, audit CB01): apply the tile clean-up of 03_monomer_copies.ordered_valid_calls to an EXISTING
03_monomer_copies.tsv, without re-running StringDecomposer.

StringDecomposer (default batch 5,000 bp) reports the tiles of every batch-overlap zone twice (identical coordinates, consecutive
monomer_index) or as partially overlapping tiles: 4.3 % (Gar) / 5.2 % (Mac) of the v4.5 rows, 99.7 % of them within 400 bp after a
5-kb multiple of array_relative_start (repro: audit coreA chk_dup.py). The rule is imported from 03_monomer_copies.py, so both
paths agree: exact duplicates -> the higher-identity row is kept; overlap > OVL_FRAC (0.5) of the shorter tile -> the higher-identity
tile is kept; smaller overlap -> the later tile's start is trimmed to the previous end.

Per array (arrays are contiguous in the 03 tables; checked, two-pass fallback by array_id otherwise), every column is kept:
  - monomer_index and monomer_id are renumbered, boundary_copy recomputed (first / last kept copy);
  - a trimmed row gets genomic_start_0based / array_relative_start / observed_length_bp / sequence / sequence_sha1 / AT_percent /
    GC_percent updated (the sequence is oriented: a '-' strand copy loses the bases at ITS end) and WARN_overlap_trimmed added to
    qc_flag; best_alignment_identity, matches / mismatches / insertions / deletions, divergence and complete_or_truncated stay the
    StringDecomposer values of the untrimmed tile (they cannot be recomputed without the class consensus).
Outputs: <out>.tsv (same columns) and <out>.dedup_log.tsv (array_id, rows_in, rows_out, dup_dropped, overlap_dropped,
overlap_trimmed). Memory: one array at a time (the largest v4.5 array is ~100k rows). Needs edlib importable (03 imports it).
usage: 03b_dedup_copies.py <03_monomer_copies.tsv> <out prefix>"""
import sys, os, csv, hashlib, importlib.util, collections
csv.field_size_limit(1 << 30)
HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('monomer_copies_03', os.path.join(HERE, '03_monomer_copies.py'))
_m03 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_m03)
ordered_valid_calls = _m03.ordered_valid_calls


def dedup_array(aid, rows, col):
    """rows: list of field lists of one array (same order as the file); returns (kept rows, stats). Coordinates are taken from
    array_relative_start / array_relative_end (half-open) and handed to ordered_valid_calls as inclusive SD tiles."""
    tiles = [[str(i), r[col['monomer_class']], r[col['array_relative_start']], str(int(r[col['array_relative_end']]) - 1),
              r[col['best_alignment_identity']]] for i, r in enumerate(rows)]
    arr_len = max(int(r[col['array_relative_end']]) for r in rows) + 1
    st = {}
    kept = ordered_valid_calls(tiles, arr_len, st)
    out = []
    for mi, t in enumerate(kept):
        r = list(rows[int(t[0])]); s0 = int(t[2]); old_s0 = int(r[col['array_relative_start']])
        if s0 != old_s0:                                                   # trimmed: drop the first (s0 - old_s0) genomic bases
            cut = s0 - old_s0; seq = r[col['sequence']]
            seq = seq[cut:] if r[col['strand']] == '+' else seq[:len(seq) - cut]   # oriented sequence of a '-' copy ends where the genome starts
            r[col['sequence']] = seq; r[col['sequence_sha1']] = hashlib.sha1(seq.encode()).hexdigest()[:16]
            r[col['array_relative_start']] = str(s0)
            r[col['genomic_start_0based']] = str(int(r[col['genomic_start_0based']]) + cut)
            r[col['observed_length_bp']] = str(len(seq))
            nN = seq.count('N'); atp = 100 * (seq.count('A') + seq.count('T')) / max(1, len(seq) - nN)
            r[col['AT_percent']] = f'{atp:.2f}'; r[col['GC_percent']] = f'{100 - atp:.2f}'
            q = r[col['qc_flag']]; flags = [] if q == 'PASS' else q.split(';')
            if 'WARN_overlap_trimmed' not in flags: flags.append('WARN_overlap_trimmed')
            r[col['qc_flag']] = ';'.join(flags)
        r[col['monomer_index']] = str(mi + 1); r[col['monomer_id']] = f'{aid}_M{mi + 1:05d}'
        r[col['boundary_copy']] = 'yes' if (mi == 0 or mi == len(kept) - 1) else 'no'
        out.append(r)
    return out, st


def main():
    src, pre = sys.argv[1:3]
    with open(src, newline='') as fh: hdr = fh.readline().rstrip('\n').split('\t')
    col = {c: i for i, c in enumerate(hdr)}; ai = col['array_id']
    tot = collections.Counter(); narr = 0
    with open(pre + '.tsv', 'w', newline='') as out, open(pre + '.dedup_log.tsv', 'w', newline='') as lg:
        w = csv.writer(out, delimiter='\t', lineterminator='\n'); w.writerow(hdr)
        lw = csv.writer(lg, delimiter='\t', lineterminator='\n')
        lw.writerow(['array_id', 'rows_in', 'rows_out', 'dup_dropped', 'overlap_dropped', 'overlap_trimmed'])

        def flush(aid, rows):
            nonlocal narr
            kept, st = dedup_array(aid, rows, col); w.writerows(kept); narr += 1
            lw.writerow([aid, len(rows), len(kept), st['dup_dropped'], st['overlap_dropped'], st['overlap_trimmed']])
            tot['in'] += len(rows); tot['out'] += len(kept)
            for k in ('dup_dropped', 'overlap_dropped', 'overlap_trimmed'): tot[k] += st[k]

        seen = set(); cur = None; rows = []; contiguous = True
        with open(src, newline='') as fh:
            fh.readline()
            for line in fh:
                r = line.rstrip('\n').split('\t'); aid = r[ai]
                if aid != cur:
                    if cur is not None: flush(cur, rows); seen.add(cur)
                    if aid in seen: contiguous = False; break
                    cur = aid; rows = []
                rows.append(r)
            if contiguous and cur is not None: flush(cur, rows)
        if not contiguous:
            # two-pass fallback: arrays are not contiguous -> pass 1 records the line offsets per array, pass 2 reads each array
            print(f'    arrays not contiguous in {src} (array {aid} re-appears): two-pass mode by array_id', flush=True)
            out.seek(0); out.truncate(); w.writerow(hdr); lg.seek(0); lg.truncate()
            lw.writerow(['array_id', 'rows_in', 'rows_out', 'dup_dropped', 'overlap_dropped', 'overlap_trimmed'])
            tot.clear(); narr = 0
            offs = collections.OrderedDict()
            with open(src, 'rb') as fh:
                fh.readline(); pos = fh.tell()
                for line in fh:
                    offs.setdefault(line.split(b'\t', ai + 1)[ai].decode(), []).append(pos); pos += len(line)
            with open(src, 'rb') as fh:
                for aid, ps in offs.items():
                    rows = []
                    for p in ps: fh.seek(p); rows.append(fh.readline().decode().rstrip('\n').split('\t'))
                    flush(aid, rows)
    print(f'    {narr} arrays; rows in {tot["in"]:,} -> out {tot["out"]:,}; dropped: {tot["dup_dropped"]:,} exact duplicates, '
          f'{tot["overlap_dropped"]:,} large overlaps; trimmed {tot["overlap_trimmed"]:,}; wrote {pre}.tsv and {pre}.dedup_log.tsv')


if __name__ == '__main__':
    main()
