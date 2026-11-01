"""gaps.py  (SatChain v4.7): assembly gaps as hard stops.
A gap is an interval of the assembly-gap BED (detgaps; N placeholders, e.g. 100 or 200 N) or a zero-length row (two contigs joined without
Ns: the sequence across the join is not observed either). Optionally N runs found directly in the FASTA are added (--scan-n), so a catalogue built
on a genome whose BED is missing still stops at Ns.
between(c, a, b): is there a gap (or join point) strictly inside (a, b) or overlapping [a, b)?  Used by 56 (a run ends at a gap), 59 (arrays
split into gap-free segments), 59e and 61 (no interruption / chain link across a gap), 61b (JSD windows never touch a gap).
"""
import bisect, collections, re


class Gaps:
    def __init__(self, bed=None, fasta=None, min_n=10, include_joins=True):
        self.iv = collections.defaultdict(list)
        if bed:
            for l in open(bed):
                if not l.strip() or l.startswith(('#', 'track')): continue
                f = l.split('\t'); c, a, b = f[0], int(f[1]), int(f[2])
                if a == b and (not include_joins or a == 0): continue          # a 0-0 row is the chromosome start, not a join
                self.iv[c].append((a, b))
        if fasta: self.scan_fasta(fasta, min_n)
        for c in self.iv: self.iv[c] = sorted(set(self.iv[c]))
        self.starts = {c: [x[0] for x in v] for c, v in self.iv.items()}

    def scan_fasta(self, fasta, min_n):
        name, buf, pos = None, [], 0
        def flush():
            s = ''.join(buf).upper()
            for m in re.finditer('N{%d,}' % min_n, s): self.iv[name].append((m.start(), m.end()))
        for l in open(fasta):
            if l.startswith('>'):
                if name: flush()
                name, buf = l[1:].split()[0], []
            else: buf.append(l.strip())
        if name: flush()

    def between(self, c, a, b, closed=False):
        """True when a gap overlaps [a, b) or a join point lies strictly inside (a, b).
        closed=True (2026-11-01, audit HB07): a zero-length join at a or b also counts. Use it for the stretch BETWEEN two
        neighbouring copies / runs (last_end, next_start): a join exactly on a tile boundary separates the tiles, and the open test
        let a chain continue across it. Keep closed=False for the tile itself (a join at its own edge does not make it a gap copy)."""
        v = self.iv.get(c)
        if not v or b < a or (b == a and not closed): return False
        i = bisect.bisect_left(self.starts[c], b + (1 if closed else 0)) - 1
        while i >= 0:
            s, e = v[i]
            if s == e:
                if (a <= s <= b) if closed else (a < s < b): return True
            elif e > a and s < b: return True
            if e <= a and s < a - 10_000_000: break
            i -= 1
        return False

    def gaps_in(self, c, a, b):
        """all gap / join rows of chromosome c that between(c, a, b, closed=True) would see (for reports)"""
        return [(s, e) for s, e in self.iv.get(c, []) if (s == e and a <= s <= b) or (s != e and e > a and s < b)]

    def count(self): return sum(len(v) for v in self.iv.values())