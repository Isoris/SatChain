#!/usr/bin/env python3
"""
faidx_lite.py -- random access to a FASTA without loading it.

Four scripts previously did dict(read_fasta(genome)), which holds the entire
assembly in RAM. On a ~1 Gb genome that is ~1 GB of Python string plus overhead,
which does not fit a 4-5 GB laptop alongside numpy working arrays.
"""
import os, re


class Fai:
    # COMPARISON RECOMMENDATION [DESIGN_ONLY]: Verify .fai is current for this exact assembly checksum.
    def __init__(self, path):
        self.path = path
        self.idx = {}
        for name, off, ln in self._index(path):
            self.idx[name] = (off, ln)

    @staticmethod
    # COMPARISON RECOMMENDATION [DESIGN_ONLY]: Test FASTA index geometry against samtools faidx; avoid assumptions about wrapped lines.
    def _index(path):
        recs = []
        with open(path, 'rb') as fh:
            name, start, pos = None, None, 0
            for line in fh:
                if line.startswith(b'>'):
                    if name is not None:
                        recs.append((name, start, pos - start))
                    name = line[1:].split()[0].decode()
                    start = pos + len(line)
                pos += len(line)
            if name is not None:
                recs.append((name, start, pos - start))
        return recs

    # COMPARISON RECOMMENDATION [DESIGN_ONLY]: Maintain canonical contig IDs without collapsing Gar and Mac namespaces.
    def names(self):
        return list(self.idx)

    # COMPARISON RECOMMENDATION [DESIGN_ONLY]: Validate 0-based/1-based and half-open extraction with edge-case tests.
    def get(self, name):
        off, ln = self.idx[name]
        with open(self.path, 'rb') as fh:
            fh.seek(off)
            return fh.read(ln).replace(b'\n', b'').replace(b'\r', b'').decode().upper()

    # COMPARISON RECOMMENDATION [DESIGN_ONLY]: Use source sequence bounds as invariant for array and monomer calls.
    def length(self, name):
        return len(self.get(name))

    # COMPARISON RECOMMENDATION [DESIGN_ONLY]: Keep contig-length map for window denominators and IGV export.
    def sizes(self):
        """lengths without holding any sequence"""
        out = {}
        for n, (off, ln) in self.idx.items():
            with open(self.path, 'rb') as fh:
                fh.seek(off)
                blk = fh.read(ln)
            out[n] = len(blk) - blk.count(b'\n') - blk.count(b'\r')
            del blk
        return out
