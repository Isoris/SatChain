#!/bin/bash
# Step 1-2 (2026-10-22): ORFs of the MacSat03-interleaved element (EDTA Mac TE_00003409_INT, 8,164 bp) and their nearest protein domains.
# No HMMER (user): the reference set is the consensus sequence stored in every profile of TEsorter's REXdb metazoa v3, GyDB 2 and Kapitonov LINE
# databases (CONS column of each match state), searched with BLASTP. Output: domains/orfs.faa, domains/orf_vs_refdomains.tsv (hits per ORF segment).
set -eo pipefail   # no -u: conda activate scripts use unset variables
source ~/mambaforge/etc/profile.d/conda.sh; conda activate assembly
W=/mnt/f/data/satellite_FISH/LTR_MacSat03_2026-10-22
DB=$(ls -d ~/mambaforge/envs/EDTA2.2/lib/python3.10/site-packages/TEsorter/database)
python3 - "$DB" "$W/domains/refdomains.faa" <<'PY'
import sys, os
db, out = sys.argv[1], open(sys.argv[2], 'w'); n = 0
for f in ('REXdb_protein_database_metazoa_v3.hmm', 'GyDB2.hmm', 'Kapitonov_et_al.GENE.LINE.hmm'):
    name, seq, inhmm = None, [], False
    for l in open(os.path.join(db, f)):
        if l.startswith('NAME'): name = l.split(None, 1)[1].strip(); seq = []; inhmm = False
        elif l.startswith('HMM '): inhmm = True
        elif l.startswith('//'):
            if name and seq: out.write('>' + f.split('.')[0].split('_')[0] + '|' + name + '\n' + ''.join(seq) + '\n'); n += 1
            name = None
        elif inhmm:
            p = l.split()
            if len(p) >= 26 and p[0].isdigit():          # match-state line: idx, 20 scores, MAP, CONS, RF, MM, CS
                c = p[22]
                seq.append(c.upper() if c not in '-.' else 'X')
print('reference domain consensus sequences:', n, file=sys.stderr)
PY
getorf -sequence $W/seq/TE_00003409_library.fa -outseq $W/domains/orfs.faa -minsize 450 -find 1 >/dev/null 2>&1
echo "ORFs >= 150 aa: $(grep -c '>' $W/domains/orfs.faa)"
makeblastdb -in $W/domains/refdomains.faa -dbtype prot -out $W/domains/refdomains >/dev/null
blastp -query $W/domains/orfs.faa -db $W/domains/refdomains -evalue 1e-5 -max_target_seqs 50 -num_threads 2 \
  -outfmt '6 qseqid qstart qend sseqid sstart send pident length evalue bitscore slen' > $W/domains/orf_vs_refdomains.tsv
echo "hits: $(wc -l < $W/domains/orf_vs_refdomains.tsv)"
