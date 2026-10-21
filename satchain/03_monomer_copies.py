#!/usr/bin/env python3
"""
03_monomer_copies.py  --  iterative multi-consensus decomposition.

The naive form of this step, one consensus per array, fails whenever an array is
heterogeneous. If an array contains

    A A A B A A C B A A          or        A A A A A A | B B B B B B

then a single consensus X lands between A and B and is neither. Every copy then
scores badly against it. That is exactly what produced the 66-69% arrays in the
first pass, while homogeneous arrays reached 90-96%.

So the loop is:

    approximate period (prior, NOT a fixed copy length)
      -> provisional indel-aware boundaries (StringDecomposer, one seed)
      -> extract copies, normalise orientation and circular phase
      -> cluster copies into monomer CLASSES
      -> consensus per class  => a monomer LIBRARY, not a single monomer
      -> re-decompose the whole array against the library
      -> repeat until boundaries and class assignments stop moving

Clustering uses identity AND coverage AND length AND completeness, not identity
alone, so that:
  - small copy-specific indels stay inside a class
  - a recurrent indel shared by many copies can define a class
  - terminal fragments are labelled truncated, never promoted to a class

Convergence currently uses repeated boundary stability. Library similarity is
reported separately but does not independently establish convergence. Class assignment stability is NOT implemented as
an independent criterion. `unresolved_fraction` is the share of array positions
not covered by a decomposer interval; it is not a confidence-filtered fraction.
"""
import argparse, csv, os, sys, subprocess, tempfile, collections, shutil, hashlib, gc
import numpy as np
import edlib

# --- clustering
CLASS_ID      = 90.0    # % identity to join a monomer class
CLASS_COV     = 0.85    # min alignment coverage of the shorter sequence
LEN_TOL       = 0.15    # max length difference within a class
MIN_CLASS     = 8       # copies needed before a class earns a library entry
MAX_CLASSES   = 12      # per array; beyond this the array is flagged, not split further
# --- orientation
ORIENT_MARGIN = 3.0
PALIN_ID      = 80.0
# --- iteration
MAX_ITER      = 5
BOUNDARY_EPS  = 2.0     # mean |delta start| in bp
LIB_EPS       = 0.98    # library consensus similarity, reported
STABLE_ROUNDS = 2       # consecutive rounds with boundary shift < BOUNDARY_EPS
MAX_BOOT_BP   = 150000
CLUSTER_SAMPLE = 600    # copies sampled for class discovery

COMP = str.maketrans('ACGTN', 'TGCAN')
# COMPARISON RECOMMENDATION [DESIGN_ONLY]: Keep reverse-complement primitive; validate orientation on palindromic and IUPAC-containing copies.
def rc(s): return s.translate(COMP)[::-1]
# COMPARISON RECOMMENDATION [DESIGN_ONLY]: Keep sequence GC as descriptive; separately report masked/ambiguous base fraction.
def gc_(s): return (s.count('G') + s.count('C')) / max(1, len(s))


# COMPARISON RECOMMENDATION [DESIGN_ONLY]: Preserve ordered_valid_calls behavior while comparing same-grain outputs, provenance and boundary conditions.
def ordered_valid_calls(calls, array_length):
    """Preserve observed genomic order; fail rather than silently misplace copies."""
    ordered = sorted(calls, key=lambda f: (int(f[2]), int(f[3])))
    for f in ordered:
        start, end = int(f[2]), int(f[3])  # StringDecomposer: inclusive end
        if start < 0 or end < start or end >= array_length:
            raise ValueError(f'invalid StringDecomposer interval [{start},{end}] in length={array_length}')
    return ordered


# COMPARISON RECOMMENDATION [DESIGN_ONLY]: Preserve nearest_boundary_shift behavior while comparing same-grain outputs, provenance and boundary conditions.
def nearest_boundary_shift(previous, current):
    """None means shift cannot be assessed (empty/insufficient calls)."""
    if not previous or not current:
        return None
    ps = np.sort(np.asarray(previous, dtype=np.int64))
    cs = np.asarray(current, dtype=np.int64)
    j = np.searchsorted(ps, cs)
    left = ps[np.maximum(j - 1, 0)]
    right = ps[np.minimum(j, len(ps) - 1)]
    return float(np.minimum(np.abs(cs - left), np.abs(cs - right)).mean())


# COMPARISON RECOMMENDATION [DESIGN_ONLY]: Keep only if required; compare against indexed sequence extraction with source checksums.
def read_fasta(path):
    n, s = None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith('>'):
            if n: yield n, ''.join(s).upper()
            n = line[1:].split()[0]; s = []
        else: s.append(line)
    if n: yield n, ''.join(s).upper()


# COMPARISON RECOMMENDATION [DESIGN_ONLY]: Keep StringDecomposer placement as primary; compare its boundaries with independent builtin/TRASH/ArraySplitter results at the same array_id.
def run_sd(seq_fa, mono_fa, outdir, second_best=False):
    os.makedirs(outdir, exist_ok=True)
    cmd = ['stringdecomposer', seq_fa, mono_fa, '-o', outdir, '-t', '1']
    if second_best: cmd.append('--second-best')
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    return [l.rstrip('\n').split('\t')
            for l in open(os.path.join(outdir, 'final_decomposition.tsv'))]


# COMPARISON RECOMMENDATION [DESIGN_ONLY]: Keep rotation/strand-aware recognition; benchmark edlib normalization, unequal length, true cyclic offset/strand estimates.
def circ_match(a, b):
    """rotation- and strand-invariant identity + coverage of a against b"""
    best = (0.0, 0.0, '+')
    for st, bb in (('+', b), ('-', rc(b))):
        r = edlib.align(a, bb + bb, mode='HW')
        if r['editDistance'] < 0: continue
        pid = 100.0 * (1 - r['editDistance'] / len(a))
        cov = min(1.0, len(a) / max(1, len(bb)))
        if pid > best[0]: best = (pid, cov, st)
    return best


# COMPARISON RECOMMENDATION [DESIGN_ONLY]: Keep identity+coverage+length gates; fan out threshold sensitivity and alternative linkage models.
def same_class(a, b):
    if abs(len(a) - len(b)) > LEN_TOL * max(len(a), len(b)):
        return False
    pid, cov, _ = circ_match(a, b)
    return pid >= CLASS_ID and cov >= CLASS_COV


# COMPARISON RECOMMENDATION [DESIGN_ONLY]: Keep local per-array classes; test overflow MAX_CLASSES=12 without forced merging and compare greedy vs hierarchical clustering.
def cluster_copies(copies):
    """greedy on frequency-ordered seeds; returns list of member-index lists"""
    order = sorted(range(len(copies)), key=lambda i: -len(copies[i]))
    reps, lab = [], [-1] * len(copies)
    for i in order:
        hit = -1
        for k, ri in enumerate(reps):
            if same_class(copies[i], copies[ri]):
                hit = k; break
        if hit < 0:
            if len(reps) >= MAX_CLASSES:
                hit = int(np.argmin([edlib.align(copies[i], copies[r], mode='NW')['editDistance']
                                     for r in reps]))
            else:
                reps.append(i); hit = len(reps) - 1
        lab[i] = hit
    groups = collections.defaultdict(list)
    for i, l in enumerate(lab): groups[l].append(i)
    return [g for _, g in sorted(groups.items(), key=lambda kv: -len(kv[1]))]


# COMPARISON RECOMMENDATION [DESIGN_ONLY]: REVIEW: first explicitly phase-register rotations/strand per member, then align indels; test against varied-start synthetic copies.
def consensus_of(members, ref):
    L = len(ref)
    cnt = [collections.Counter() for _ in range(L)]
    for m in members:
        r = edlib.align(m, ref, task='path', mode='NW')
        ni = edlib.getNiceAlignment(r, m, ref)
        ci = 0
        for x, y in zip(ni['query_aligned'], ni['target_aligned']):
            if y == '-': continue
            if ci < L and x != 'N': cnt[ci][x] += 1      # N never votes
            ci += 1
    return ''.join(cnt[i].most_common(1)[0][0]
                   if cnt[i] and cnt[i].most_common(1)[0][0] != '-' else ref[i]
                   for i in range(L))


COLUMNS = ['assembly_id','species','haplotype','array_id','monomer_id','chrom',
           'genomic_start_0based','genomic_end','array_relative_start','array_relative_end',
           'monomer_index','expected_period_bp','observed_length_bp',
           'strand','strand_status','reverse_complemented','phase_offset_bp',
           'family_id','subfamily_id','monomer_class','consensus_version',
           'best_alignment_identity','second_best_family_identity','orientation_margin',
           'class_margin','divergence','matches','mismatches','insertions','deletions',
           'complete_or_truncated','boundary_copy','AT_percent','GC_percent',
           'sequence_sha1','source_tool','tool_version','run_id','qc_flag','sequence']


# COMPARISON RECOMMENDATION [DESIGN_ONLY]: Preserve iterative per-array library bootstrap; compare fixed seed vs multiple seeds and validate convergence, skipped-array inventory, complete physical paths.
# FIX v4.5.1: the per-array bootstrap is a pure function of (array row, genome),
# so arrays are distributed over a process pool (--threads workers, SD -t 1
# each). The parent writes rows in array order and checkpoints per array, as
# before. A failure in one array is recorded as FAILED:<reason> in
# 03_array_status.tsv instead of aborting the whole step.
_W = {'fai': None, 'chrom': None, 'seq': None, 'tag': None, 'run_id': None}


def _worker_init(genome, tag, run_id):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from faidx_lite import Fai
    _W['fai'] = Fai(genome); _W['tag'] = tag; _W['run_id'] = run_id


# COMPARISON RECOMMENDATION [DESIGN_ONLY]: Cache by chromosome; verify faidx index, coordinate convention and checksum on per-array extraction.
def _chromseq(c):
    if _W['chrom'] != c:
        _W['seq'] = None; _W['chrom'] = c; _W['seq'] = _W['fai'].get(c)
    return _W['seq']


def decompose_array(r):
    """one array -> dict(aid, status, iters, cons, mono, nlib); never raises"""
    aid = r['array_id']; tag = _W['tag']; run_id = _W['run_id']
    out = {'aid': aid, 'status': 'NOT_PROCESSED', 'iters': [], 'cons': [], 'mono': [], 'nlib': 0}
    try:
        chrom = r['chrom']; st = int(r['start_0based']); en = int(r['end'])
        per = int(r['period_mode_bp'])
        if per <= 0 or en <= st:
            out['status'] = 'INVALID_PERIOD_OR_INTERVAL'; return out
        full = _chromseq(chrom)[st:en]
        boot = full[:MAX_BOOT_BP]
        if len(boot) < 3 * per:
            out['status'] = 'SKIPPED_FEWER_THAN_THREE_PERIODS'; return out
        tmp = tempfile.mkdtemp(prefix='sd_')
        try:
            # SATELLITE v4.2 SOURCE CLARIFICATION: seed is NOT the TRF motif.
            # It comes from the Step-02 lag-scan period, selecting an actual
            # genome window with high identity to its neighbor one period away.
            # StringDecomposer then infers flexible boundaries using sequence
            # alignment; the seed interval does NOT enforce uniform copy length.
            # (seed alternatives: scripts/29_seed_ensemble_pilot.py)
            # ---- seed: best next-copy window; period is a PRIOR, not a fixed length
            arr = np.frombuffer(boot.encode(), dtype=np.uint8)
            # FIX v4.5.1: a seed may not contain N (assembly gaps inside the
            # scan interval produced all-N consensuses on the Gar run)
            isN = np.concatenate(([0], np.cumsum(arr == ord('N'))))
            bi, bv = -1, -1.0
            for i in range(0, min(len(boot) - 2 * per, 60000), max(1, per // 4)):
                if isN[i + 2 * per] - isN[i] > 0: continue
                v = (arr[i:i + per] == arr[i + per:i + 2 * per]).mean()
                if v > bv: bv, bi = v, i
            if bi < 0:
                out['status'] = 'SKIPPED_N_MASKED_NO_CLEAN_SEED'; return out
            library = [('C01', boot[bi:bi + per])]

            fa = os.path.join(tmp, 'arr.fa')
            with open(fa, 'w') as seq_out:
                seq_out.write(f'>{aid}\n{full}\n')

            prev_starts, prev_lab = None, None
            prev_libsim = float('nan')
            stable = 0
            conv = False
            for it in range(MAX_ITER):
                mf = os.path.join(tmp, 'lib.fa')
                with open(mf, 'w') as fh:
                    for nm, s in library: fh.write(f'>{nm}\n{s}\n')
                dec = ordered_valid_calls(
                    run_sd(fa, mf, os.path.join(tmp, 'sd'), second_best=True), len(full))

                starts, copies, labs, ids = [], [], [], []
                for f in dec:
                    idn = float(f[4])
                    s0, e0 = int(f[2]), int(f[3])
                    mono = f[1]; revf = mono.endswith("'")
                    cpy = full[s0:e0 + 1]
                    if revf: cpy = rc(cpy)
                    starts.append(s0); copies.append(cpy)
                    labs.append(mono.rstrip("'")); ids.append(idn)
                iv = sorted(zip(starts, [s + len(c) for s, c in zip(starts, copies)]))
                mg = []
                for s_, e_ in iv:
                    if mg and s_ <= mg[-1][1]: mg[-1][1] = max(mg[-1][1], e_)
                    else: mg.append([s_, e_])
                covered = sum(e_ - s_ for s_, e_ in mg)
                unres = max(0.0, 1.0 - covered / max(1, len(full)))

                shift = float('nan'); stab = float('nan')
                if prev_starts is not None:
                    dshift = nearest_boundary_shift(prev_starts, starts)
                    if dshift is not None: shift = dshift
                    stab = prev_libsim
                    # boundaries are what every downstream step consumes, so
                    # convergence is declared on boundary stability, not on the
                    # library churning between an equivalent set of classes
                    stable = stable + 1 if np.isfinite(shift) and shift < BOUNDARY_EPS else 0
                    conv = stable >= STABLE_ROUNDS
                out['iters'].append([aid, it, len(library), len(copies),
                                     round(float(np.mean(ids)), 2) if ids else 'NA',
                                     ('NA' if shift != shift else round(shift, 2)),
                                     ('NA' if stab != stab else round(stab, 4)),
                                     round(unres, 4), int(conv)])
                if conv or it == MAX_ITER - 1:
                    break

                # ---- rebuild the library: cluster copies into monomer CLASSES
                good = [c for c, i2 in zip(copies, ids)
                        if i2 >= 60 and abs(len(c) - per) <= LEN_TOL * per * 2
                        and c.count('N') <= 0.05 * len(c)]
                if len(good) < MIN_CLASS:
                    break
                samp = good if len(good) <= CLUSTER_SAMPLE else \
                       [good[i] for i in np.linspace(0, len(good) - 1, CLUSTER_SAMPLE).astype(int)]
                groups = cluster_copies(samp)
                newlib = []
                for gi, g in enumerate(groups, 1):
                    if len(g) < MIN_CLASS: continue
                    mem = [samp[i] for i in g]
                    ref = mem[len(mem) // 2]
                    newlib.append((f'C{gi:02d}', consensus_of(mem, ref)))
                if not newlib:
                    break
                # library similarity: how much did the consensus set actually move?
                if len(newlib) == len(library):
                    sims = []
                    for (_, o), (_, n2) in zip(library, newlib):
                        d = edlib.align(o, n2 + n2, mode='HW')['editDistance']
                        sims.append(max(0.0, 1 - d / len(o)))
                    prev_libsim = float(np.mean(sims))
                else:
                    prev_libsim = 0.0
                # High library similarity alone is not evidence that the
                # source-coordinate boundaries are stable. Keep iterating
                # until the boundary criterion is satisfied or MAX_ITER hits.
                library = newlib
                prev_starts, prev_lab = starts, labs

            # ---- final pass, library fixed
            mf = os.path.join(tmp, 'lib_final.fa')
            with open(mf, 'w') as fh:
                for nm, s in library: fh.write(f'>{nm}\n{s}\n')
            dec = ordered_valid_calls(
                run_sd(fa, mf, os.path.join(tmp, 'sdf'), second_best=True), len(full))
            libd = dict(library)
            classcount = collections.Counter()
            for f in dec: classcount[f[1].rstrip("'")] += 1
            ntot = sum(classcount.values())

            for nm, cs in library:
                palin = 100 * (1 - edlib.align(cs, rc(cs) + rc(cs), mode='HW')['editDistance'] / len(cs))
                out['cons'].append([aid, chrom, st, en, nm, f'{aid}_{nm}_v{it+1}',
                                    len(cs), round(gc_(cs) * 100, 2),
                                    classcount.get(nm, 0),
                                    round(100 * classcount.get(nm, 0) / max(1, ntot), 2),
                                    round(palin, 2),
                                    'palindromic' if palin >= PALIN_ID else 'not_palindromic', cs])

            for mi, f in enumerate(dec):
                mono, s0, e0, idn = f[1], int(f[2]), int(f[3]), float(f[4])
                sb = float(f[6]) if len(f) > 6 and f[6] not in ('None', '') else -1.0
                strand = '-' if mono.endswith("'") else '+'
                cls = mono.rstrip("'")
                cons = libd[cls]
                raw = full[s0:e0 + 1]
                oriented = rc(raw) if strand == '-' else raw
                df = edlib.align(raw, cons, mode='NW')['editDistance']
                dr = edlib.align(rc(raw), cons, mode='NW')['editDistance']
                Lm = max(len(raw), len(cons))
                omargin = abs((1 - df / Lm) - (1 - dr / Lm)) * 100
                palin = 100 * (1 - edlib.align(cons, rc(cons) + rc(cons), mode='HW')['editDistance'] / len(cons))
                sstat = ('palindromic' if palin >= PALIN_ID else
                         'resolved' if omargin >= ORIENT_MARGIN else 'ambiguous')
                al = edlib.align(oriented, cons, task='path', mode='NW')
                ni = edlib.getNiceAlignment(al, oriented, cons)
                m = ni['matched_aligned']
                L = e0 - s0 + 1
                # v4.5.1 (SAT-017, AUD L81): per-copy rotation offset = start of
                # the best infix placement of the copy in the doubled consensus
                # (0 = in phase with the class consensus)
                loc = edlib.align(oriented, cons + cons, mode='HW', task='locations')['locations']
                phase = int(loc[0][0]) % len(cons) if loc else 'NA'
                trunc = 'complete' if abs(L - len(cons)) <= LEN_TOL * len(cons) else 'truncated'
                nN = oriented.count('N')
                atp = 100 * (oriented.count('A') + oriented.count('T')) / max(1, len(oriented) - nN)   # ACGT only
                out['mono'].append([
                    tag, ('C. gariepinus' if tag.startswith('Gar') else 'C. macrocephalus'),
                    tag, aid, f'{aid}_M{mi+1:05d}', chrom,
                    st + s0, st + e0 + 1, s0, e0 + 1, mi + 1, per, L,
                    strand, sstat, 'yes' if strand == '-' else 'no', phase,
                    'unassigned', 'unassigned', cls, f'{aid}_{cls}_v{it+1}',
                    round(idn, 2), (round(sb, 2) if sb >= 0 else 'NA'),
                    round(omargin, 2), (round(abs(idn - sb), 2) if sb >= 0 else 'NA'),
                    round(100 - idn, 2), m.count('|'), m.count('.'),
                    ni['target_aligned'].count('-'), ni['query_aligned'].count('-'),
                    trunc, 'yes' if (mi == 0 or mi == len(dec) - 1) else 'no',
                    round(atp, 2), round(100 - atp, 2),
                    hashlib.sha1(oriented.encode()).hexdigest()[:16],
                    'stringdecomposer', '1.1.2', run_id,
                    ('WARN_N_bases' if nN > 0.05 * len(oriented) else
                     'PASS' if idn >= 70 else 'WARN_low_identity'), oriented])
            out['status'] = 'DECOMPOSED' if dec else 'NO_MONOMER_CALLS'
            out['nlib'] = len(library)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
            gc.collect()
    except Exception as ex:                       # per-array failure, step continues
        out['status'] = 'FAILED:' + type(ex).__name__ + ':' + str(ex).replace('\t', ' ').replace('\n', ' ')[:120]
        out['cons'] = []; out['mono'] = []
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--genome', required=True)
    ap.add_argument('--arrays', required=True)
    ap.add_argument('--tag', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--min-array-bp', type=int, default=5000)
    a = ap.parse_args()
    run_id = os.environ.get('RUN_ID', 'na')

    with open(a.arrays) as array_handle:
        all_arrays = list(csv.DictReader(array_handle, delimiter='\t'))
    arrays = [r for r in all_arrays if int(r['array_length_bp']) >= a.min_array_bp]
    # Independent completion ledger, including calls ineligible for decomposition.
    status = {r['array_id']: ('SKIPPED_BELOW_MIN_ARRAY_BP'
              if int(r['array_length_bp']) < a.min_array_bp else 'NOT_PROCESSED')
              for r in all_arrays}
    arrays.sort(key=lambda r: (r['chrom'], int(r['start_0based'])))

    monp = os.path.join(a.out, '03_monomer_copies.tsv')
    conp = os.path.join(a.out, '03_array_consensus.tsv')
    itp  = os.path.join(a.out, '03_decomposition_iterations.tsv')
    ckpt = os.path.join(a.out, '.03_done_arrays.txt')
    # FIX v4.5.1 (audit #1): the per-array checkpoint is bound to the content
    # it depends on: the array intervals/periods it was computed for (columns
    # of the 02 table that this step consumes) and this script. A different
    # fingerprint means the old rows describe other loci -> start fresh and
    # park the old outputs as *.stale.
    fp_h = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read())
    for r in arrays:
        fp_h.update(f"{r['array_id']}\t{r['chrom']}\t{r['start_0based']}\t{r['end']}\t{r['period_mode_bp']}\n".encode())
    fp_tag = '#fingerprint=' + fp_h.hexdigest()
    done = set()
    if os.path.exists(ckpt):
        lines = [l.rstrip('\n') for l in open(ckpt)]
        head = lines[0] if lines else ''
        if head != fp_tag or not all(os.path.exists(p) for p in (monp, conp, itp)):
            why = 'inputs/script changed' if head != fp_tag else 'an output table is missing'
            print(f'    step-03 checkpoint invalid ({why}); previous outputs parked as *.stale, restarting')
            for p in (monp, conp, itp, ckpt, os.path.join(a.out, '03_array_status.tsv')):
                if os.path.exists(p): os.replace(p, p + '.stale')
        else:
            done = {l for l in lines[1:] if l}
            print(f'    resuming: {len(done)} arrays already decomposed')
    if done:
        # FIX v4.5.1 (audit #7): drop rows of arrays that were only partially
        # written before a crash (present in the tables, absent from the ckpt)
        for p in (monp, conp, itp):
            tmp = p + '.trim'
            with open(p, newline='') as src, open(tmp, 'w', newline='') as dst:
                hdr = src.readline(); dst.write(hdr)
                col = hdr.rstrip('\n').split('\t').index('array_id')
                for line in src:
                    if line.split('\t')[col] in done: dst.write(line)
            os.replace(tmp, p)
        for aid in done:
            if aid in status: status[aid] = 'DONE_PREVIOUS_CHECKPOINT'
    # a fresh step-03 result invalidates the step-04 family join of the old one
    if not done:
        for p in ('03_monomer_copies.pre_family_v3.tsv', '04_class_family_assignments.tsv', '04_satellite_families.tsv',
                  '04_array_dominant_family.tsv', '04_satellite_family_consensuses.fa', '04_consensus_identity_matrix.npy'):
            q = os.path.join(a.out, p)
            if os.path.exists(q): os.replace(q, q + '.stale')
    # previous per-array statuses survive a resume (FAILED arrays are retried)
    prevstat = os.path.join(a.out, '03_array_status.tsv')
    if done and os.path.exists(prevstat):
        for r in csv.DictReader(open(prevstat), delimiter='\t'):
            if r['array_id'] in done and r['status'].startswith(('DECOMPOSED', 'NO_MONOMER', 'SKIPPED', 'INVALID')):
                status[r['array_id']] = r['status']
    mode = 'a' if done else 'w'
    mfh = open(monp, mode, newline=''); cfh = open(conp, mode, newline='')
    ifh = open(itp, mode, newline='')
    mw = csv.writer(mfh, delimiter='\t', lineterminator='\n')
    cw = csv.writer(cfh, delimiter='\t', lineterminator='\n')
    iw = csv.writer(ifh, delimiter='\t', lineterminator='\n')
    if not done:
        mw.writerow(COLUMNS)
        cw.writerow(['array_id','chrom','start_0based','end','monomer_class',
                     'consensus_version','consensus_length_bp','consensus_GC_percent',
                     'n_copies_in_class','pct_of_array_copies',
                     'self_vs_revcomp_identity','palindromy_call','consensus_sequence'])
        iw.writerow(['array_id','iteration','n_library_monomers','n_copies',
                     'mean_identity','mean_boundary_shift_bp','library_similarity',
                     'unresolved_fraction','converged'])
    ck = open(ckpt, 'a')
    if not done: ck.write(fp_tag + '\n'); ck.flush()
    todo = [r for r in arrays if r['array_id'] not in done]
    nrows = narr = nfail = 0
    print(f'    {len(todo)} arrays to decompose on {max(1, a.threads)} worker(s)', flush=True)
    import contextlib
    # FIX v4.5.1 (audit #8): ProcessPoolExecutor raises BrokenProcessPool when
    # a worker is killed (OOM) instead of hanging like multiprocessing.Pool.imap
    from concurrent.futures import ProcessPoolExecutor
    if a.threads <= 1:                       # in-process: no pool, no pickling
        _worker_init(a.genome, a.tag, run_id)
        pool_cm = contextlib.nullcontext(None)
        def run_all(): return map(decompose_array, todo)
    else:
        pool_cm = ProcessPoolExecutor(a.threads, initializer=_worker_init, initargs=(a.genome, a.tag, run_id))
        def run_all(): return pool_cm.map(decompose_array, todo, chunksize=1)
    try:
        with pool_cm:
            for ai, res in enumerate(run_all(), 1):
                aid = res['aid']
                iw.writerows(res['iters']); cw.writerows(res['cons']); mw.writerows(res['mono'])
                nrows += len(res['mono'])
                if res['status'] in ('DECOMPOSED', 'NO_MONOMER_CALLS'): narr += 1
                if res['status'].startswith('FAILED'):
                    nfail += 1; print(f'    !! {aid}: {res["status"]}', flush=True)
                status[aid] = res['status']
                mfh.flush(); cfh.flush(); ifh.flush()
                # checkpoint only final outcomes; FAILED arrays are retried on resume
                if not res['status'].startswith('FAILED'):
                    ck.write(aid + '\n'); ck.flush()
                if ai % 5 == 0 or ai == len(todo):
                    print(f'    {ai}/{len(todo)} arrays, {nrows:,} copies, '
                          f'{res["nlib"]} classes on last array', flush=True)
    finally:
        mfh.close(); cfh.close(); ifh.close(); ck.close()
        # Status is a separate non-destructive product. Full input universe retained.
        with open(os.path.join(a.out, '03_array_status.tsv'), 'w', newline='') as sf:
            sw = csv.writer(sf, delimiter='\t', lineterminator='\n')
            sw.writerow(['array_id','chrom','start_0based','end','status','run_id'])
            for r in all_arrays:
                sw.writerow([r['array_id'],r['chrom'],r['start_0based'],r['end'],
                             status[r['array_id']],run_id])
    print(f'    wrote {nrows:,} monomer copies across {narr} arrays this pass; {nfail} arrays FAILED (see 03_array_status.tsv)')
    print(f'    per-iteration trace: {itp}')


if __name__ == '__main__':
    main()
