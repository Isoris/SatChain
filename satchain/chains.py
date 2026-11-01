"""chains.py  (SatChain v4.7, 2026-11-01; audit D02 / D23 / D17 / D12): ONE chain builder for 61, 61b, 64 and the benchmarks.

A chain = consecutive catalogue runs of one array and one family, in genomic order, linked when
  * the array key is the same: array_id_v46 when the column exists, else array_id with a 59 '~g<segment>' suffix stripped (D12),
  * the family key is the same: family_v46 when the column exists, else family (D02),
  * 0 <= gap <= max_gap (two runs that touch are linked, D23),
  * no assembly gap / contig join lies between the two runs (gaps.Gaps.between; D17), same chromosome.
Chains with fewer than min_runs runs are dropped; ONE default, MIN_RUNS = 2 (before 2026-11-01: 61 keyed on `family` with min_runs 4,
61b ignored the family with min_runs 2, 64 had no gap check, so the three scripts produced different junction sets from one catalogue).
Every link of a kept chain is one junction record: left run (OUT edge = left_end), right run (IN edge = right_start0), spacer_bp,
family pair.

usage:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from chains import build_chains
    chains, junctions, stats = build_chains(rows, max_gap=50000, min_runs=2, gaps=Gaps(bed), keep=None, satellite_only=True)
chain    : dict(chain_id, array_id, family, chrom, runs=[catalogue rows], start0, end, n_runs)
junction : dict(chain_id, array_id, family, chrom, gap_index, left_run, right_run, left_end, right_start0, spacer_bp, family_left,
                family_right, left=<row>, right=<row>)
stats    : dict(arrays, runs, chains, junctions, runs_in_chains, links_refused_at_gaps, links_refused_family, links_refused_distance,
                links_refused_chrom, chains_dropped_min_runs)
"""
import collections

MIN_RUNS = 2
MAX_GAP = 50000


def array_key(r):
    """array_id_v46 when present, else array_id without the 59 gap-segment suffix '~g<n>' (D12)."""
    v = r.get('array_id_v46')
    return v if v else r['array_id'].split('~g')[0]


def family_key(r):
    """family_v46 when the column exists (even when empty), else family (D02)."""
    return r['family_v46'] if 'family_v46' in r else r.get('family', '')


def _between(gaps, c, a, b):
    """gap or join in the stretch between two runs; closed = a join exactly on a run edge separates the runs (gaps.py HB07, 2026-11-01)"""
    try: return gaps.between(c, a, b, closed=True)
    except TypeError: return gaps.between(c, a, b)


def group_runs(rows, keep=None, satellite_only=True):
    """array key -> catalogue rows sorted by (chrom, start0). keep = set of array ids (array_id or array key)."""
    by = collections.defaultdict(list)
    for r in rows:
        if satellite_only and r.get('element_class', 'satellite') != 'satellite': continue
        if keep and r['array_id'] not in keep and array_key(r) not in keep: continue
        by[array_key(r)].append(r)
    for v in by.values(): v.sort(key=lambda r: (r['chrom'], int(r['start0'])))
    return by


def build_chains(rows, max_gap=MAX_GAP, min_runs=MIN_RUNS, gaps=None, keep=None, satellite_only=True):
    by = group_runs(rows, keep, satellite_only)
    chains, junctions = [], []; st = collections.Counter()
    for aid, rs in by.items():
        st['arrays'] += 1; st['runs'] += len(rs); cur = [rs[0]]

        def flush(cur):
            if len(cur) < min_runs: st['chains_dropped_min_runs'] += 1; return
            cid = f'{aid}#{len(chains)}'
            ch = dict(chain_id=cid, array_id=aid, family=family_key(cur[0]), chrom=cur[0]['chrom'], runs=cur,
                      start0=int(cur[0]['start0']), end=int(cur[-1]['end']), n_runs=len(cur))
            chains.append(ch)
            for i, (p, r) in enumerate(zip(cur, cur[1:])):
                junctions.append(dict(chain_id=cid, array_id=aid, family=ch['family'], chrom=ch['chrom'], gap_index=i,
                                      left_run=p['run_id'], right_run=r['run_id'], left_end=int(p['end']), right_start0=int(r['start0']),
                                      spacer_bp=int(r['start0']) - int(p['end']), family_left=family_key(p), family_right=family_key(r),
                                      left=p, right=r))
        for r in rs[1:]:
            p = cur[-1]; gap = int(r['start0']) - int(p['end']); ok = True
            if r['chrom'] != p['chrom']: st['links_refused_chrom'] += 1; ok = False
            elif not (0 <= gap <= max_gap): st['links_refused_distance'] += 1; ok = False
            elif family_key(r) != family_key(p): st['links_refused_family'] += 1; ok = False
            elif gaps is not None and _between(gaps, r['chrom'], int(p['end']), int(r['start0'])): st['links_refused_at_gaps'] += 1; ok = False
            if ok: cur.append(r)
            else: flush(cur); cur = [r]
        flush(cur)
    st['chains'] = len(chains); st['junctions'] = len(junctions); st['runs_in_chains'] = sum(c['n_runs'] for c in chains)
    for k in ('links_refused_at_gaps', 'links_refused_family', 'links_refused_distance', 'links_refused_chrom', 'chains_dropped_min_runs'): st.setdefault(k, 0)
    return chains, junctions, dict(st)
