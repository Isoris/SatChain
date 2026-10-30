# SatChain: where a satellite array comes from, step by step

This document fixes what an "array" is in SatChain and traces every array from the first candidate to the final object. The definition
drifted between pipeline versions, and several results had to be withdrawn because of it (section 7). The rule adopted here is that
**sequence decides what is tandem, and topology decides what is one array**.

## 1. One word, one object

| Term | Definition | Step | Use it for |
|---|---|---|---|
| **candidate interval** | Region with tandem signal: a v4.5 scan interval (`<G>_SATnnnnn`) or a scan-57 candidate (`<G>_NEWnnnnn`) | 02 (v4.5); 57-58-58b | **Search space only.** Never satellite bp, never an array extent |
| **copy** | One StringDecomposer tile: C (complete, identity >= conf), W (weak) or J (junk) | 03 | Building blocks; a J tile is not satellite |
| **run** | Maximal stretch called array by the copy-chain HMM, trimmed to its first and last C copy, >= `min_c` C copies, passed the per-run sequence check | 56, val56, 59 | The tandem stretch of record; satellite bp = sum of run bp |
| **node (state)** | One element of the path along the chromosome: a copy class (k-mer verified, rank-named V1, V2, ... so arrays share one state space), another family's copy (OF), an uncertain copy (U), or an annotation state (TE:<superfamily>, 5S, tRNA, SINE, CDS/exon, GAP). A TE is one node, never masked | 53, 53b | The alphabet of topology |
| **topology** | The transition graph of a path: composition pi, transition counts, entropy rate H(next given current), memory MI, structure index S, second-order test | 53b | **What an array is** (section 3) |
| **segment** | Part of a run between accepted change points (two feature blocks agree, at least one label-free) | 53 | Sub-topologies of one array |
| **interruption** | Gap <= 50 kb between consecutive runs, classified TE / missed TE / other / unannotated | 59e | Interruption nodes of a chain |
| **chain** | Consecutive runs of one family, separated by interruptions, whose topologies are compatible (section 3, rule 2) | 61 (+ topology test) | Interrupted and composite arrays |
| **junction** | A run edge facing an interruption: OUT (satellite -> interruption) or IN (interruption -> satellite), array orientation | 61 | Boundary precision, structural ancestry |
| **cycle** | A recurring closed walk in the topology: self-loop (monomeric), cycle over n monomer nodes (HOR), cycle through an interruption node (composite unit) | 62, 65, 61 | The unit of an array at every scale |
| **family / lineage** | Majority family of confident copies; lineage A / B / centromeric minisatellite from copy-level k-mer clustering | 59, 55d | Naming, comparison |
| **putative centromeric array** | Largest run cluster of a chromosome (L1) | 31 / 59 | Location call, not function |

`array_id` in the catalogue is the identifier of the **candidate interval** a run came from: a container label, not an array.

## 2. The steps, in order

**Logic.** The array is inferred at the end, from the copy chain. Copies are observations, not array calls. They are computed only where
the genome-wide sweep finds tandem signal, so the sweep decides what can become an array at all:

```
genome-wide sweep (57-58-58b)  ->  copies = observations (03)  ->  array inference on the copy chain (56 HMM + val56)
      where to look                     what is there                 what is an array (runs)
                                                   ->  topology (53/53b) decides one array or several (section 3)
```

### Step 1 - genome-wide discovery sweep (where to look)
- **57 (the sweep).** Every chromosome, no annotation, no threshold. 2-kb windows (units 2 bp - ~1 kb) and 10-kb windows (units 1-5 kb),
  step half a window. Per window: k-mer (k = 10) next-recurrence histogram minus that of a base-shuffled copy of the window (excess), its
  lag, harmonic support at 2x the lag, recurring k-mer fraction. About 30 minutes on 3 cores for both haplotypes (1.86 Gb).
- **58 (calibration).** Null = ordinary windows outside known tandem repeats (TEs allowed); tiers p95 / p99 / p99.9 per scale; tandem windows
  with agreeing lags (+-10 % or 2x / 3x multiples) merged into candidate regions. Recall of validated runs at p99.9: 0.91 (Gar) / 0.86 (Mac).
- **58b (3-copy test).** A candidate is kept as `new_array` only with >= 3 consecutive copies (3-copy chain excess vs a dinucleotide-block
  shuffle); duplicated pairs and weak candidates are rejected; TRF agreement reported.
- **What it recovered.** 2,244 (Gar) / 2,138 (Mac) new runs, 1.90 / 1.93 Mb, outside every earlier array, mostly 7-39 bp minisatellites the
  old 40-1,000 bp lag scan never examined; TRF supports 70-76 % of the 61-1,000 bp candidates the test keeps and 5-10 % of the pairs it rejects.
- **Legacy seeds.** The v4.5 step-02 intervals (10-kb windows, lag identity >= 0.34 at lags 40 - max period, merged, reconciled with TRF)
  are kept as additional search regions. Their extents are never array extents: they are window unions that swallow adjacent TE and
  non-tandem sequence (satellite bp inside them 4-6x the real satellite: Gar 56.1 vs 13.5 Mb, Mac 73.6 vs 9.6 Mb), and the 0.34 cut lies on
  a cliff (real satellite at 0.30-0.34 lost).
- Guarantees: tandem signal somewhere in each region. Not: where an array starts or ends (steps 3 and 7 decide).

### Step 2 - copies (03)
- Iterative multi-consensus StringDecomposer: tile against a library of class consensuses, rebuild the library until boundaries stop moving;
  per copy: class, strand, identity, completeness.
- Not guaranteed: every tile is satellite. StringDecomposer tiles the whole interval; non-tandem sequence becomes low-identity tiles (35 %
  of MacSat03-centromere calls, 4,600 / 13,100, mixed strands, no similarity to neighbours). Downstream: runs (step 3) or an identity filter.

### Step 3 - array inference on the copy chain (56 -> val56)
- 56 (the Markov model of discovery): path of C / W / J / G observations; two-state HMM (array / background), Baum-Welch per subgenome;
  run = posterior >= `post` (swept 0.5 / 0.8 / 0.95), trimmed to first and last C, >= `min_c` C copies; short TE or degraded stretches stay
  inside when the chain continues. It replaces the fixed "N copies at X % identity" rules of v4.5 step 36.
- val56: mean 1-kb tandem signal of the run vs the background quantile (p80 / p90 / p95; p90 kept).
- What it recovered inside the old intervals: HMM-only bases are tandem (median 1-kb signal 0.157 / 0.116 vs 0.000 for in-array background);
  Gar adds 1.8 / 1.5 / 1.0 Mb at p80 / p90 / p95 and keeps 98-94 % of v4.5 bases. (Before SatChain, the v4.5.3 step-51 rescue had added
  0.97 / 1.52 Mb by re-decomposing uncovered 1-kb bins with tandem signal.)
- Guarantees: tandem in sequence, junk excluded (MacSat03 centromeres: 98 % of validated copies kept, 0 % of false ones).
  Not: base-pair edges (the HMM edge sits 9-154 bp inside the satellite edge; step 61 fixes it).

### Step 4 - catalogue (59, 59b, 59c, 59d)
- Per run: gates G1 GC / G2 identity / G3 strand coherence, family, lineage, unit class, tier; tandem gene arrays apart (59b); TE layer with
  intrinsic labels separated (59c); names (59d). Guarantees one row per run. Not: which rows are one array.

### Step 5 - topology of each run (53 -> 53b)
- 53: windows inside each run; feature blocks comp / kmer / period (label-free) and label / topo; change points accepted only when two
  blocks agree and one is label-free -> segments, each profiled (homogeneous monomeric, variant mixed, alternating, interleaved families,
  unresolved). Class merges: classes that no k-mer (k = 5-20) separates are merged before any topology is built.
- 53b: each run and each segment becomes a path of nodes (section 1); per unit: pi, transitions, entropy, entropy rate, MI, structure index S;
  second-order test (H1 - H2 vs simulations of the unit's own first-order chain, effect-size floor).
- Guarantees: a sequence-checked description of order inside every run, with annotations as nodes.

### Step 6 - interruptions (59e)
- Gaps <= 50 kb between consecutive runs of one `array_id`; EDTA labels intrinsic to the flanking runs ignored unless the feature is a whole
  element in the gap (>= 80 % inside, >= 500 bp); class TE insertion / missed TE / other / unannotated.

### Step 7 - chains, junctions, composite units (61 + topology compatibility)
- Chain candidates: consecutive runs of one family within one `array_id`, gaps <= `max_gap`, >= `min_runs`.
- Edge refinement: every edge facing an interruption re-aligned (local, rotation-aware) to the array's own dominant-class consensus ->
  exact edge, monomer phase, `delta_hmm`.
- Junction homology: junction contexts aligned to the medoid junction of their type -> flank identity (homologous or not), offset (bp).
- Stability score per chain and setting (flank 100/150/300 bp x cutoff 0.6/0.7/0.8 x tolerance +-10/25/50 bp): homologous share, precision,
  phase concordance, median |delta_hmm|. First genome-wide result: where junctions are homologous they are base-pair exact (median 100 %
  within +-10 bp); homologous share 40-67 % (Mac), 72-100 % (Gar).
- Composite units: IN junction -> next IN junction.

### Step 8 - cycles inside runs and chains (62, 65, 63, 64)
- 62: lag identity along the copy chain (shuffled-order null) and k-mer interval spectra (k = 6-31): size of the cycle.
- 65: monomer alphabet from Louvain communities over a swept identity threshold, hierarchical tandem-repeat mining of the node string,
  rotation merge, coverage ranking (adapted from HiCAT; every threshold reported; cross-array mode = one shared alphabet): the cycle itself.
- 63: exact-copy units per window, unit haplotypes, self-identity matrix: where cycles were copied recently.
- 64: one sequence profile per copy, run, interruption or extra object (length, GC, CpG, skews, entropy k = 1-6, k-mer spectra k = 2-8);
  separation tests with whole runs or arrays left out.

## 3. The array definition (decision): topology decides

1. **Tandem or not is a sequence question** (steps 1-3): only runs are tandem satellite.
2. **One array or several is a topology question.**
   - Inside a run: a change of topology confirmed by a label-free feature (step 5 segment boundary) separates sub-arrays; a run with one
     topology is one array.
   - Between runs of one family separated by an interruption: they are one array (a chain) when their topologies are compatible, i.e. the
     distance between their node compositions and transition matrices (Jensen-Shannon on pi and on the rows of the transition matrix, in the
     shared rank-named state space) is within the null obtained by splitting single uninterrupted runs at random points. Incompatible
     topologies = two arrays that happen to be adjacent, even inside one `array_id`.
3. **What kind of array is read from its cycles**, in one graph whose nodes are copy classes plus annotation / interruption nodes:
   - dominant self-loop (S high, MI low beyond self-loop) -> **monomeric array**;
   - a cycle over n monomer nodes that dominates the path (65 coverage, 62 lag-n peak) -> **HOR array** (n monomers, nested cycles allowed);
   - a recurring cycle that passes through an interruption node, with homologous junctions (61) -> **composite array** of composite units;
   - interruption nodes that do not recur and whose junctions are not homologous -> **interrupted array** with independent insertions.
4. **Boundaries** are reported at the resolution the evidence supports: copy-exact (run edge), base-pair (step 61 refined edge, with its
   stability score), or segment-level (step 5 change point, +-1 window).
5. **A candidate interval is never an array; the centromeric array is an L1 location call.**

## 4. Which number to quote

| Quantity | Quote | Never quote |
|---|---|---|
| Satellite bp | sum of run bp | candidate-interval / step-02 extents |
| Array length | run bp; for a chain: span, satellite bp and interruption bp separately | chain span as satellite |
| Array count | arrays per section 3 (runs, segments, chains: say which) | `array_id` count |
| Copies | confident copies inside runs | all 03 tiles |
| Boundary | step-61 refined edge with its stability score | HMM edge as an exact position |
| Unit | cycle length: monomer (median confident copy), HOR (65 / 62), composite unit (61) | consensus length |
| Differences between arrays | arrays as replicates | copies as replicates |

## 5. Orientation and phase
- Copy strand comes from step 03 (relative to the class consensus). A run's `reverse_fraction` can be inverted relative to the sequence on some
  chromosomes, so layouts and junction types use orientation from sequence (majority strand of confident copies, or step-61 orientation).
- Phase is measured on the array's own consensus.

## 6. Topology checks that protect the definition
- Classes are merged before topology when no k-mer separates them (53 class merges): a topology built on unseparable classes is noise.
- A topology change counts only with a label-free confirmation (53): wrong monomer calls cannot create array boundaries.
- Annotations enter as nodes, never as masks (53b): a TE inside an array is a step of the path, so a recurring TE node is visible as a cycle.
- Order beyond first order is claimed only when the second-order test passes against the array's own first-order simulations (53b).

## 7. Pitfalls we hit, and the rule that fixes each

| What went wrong | Effect | Rule |
|---|---|---|
| Satellite bp read from step-02 extents | satellite inflated 4-6x; artefact families (MacSat01/02, GarSat01) | satellite bp = runs only |
| 03 tiles used as monomers | 35 % false monomers in MacSat03 centromeres | runs or identity filter first |
| EDTA labels taken as TE content | satellites counted as TE | 59c / 59e intrinsic-label filter |
| A TE consensus coordinate used as unit anchor | extra boundaries 4-6 kb apart (LG04) | units from the satellite side (61) and cycles (section 3) |
| A TE consensus built from the array itself | "element" = whole repeat unit; its built-in satellite masked the real satellite | check a TE against its dispersed copies first |
| HMM edge used as exact boundary | 9-154 bp systematic shift, apparent drift | step-61 refinement |
| Family consensus used for phase / orientation | no alignment on divergent chromosome variants | array's own consensus |
| Windows mixing satellite and spacer (53 on composite chains) | every chain one segment, profiles averaged | per-object profiles (64) and interruption nodes |
| Copies from one run used as replicates | inflated separation | leave whole runs / arrays out; n = arrays |

## 8. Sweeps
56 posterior 0.5 / 0.8 / 0.95; val56 p80 / p90 / p95; 58 p95 / p99 / p99.9; 53 penalty and k-mer merge cutoffs; 61 flank x cutoff x tolerance;
62 k = 6-31, lags 1-40; 65 identity threshold 0.90-0.995; 63 K, window; 64 k = 2-8, leave-out block; topology compatibility: split-point null.

## 9. Topology is descriptive today, not yet a discovery sweep
Discovery uses the k-mer sweep (57) and the copy-chain HMM (56). The topology and second-order Markov test (53 / 53b) run on every catalogue
run after discovery (second-order structure: 88 / 459 Gar and 73 / 552 Mac arrays, L2). A genome-wide topology sweep (order structure
searched outside known runs, e.g. on 57 windows with a recurring lag) does not exist yet; it would be the next discovery module.

## 10. Open items
- Topology compatibility test between runs (section 3, rule 2) is defined here but not yet coded (inputs exist in 53b).
- 03 still writes junk tiles to its copy table; a flag in 03 would protect downstream users.
- Cross-array alphabet (65) to be run genome-wide; spacer-length variation inside composite units has no module yet.
- Every array, run, chain and unit is L1 (sequence); interpretations follow the project's L1 / L2 / L3 rules.
