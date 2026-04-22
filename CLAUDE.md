# CLAUDE.md — svaba working notes

This file captures conventions, file landmarks, and open investigations for the
svaba SV/indel caller project so future sessions can pick up quickly. Update it
as understanding changes — it's the crash-safety net, not the README.

## Project at a glance

svaba is a structural variant (SV) and indel caller that uses local assembly +
read realignment to call variants from short-read BAMs. The canonical use case
is tumor/normal somatic calling, but it also supports germline and multi-sample
modes.

Top-level layout:

- `src/svaba/` — main C++ sources. Entry point is `run_svaba.cpp`; assembly,
  realignment, breakpoint scoring, VCF output, and postprocess all live here.
  File-naming convention is **PascalCase** — e.g. `SvabaOptions.cpp`,
  `SvabaOutputWriter.h`. A few intentional exceptions: `refilter.cpp/h`,
  `run_svaba.cpp`, `svaba.cpp`, `threadpool.h`, `tovcf.cpp`, `vcf.cpp/h`.
- `src/SGA/` — String Graph Assembler sources (vendored).
- `SeqLib/` — vendored htslib/bwa/fermi-lite wrapper used for BAM I/O,
  alignment, and assembly primitives. See "Build system" for how its flags
  get set.
- `bin/`, `build/` — build artifacts; don't edit by hand.
- `R/`, `viewer/`, `tracks/` — downstream analysis/visualization helpers.
- `tests/`, `example_data/` — test fixtures.
- `scripts/` — post-processing and utility shell helpers, all kept here
  (not at the repo root): `svaba_postprocess.sh`, `combine_blacklists.sh`,
  `extract_discordants.sh`, `filter_contig_supporting_reads.sh`,
  `r2c_for_contig.sh`, `sort_bps.sh`, `svaba_cloud.sh`. Profiling
  helpers (`memprof*.sh`) live under `opt/` (the user's ad-hoc tooling
  dir).
- `somlod_maxlod_analysis.html` — deep-dive writeup of the somatic log-odds
  scoring model. See "Statistical model" below.

Scripts that used to live here and are gone, in case you're looking for them:
`sort_output.sh` and `sort_and_deduplicate_bps.sh` were subsumed into the
unified `svaba_postprocess.sh`.

## Build system

`cmake -B build && cmake --build build` defaults to `CMAKE_BUILD_TYPE =
RelWithDebInfo`, which gives you **`-O2 -g -DNDEBUG -fno-omit-frame-pointer`**
for the svaba/SeqLib/SGA C++ code. Not O0, not O3.

Critical gotcha: **the `-O2` is hardcoded across all of the vendored
submodules**, and CMake's build-type doesn't reach them:

- `SeqLib/bwa/Makefile`: `CFLAGS = -g -Wall -Wno-unused-function -O2`
  (hardcoded; doesn't read the parent CMake). Makes `libbwa.a`, owner of ~38%
  of wall-time in real runs.
- `SeqLib/fermi-lite/Makefile`: `CFLAGS = -g -Wall -O2 -Wno-unused-function`
  (same story). Makes `libfml.a`, ~27% of wall-time.
- `SeqLib/CMakeLists.txt:10` does `set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")`
  — appends `-O2` to whatever you passed, and because later `-O` wins on
  gcc/clang, this silently defeats any top-level attempt to set `-O3` via
  `CMAKE_CXX_FLAGS` unless you pass it via `CMAKE_CXX_FLAGS_<CONFIG>` (which
  is appended after, and therefore wins).
- htslib is external; whoever built yours picked its flags.

To push the submodules to `-O3 -mcpu=native` (the single biggest free-perf
knob on Apple Silicon — 65% of compute is at -O2 generic right now):

```bash
make -C SeqLib/bwa        clean
make -C SeqLib/fermi-lite clean
make -C SeqLib/bwa        -j CFLAGS="-g -O3 -mcpu=native -fno-omit-frame-pointer -Wall -Wno-unused-function"
make -C SeqLib/fermi-lite -j CFLAGS="-g -O3 -mcpu=native -fno-omit-frame-pointer -Wall -Wno-unused-function"
cd build && make -j
```

Sanity-check: `make VERBOSE=1 2>&1 | grep -oE -- "-O[0-9sg]" | sort | uniq -c`.

## Statistical model — the files that matter

For anything related to how variants are scored (LOD, somatic vs germline,
error model), the two files to read first are:

- `src/svaba/SvabaModels.cpp` / `.h` — self-contained statistical primitives.
  - `LogLikelihood(d, a, f, e_fwd, e_rev)` (~lines 11-63) is the per-sample
    two-state error model: `p_ref = (1-f)(1-e_fwd) + f*e_rev` and
    `p_alt = f(1-e_rev) + (1-f)*e_fwd`. Returns log10 likelihood of observing
    `a` alt reads out of `d`. This is the primitive every higher-level score
    is built from.
  - `SomaticLOD(...)` (~lines 70-83) — public wrapper; forwards to the
    split-error implementation.
  - `SomaticLOD_withSplitErrors(...)` (~lines 86-189) — the active somatic
    model. Enumerates sub-hypotheses: `SOM_true`, `SOM_art`, `GERM_het`,
    `GERM_hom`, `GERM_art`, `GERM_shared`. Returns
    `log10( P(somatic) / P(any non-somatic) )`.
  - Line ~79: `const double eN_fwd = std::min(e_art_fwd, 0.005);` — hard cap
    on the normal-sample forward error rate. This cap is important: it means
    even in regions where the artifact model infers a high error rate, the
    normal sample is assumed to be clean. That's the knob you'd relax if you
    want somlod to be "artifact-aware" on the normal side.
  - Lines ~138-147: `GERM_shared` free-MLE branch. This is the sub-hypothesis
    that fits a single pooled allele fraction across tumor+normal. It exists
    to catch LOH (loss of heterozygosity) germline events where tumor VAF can
    be much higher than 0.5 while normal is still ~0.5. It is also the main
    reason `somlod` asymptotes slowly as tumor alt-support grows (see below).

- `src/svaba/BreakPoint.cpp` / `.h` — per-breakpoint scoring glue. Points of
  interest:
  - `BreakPoint::score_somatic()` at ~line 975 is the entry point that sets
    `LO_s` (the somatic LOD). Calls
    `SvabaModels::SomaticLOD(scaled_alt_n, a_cov_n, scaled_alt_t, a_cov_t, error_fwd, error_rev)`
    around lines 1028-1031.
  - `SampleInfo::modelSelection()` (~lines 1562-1674) computes the per-sample
    `LO = ll_alt - ll_err` at line ~1610. These are unnormalized log10
    likelihoods — absolute value is not meaningful, the difference is.
  - `max_lod` is computed at lines ~198-200 and ~1196-1198 as the max of
    `al.LO` across samples. This is the "is this artifact or not" score; it
    does grow with additional supporting reads because it compares a variant
    hypothesis to a pure-error hypothesis with no germline branch.
  - Lines ~1072-1073: the current INDEL somatic gate only tests `somlod` —
    it does not use `maxlod` as a co-gate. This is one of the levers in the
    proposed fixes.

## The somlod / maxlod investigation (still open)

**Problem statement.** Users observe that `maxlod` grows with tumor alt
support, but `somlod` barely moves once tumor alt-support gets high. For a
clean normal (`aN=0` or `aN=1`) you'd naively expect `somlod` to also keep
climbing as tumor evidence accumulates, but it asymptotes around ~9 for
`dN ≈ 30`.

**Diagnosis.**
- There is a real statistical ceiling on somlod roughly equal to
  `dN · log10(1/(1-fT_hat))`. With 30 normal reads you can only ever rule out
  a shared germline hypothesis so hard — it's a finite amount of evidence.
- The `GERM_shared` free-MLE branch makes `somlod` *approach* that ceiling
  slowly at sub-clonal tumor VAFs: when the MLE of a pooled AF lands in a
  germline-plausible band, it provides a strong non-somatic explanation that
  the somatic hypothesis has to beat.
- You cannot simply delete `GERM_shared` because doing so makes LOH germline
  cases (e.g. `aN=15/dN=30, aT=285/dT=300`) return `somlod ≈ +36`, a false
  positive. Verified numerically while writing the analysis.
- The `eN_fwd ≤ 0.005` cap means that in truly high-artifact regions, the
  model refuses to "forgive" 2-3 alt reads in the normal as artifacts, which
  is both good (prevents somatic false positives where normal is contaminated
  by real signal) and limiting (prevents somatic calls where the artifact
  model really does explain the normal reads).

**Proposed fixes (see `somlod_maxlod_analysis.html` for the full writeup).**
- Fix 1 — disjunction gate on `GERM_shared`: only let the free-MLE pooled
  branch influence somlod when `shared_is_germline_plausible || normal_evidence > 1.0`,
  where `normal_evidence = LL_N(dN, aN, f_n_mle) - LL_N(dN, aN, 0)`. This is
  error-rate aware and is the fix I'd land first.
- Fix 2 — loosen `eN_fwd` cap in known high-artifact regions.
- Fix 3 — BIC penalty on the free-MLE branch (1 free parameter costs
  ~`0.5 log10(dN+dT)` nats ≈ bits of evidence).
- Fix 4 — joint `maxlod` + `somlod` gate for INDELs (require both above
  threshold), since `maxlod` moves freely with tumor depth.
- Fix 5 — debug dump of sub-hypothesis LLs when `somlod` is within some
  epsilon of the gate, to make future failures diagnosable.

**Related fix landed already (SvABA2.0 v3 split-coverage gate):** the
old `both_split && homlen > 0` / `one_split && homlen == 0` branching
in `BreakPoint::splitCoverage` was removed. A read is now credited as
a split-supporter iff (a) its r2c alignment scores strictly higher
than its native alignment by at least a per-sample-prefix margin
(`T_R2C_MIN_MARGIN`=10% for tumor, `N_R2C_MIN_MARGIN`=0 strict for
normal — see `src/svaba/SvabaOptions.h`), and (b) it spans at least
one breakend on the contig. Long junction homology → r2c and native
tie → read doesn't credit either sample, which is the correct
conservative behavior (rather than the old "homology=0 one_split is
fine, homology>0 you need both_split" which nuked normal support
specifically when homology was long, biasing toward somatic calls).
The repeat_seq-length padding on the buffers is also gone — same
rationale, subsumed by the comparative score gate. See the user-
facing bp-id (v3 schema) work for how to trace a specific read's
current support attribution end-to-end.

**Important correctness notes (earned the hard way):**
- Don't propose `aN >= 2` style hard count gates without an error-rate
  adjustment. In a high-artifact region, 2-3 normal alt reads can be
  genuine artifacts, and gating on raw counts overcalls the case away. Always
  reason about `normal_evidence` (LL delta against `f=0`) instead of raw `aN`.
- A LL ratio reported by `LogLikelihood` is log10. Multiplying by ~3.32 gives
  bits. Per-read surprise is `(LL_alt - LL_ref) / d`. Absolute LL values have
  no meaning — always compare two hypotheses at the same data.
- The ~9 ceiling at `dN=30` is a real statistical bound; no reformulation of
  the somatic test can push above it. `GERM_shared` changes the *slope* of
  approach, not the asymptote.

## Postprocess pipeline

Everything post-`svaba run` is orchestrated by `scripts/svaba_postprocess.sh`.
It's the unified replacement for the old `sort_output.sh` +
`sort_and_deduplicate_bps.sh` pair, which no longer exist.

Five steps per invocation, all idempotent (missing inputs log a one-liner and
continue):

1. **Merge per-thread BAMs** — `${ID}.thread*.${suffix}.bam` → `${ID}.${suffix}.bam`
   for `discordant` / `weird` / `corrected`. Single-file inputs are moved;
   no-file inputs silently skipped.
2. **`svaba postprocess`** — the C++ subcommand in `src/svaba/SvabaPostprocess.cpp`.
   For each suffix (`weird`, `corrected`, `discordant`, `contigs`):
   - `samtools sort -@ per_job_threads -m MEM` (shell out — htslib doesn't
     expose its sort as a library call). **Auto-skipped** when the BAM
     already declares `@HD SO:coordinate` — `isCoordinateSorted()`
     inspects the header via `readHeaderOnly()` and logs
     "already coordinate-sorted; skipping sort" so reruns are a no-op.
   - Native streaming dedup (only for `weird`/`corrected`/`discordant`):
     reads BAM via SeqLib::BamReader, collapses exact (qname, flag)
     duplicates at each locus, and **unions their `bi:Z` / `bz:Z` comma-token
     lists** so alt-supporting-contig evidence isn't lost when the same read
     got emitted by two overlapping assembly windows. Key function:
     `mergeCommaTokens` — boundary-aware union, mirrors
     `SvabaOutputWriter::stamp_tag`.
     Both the reader and writer have the htslib BGZF thread pool enabled
     via a new `SeqLib::BamReader::SetThreads(int)` /
     `SeqLib::BamWriter::SetThreads(int)` API that calls
     `hts_set_threads(fp, n)`. `streamDedup` accepts a `threads` parameter
     (wired from the full postprocess budget, not the per-suffix slice —
     see two-phase note below) and applies it on both sides.
     Without this pool, BGZF decompress + compress is single-threaded
     and dominates wall time (40 GB BAM → ~2 hours). With `-t 4..8`
     this typically drops to 25–40 min, a 3–5× speedup.

     **Buckets-clear gotcha** (landed alongside the thread-pool fix):
     the per-locus `idx_by_key` is a `std::unordered_map<std::string,
     size_t>`. `unordered_map::clear()` is **O(bucket_count)** — it
     memsets the entire bucket array even when size is 0, and
     `bucket_count` only grows, never shrinks. One pileup locus
     (centromere / simple repeat / HLA) inflates buckets to 100k+
     for the rest of the BAM, and every subsequent locus transition
     (hundreds of millions of them) pays that memset. Pre-fix perf
     showed 95% of main-thread CPU going to
     `__memset_evex_unaligned_erms` called from `unordered_map::clear`.
     `flushLocus` now swaps with a fresh small map when
     `bucket_count() > 256`, paying the inflated-map destructor cost
     ONCE per pileup exit rather than once per locus. Lesson: never
     trust `unordered_map::clear()` in a transient reuse pattern where
     the map briefly inflates.

     **Two-phase driver.** `svaba postprocess` runs in two phases so the
     thread budget lands where it actually helps:

       Phase 1 (PARALLEL): samtools sort across all active suffixes
         concurrently, each worker with `o.threads / n_active` threads.
         Sort is disk+CPU bound and scales linearly across files.
       Phase 2 (SERIAL): dedup + reheader + index, one BAM at a time,
         each with the full `o.threads` as its BGZF pool. Serial here
         is deliberate — running dedup in parallel across suffixes
         would oversubscribe (each BAM needs its own read+write BGZF
         pool), and BGZF parallelism has diminishing returns so
         `4 workers × 2 threads` is worse wall-clock than
         `1 worker × 8 threads` iterated four times.

     **Idempotency.** Every phase has its own auto-skip so rerunning
     the pipeline on an already-finished BAM is essentially instant:

     - Phase 1 sort: `isCoordinateSorted()` inspects `@HD SO:coordinate`
       and bypasses the sort when already done.
     - Phase 2 dedup: `hasSvabaPostprocessPg()` scans the @PG chain for
       any `svaba_postprocess` (or uniquified `.1`, `.2` variants); if
       present, dedup AND the subsequent reheader step are both skipped
       (only the `.bai` is rebuilt, which is cheap and covers the
       missing-index case). The first successful `streamDedup` stamps
       the PG line, so a second `svaba_postprocess.sh` run on the same
       outputs no-ops almost entirely.
     - The shell-layer merge step (`scripts/svaba_postprocess.sh` step 1)
       is already a no-op when per-thread `.thread*.bam` files aren't
       present (nothing to merge), so all three steps compose naturally.
   - **@PG stamp**: writes an `@PG ID:svaba_postprocess PN:svaba VN:<ver> CL:<argv> PP:<prev chain tail>`
     line into the output header. For dedup-eligible suffixes this is free
     (done in the writer during dedup); for `contigs` (no dedup) it's done
     via `samtools reheader`. ID is auto-uniquified if `svaba postprocess`
     has been run before on the same BAM.
   - **BAI index** via `sam_index_build(fn, 0)` directly from htslib — one
     less subprocess.
   - Intermediate filenames use a `.postprocess.*.tmp.bam` suffix and are
     renamed over on success / unlinked on failure. End state per suffix is
     exactly `${ID}.${suffix}.bam(.bai)`; no `.sorted` / `.deduped` flotsam.
3. **Sort + dedup + PASS-filter `bps.txt.gz`** — external-memory sort via
   GNU sort (`gsort` on macOS from homebrew coreutils). Sort keys:
   `chr1(V), pos1(n), strand1, chr2(V), pos2(n), strand2, maxlod(gr)` —
   maxlod descending so the "best SV per junction" survives the dedup.
   Produces three files: `.bps.sorted.txt.gz`, `.bps.sorted.dedup.txt.gz`
   (one row per unique breakpoint pair), `.bps.sorted.dedup.pass.txt.gz`
   (col 32 == "PASS" only).
   - Column positions hard-coded: `col 30 = cname (contig_and_region)`,
     `col 32 = confidence`, `col 38 = maxlod`. These come from
     `BreakPoint::toFileString` — change there and the script breaks.
4. **Filter `r2c.txt.gz` to PASS contigs (+ PASS-somatic subset)** —
   resolves PASS cnames from `bps.txt.gz` (col 32 == "PASS") and the
   PASS-somatic subset (also col 37, `somlod`, >= 1) in one pass. Then
   one pass over `r2c.txt.gz` writes two outputs via awk pipe-to-command
   (gzip compressors run in parallel as child processes):
     - `${ID}.r2c.pass.txt.gz`           (all PASS)
     - `${ID}.r2c.pass.somatic.txt.gz`   (PASS ∩ somlod >= 1)
   Both are cname-keyed; contig and read rows survive together since
   they share col 2. Either is suitable for `bps_explorer.html`'s
   alignments sub-panel; the somatic subset is the lighter load when
   you only care about the tumor-specific calls.
5. **Optional split-by-source** — `--split-by-source` (or env
   `SPLIT_BY_SOURCE=1`) demuxes the deduped BAMs by the first 4 chars of
   each QNAME into `${ID}.${suffix}.${prefix}.bam`.

CLI: `scripts/svaba_postprocess.sh -t THREADS -m MEM [other flags] <ID>`. Flags:
`-t/--threads`, `-m/--mem`, `--sort-buffer`, `--split-by-source`,
`--input-dir`, `--output-dir`, `--svaba`, `--keep-tmp`, `--skip-bam`,
`--skip-dedup`, `--skip-bps`, `--skip-r2c`, `--skip-split`, `-h/--help`.
`--skip-dedup` maps to `--sort-only` on the C++ CLI: keeps sort + @PG +
index but skips the dedup pass. Combined with the C++'s auto-skip of
sort when the BAM header already declares `@HD SO:coordinate`
(`isCoordinateSorted()` in `SvabaPostprocess.cpp`), a rerun on
already-postprocessed files is effectively instant — useful for
refreshing just the index or PG stamp. Env fallbacks for
backward compat: `THREADS`, `MEM`, `BUFFER_SIZE`, `SVABA`, `SAM`,
`INPUT_DIR`, `OUTPUT_DIR`, `SPLIT_BY_SOURCE`, `KEEP_TMP`.

Gotcha: `zcat` on macOS is BSD zcat (looks for `.Z`), not GNU. This script
uses `gzip -dc` everywhere for portability — if you add a new decompression
step, do the same or it'll fail on Mac.

## `svaba tovcf` subcommand

Standalone converter: pre-deduplicated `bps.txt.gz` → two VCFv4.5 files
(`${id}.sv.vcf.gz` + `${id}.indel.vcf.gz`). Lives in `src/svaba/tovcf.cpp`;
dispatch wired into `src/svaba/svaba.cpp`. Runs the VCFFile parser with
`skip_dedup=true` because the postprocess pipeline has already sorted
and deduplicated the bps.txt.gz upstream — the internal interval-tree
dedup would just be wasted work and can spuriously dup things via its
SV-pileup blacklist.

Design calls made (see `docs/vcf-design-decisions.md` if created; for
now this is the reference):

- **VCF spec:** declares `##fileformat=VCFv4.5` (latest formal spec,
  Oct 2024). Backwards-compatible at the record level with anything
  that accepts 4.2+.
- **File structure:** one SV VCF + one indel VCF per sample-set. Every
  record (somatic + germline) goes into both files; somatic rows carry
  the `SOMATIC` flag INFO so `bcftools view -f .,PASS -i 'INFO/SOMATIC'`
  peels them off cleanly. Replaces the older 4-file split
  (sv.somatic / sv.germline / indel.somatic / indel.germline).
- **SV representation:** intrachrom events with unambiguous orientation
  become single-record symbolic alleles — `+/+` or `-/-` → `<DEL>`,
  `-/+` → `<DUP>`, `+/-` → `<INV>`. Inter-chrom and anything else falls
  through to paired BND records with mate-bracket notation. Override
  with `--always-bnd` to force BND for every SV (useful if downstream
  tooling gets confused by symbolic alleles). Classification lives in
  `classify_symbolic_kind()` in `vcf.cpp`.
- **EVENT grouping:** both BND records of a pair get `EVENT=<bp_id>`
  (taken from col 52 of bps.txt.gz, the v3 per-BP identifier). This
  uses the same namespace as `r2c.txt.gz`'s `split_bps`/`disc_bps` and
  the BAM `bi:Z` tag, so a user can follow a single variant across all
  svaba outputs with one key.
- **QUAL column:** defaults to `.` (missing). QUAL was historically the
  Phred of `Σ per-sample LO`, which is strongly correlated with
  `INFO/MAXLOD` but can mislead users into filtering on QUAL when they
  should be filtering on SOMLOD / MAXLOD / FILTER. `.` is VCF-spec-valid
  missing and makes downstream tools fall through to the INFO fields
  where the canonical scores live. Override with `--qual maxlod`
  (writes `round(10 * maxlod)` capped at 99) or `--qual sum` (legacy).
- **SVCLAIM:** per VCF 4.5. `J` for pure assembly-only SVs, `DJ` for
  ASDIS (both assembly + discordant evidence), `D` for discordant-only.
  svaba is junction-native so most calls end up `J` or `DJ`.

Flag surface:
```
svaba tovcf -i BPS.txt.gz -b BAM -a ID [options]
  --sv-out FILE           override SV path (default ${ID}.sv.vcf.gz)
  --indel-out FILE        override indel path (default ${ID}.indel.vcf.gz)
  --plain                 write plain .vcf (no bgzip)
  --always-bnd            force BND for every SV
  --qual MODE             missing (default) | maxlod | sum
  --include-nonpass       include FILTER != PASS records
  --dedup                 re-run legacy interval-tree dedup
  -v / --verbose
  -h / --help
```

Gotcha: the BAM is required only for the chromosome name/length table
(used by contig `##contig` lines in the VCF header and by the
BreakPoint parser to turn chrom-name strings into chr IDs). No reads are
actually read from it — any BAM that shares the reference is fine.

Not-yet-done on this subcommand: bgzip-proper (current `--plain=false`
output is plain gzip, which bcftools accepts but tabix doesn't index
correctly). For now, pipe through `bcftools sort -Oz` + `tabix -p vcf`
after. If this becomes painful, revisit with htslib's `hts_open` +
`vcf_write` path instead of ogzstream.

## BreakPoint IDs (v3 schema)

Every BP gets a unique stable identifier of the form `bpTTTNNNNNNNN`
where `TTT` is the 3-digit zero-padded worker thread ID and
`NNNNNNNN` is an 8-digit zero-padded per-thread running counter.
Generation is lock-free: the counter lives on `svabaThreadUnit` (see
`next_bp_id()` in `src/svaba/SvabaThreadUnit.h`). Assignment happens
exactly once per BP in `SvabaRegionProcessor::process()` right before
the pointer is pushed into `unit.m_bps`. Because `AlignedContig`
holds BPs via `shared_ptr`, the id mutation is immediately visible
through every reference (global, multi-map, indel breaks).

The id lands as the 52nd core column of `bps.txt.gz` (right after
`flipped`, before per-sample blocks — this is the v3 schema; v2 had
51 cols, LEGACY had 41). It's also carried into `r2c.txt.gz` (see
next section) so a read's support attribution is unambiguously
linked to the exact BP row in bps.txt, eliminating the old "which BP
on this contig did this read actually support?" puzzle.

`svaba_bps_cols` (from `scripts/svaba_local_function.sh`) documents
the full layout; column 52 is the bp_id field.

**BAM aux tags `bi:Z` / `bz:Z` (v3).** The two aux tags svaba stamps
on weird/discordant/corrected BAM outputs now live in *different*
identifier namespaces — choose the right one for the join you want:

- `bi:Z` — comma-joined list of **bp_ids** this read supports as
  ALT. Matches the per-BP resolution of `r2c.txt.gz`'s `split_bps`
  / `disc_bps` columns and `bps.txt.gz`'s col 52. Pre-v3 this
  carried cnames (contig-level), which couldn't disambiguate a
  contig that hosted multiple BPs (global + multi + indel). To pull
  every ALT-supporting read for a specific variant row:
  `samtools view corrected.bam | grep bi:Z:bp00100000042`. The tag
  is populated in `SvabaRegionProcessor::process` at the BP
  finalization loop (`tag_with_bp_id` lambda), mirrored into
  `svabaThreadUnit::alt_bp_ids_by_name` for the corrected-BAM
  restamp path in `SvabaOutputWriter::writeUnit`.
- `bz:Z` — comma-joined list of **cnames** this read r2c-aligned
  to. Stays cname-keyed because "which contigs did this read align
  to" is a contig-level concept. Populated inside the r2c loop in
  `SvabaRegionProcessor` alongside `svabaRead::AddR2C(...)`,
  mirrored into `svabaThreadUnit::all_cnames_by_name`. Superset of
  `bi:Z` in the sense that every ALT-supporter r2c-aligned to its
  contig, but uses a different key namespace, so join back to
  bps.txt by col 30 (cname) here vs col 52 (bp_id) for `bi:Z`.

`scripts/r2c_for_contig.sh` defaults to `bz:Z` to pull all r2c'd
reads for a contig; set `TAG=bi` and pass a bp_id instead of a
cname if you want the ALT-supporter subset for a specific variant.

## r2c TSV format (re-plot-able alignments)

`${ID}.r2c.txt.gz` is the structured, re-plot-able alignment dump. It
replaces the old pre-rendered `${ID}.alignments.txt.gz` ASCII output,
which has been removed entirely — anything that lived implicitly in the
ASCII art (fragment orientation, leading/trailing soft-clip bases, etc.)
is now available as an explicit field in the TSV, and
`viewer/r2c_explorer.html` re-plots it in-browser on demand.

**Per-thread emission + postprocess merge.** Each svaba worker writes its
own `${ID}.thread${N}.r2c.txt.gz` during the run (stream lives in
`svabaThreadUnit::r2c_out_`; opened in the ctor, closed in the dtor,
gated on `opts.dump_alignments`). The write happens in
`SvabaOutputWriter::writeUnit` **before** `writeMutex_` is acquired —
each thread's gzip stream is independent, so deflate runs in parallel
across all workers. The first worker (threadId == 1; the worker pool in
`threadpool.h` numbers threads 1..N, so there is no thread 0) writes the
column-header line on open; other threads start with data.
`scripts/svaba_postprocess.sh` step 1 merges the per-thread files via
`cat`: gzip is concatenation-safe per RFC 1952, and the postprocess step
numerically sorts `.threadN.r2c.txt.gz` so thread 1 is first in the cat,
which means the merged file has exactly one header at the top. This is
the same architectural pattern as the per-thread BAM writers; for the
rationale, see "Perf notes" — in short, a mutex-shared gzip stream
serializes `deflate()` across all threads, losing ~15/16 of compression
parallelism at `-p 16`.

Gotcha: the header-writing branch is keyed to `threadId == 1` in
`SvabaThreadUnit::svabaThreadUnit`. Older revisions checked
`threadId == 0`; that never fired because the pool hands out 1..N,
which silently produced headerless `r2c.txt.gz` files. If you ever
change the worker-numbering base, that line must change with it.

Schema (documented by `AlignedContig::r2cTsvHeader()` in
`src/svaba/AlignedContig.cpp`) — 21 tab-separated columns, record-type
discriminated:

```
record_type  contig_name  contig_len  contig_seq  frags  bps  n_reads
             # contig-only fields above — NA on read rows

             read_id  read_chrom  read_pos  read_flag  r2c_cigar  r2c_start
             r2c_rc  r2c_nm  support  split_bps  disc_bps  r2c_score
             native_score  read_seq
             # read-only fields — NA on the contig row
```

`split_bps` / `disc_bps` are comma-joined `bp_id` lists — the unambiguous
per-BP attribution of each read. The categorical `support` field
(`split` / `disc` / `both` / `none`) is derived from these (whichever is
non-NA) and kept for grep-friendliness. The `bps` field on the contig
row also carries each BP's id as the 2nd subfield so viewers can join
back without a second file. On older r2c files without these columns
(pre-v3 emitter) the viewer falls back to the categorical support only.

One `contig` row per variant-bearing contig, followed by one `read` row per
r2c-aligned read. `contig_name` is the shared join key (same value as
bps.txt's col-30 `cname`), so grouping reads back to their contig is O(n).

Format details:
- `frags`: `|`-separated per fragment; within a frag, `:`-separated
  `chr:pos:strand:cigar:mapq:cpos_break1:cpos_break2:gpos_break1:gpos_break2:flipped`.
- `bps`: `|`-separated; within a bp, `:`-separated
  `kind:chr1:pos1:strand1:chr2:pos2:strand2:span:insertion`. `kind ∈ {global, multi, indel}`.
  Insertion is `.` when absent so the field count stays fixed at 9.
- All cell values are TSV-escaped (tab/CR/LF → space) via `r2cEscape`.
- Per-sample support classification (`split`/`disc`/`both`/`none`) and the
  scores are computed **identically** to the ASCII emitter — same
  `r2c_score > native_score` gate, same `corrected_native_cig` precedence,
  same `split_supporters`/`disc_supporters` sets. The two emitters share
  enough code to prevent drift.

## Blacklist-aware region pruning

`run_svaba.cpp` drops any queued region that is 100% covered by the
blacklist BEFORE the region is sent to a worker thread. Lives right after
`loader.countJobs(regionsToRun)`, before `sc.total_regions_to_process` is
set. Rule: for each `r` in `regionsToRun`, if
`sc.blacklist.FindOverlapWidth(r, true) >= r.Width()`, drop it.

This was a measured win. Previously each fully-blacklisted chunk (e.g. a
25 kb window on a decoy/alt/random contig) paid the full pipeline cost:
`QueryRegion` on the ref, `walker->SetRegion` + `readBam` (which
decompresses every BGZF block overlapping the region, parses each
bam1_t, allocates an `svabaRead`) — only to have `sc.blacklist.CountOverlaps`
drop every read at `SvabaBamWalker.cpp:181-182`. Now those regions never
hit a thread.

Safe because `sc.blacklist` has had `MergeOverlappingIntervals()` +
`CreateTreeMap()` called, so `FindOverlapWidth` can't double-count and
wrongly drop a partially-callable region.

The per-read and per-BP blacklist checks at `SvabaRegionProcessor.cpp:74,
818` still run for regions that **partially** overlap — this prune only
short-circuits the 100%-covered case.

Pruned regions don't get a `runtime.txt` row. That's intentional and
actually makes the runtime file cleaner (only regions that did work).

## Options surface

`src/svaba/SvabaOptions.h`:

- `dump_weird_reads` is **compile-time-only** (`static constexpr bool`).
  Flip the literal in the header and rebuild to enable; there is
  deliberately no CLI path. Because it's `static constexpr false` by
  default, the compiler dead-code-eliminates every
  `if (sc.opts.dump_weird_reads) { ... }` branch at -O2+.
- `dump_discordant_reads`, `dump_corrected_reads`, and `dump_alignments`
  default **false** and are all enabled together by the single
  `--dump-reads` runtime flag (switch-case 1800 in `SvabaOptions.cpp`).
  They control, respectively:
    - `${ID}.discordant.bam`
    - `${ID}.corrected.bam`
    - `${ID}.r2c.txt.gz` (per-thread; merged by the postprocess step)
  The fields are kept separate so individual callsites can key off
  their own concern (e.g. `svabaThreadUnit` gates `r2c_out_` opening
  on `dump_alignments` only), but there's no runtime path to toggle
  them individually — all three flip as a unit under `--dump-reads`.
- Without `--dump-reads`, svaba produces the lean output set only:
  `bps.txt.gz`, VCFs, `contigs.bam`, `runtime.txt`, `discordant.txt.gz`
  (cluster-level, tiny). No per-thread `r2c.txt.gz`, no `corrected.bam`,
  no `discordant.bam`. This is the production default because the gated
  outputs can run to tens of gigabytes on deep samples.
- **`alignments.txt.gz` is gone.** The pre-rendered ASCII viewer output
  was replaced in full by `r2c.txt.gz` (same information, not
  pre-formatted). `AlignedContig::printToAlignmentsFile` and
  `BreakPoint::printDeletionMarksForAlignmentsFile` were removed; the
  surviving `AlignmentFragment::printToAlignmentsFile` is kept only
  because one `std::cerr` debug print in `BreakPoint.cpp` still calls
  it. `viewer/alignments_viewer.html` still exists and still works on
  old `.alignments.txt.gz` files from previous runs, but new runs don't
  produce that file — use the r2c sub-panel in `viewer/bps_explorer.html`
  instead.

## Viewer suite (`viewer/`)

Entry point is `viewer/index.html`, a card grid pointing at the
sub-viewers. All client-side, no server required.

- **`bps_explorer.html`** — primary viewer. Sortable table of bps rows,
  numeric filters (somlod/maxlod/qual/span/etc.), chip filters (counts
  are live — they reflect the current filter, not the full dataset),
  per-sample detail panel, small histograms for somlod/maxlod/span,
  IGV-navigation links in the IGV1/IGV2 columns (fires
  `fetch('http://localhost:60151/goto?locus=…', {mode:'no-cors'})`;
  requires IGV running with port 60151 enabled). r2c re-plot sub-panel
  was removed — that capability now lives in the standalone
  `r2c_explorer.html` below.
- **`r2c_explorer.html`** — standalone re-plotter for the structured
  r2c TSV (emitted by `svaba run --dump-reads`, or filtered to
  PASS / PASS-somatic by `scripts/svaba_postprocess.sh`). Upload an
  `.r2c.txt.gz` / `.r2c.pass.txt.gz` / `.r2c.pass.somatic.txt.gz`,
  type or pick a contig name in the search box (browser `<datalist>`
  autocomplete, capped at 5000 entries), and get the alignment plot
  rendered in-browser: 10bp ruler, contig sequence, per-BWA-fragment
  summaries with `>/<` orientation rendering, per-breakpoint summary
  rows, and per-read gap-expanded CIGAR rendering with lowercase
  leading/trailing soft-clip bases placed at `start - lead_clip_len`
  (exact mirror of the old `AlignedContig::printToAlignmentsFile`).
  Prev/Next/Random navigation, imperfect-only and per-support-kind
  toggles, click-to-hide source-prefix legend. Load-once, explore
  many contigs — no need to load a bps.txt.gz first.
- **`alignments_viewer.html`** — legacy ASCII-plot viewer for
  `alignments.txt.gz` outputs produced by svaba runs *before* the r2c
  migration. New runs no longer produce `alignments.txt.gz`;
  `r2c_explorer.html` is the replacement. Kept in the tree so historical
  outputs still render.
- **`runtime_explorer.html`** — explorer for `runtime.txt` (the per-region
  timing TSV produced by `svabaTimer::logRuntime`, schema in
  `SvabaUtils.cpp`). Sortable table, region-range filter,
  runtime/contigs/discordant filters, IGV-click on the region column,
  prominent runtime histogram (defaults to log10 because runtime
  distributions are always long-tailed — see "Perf notes"). 17-column
  schema hardcoded from `SvabaUtils.cpp::svabaTimer::header`.
- **`comparison.html`** — side-by-side of two bps runs.
- **`bps_viewer.html`** — legacy light-theme viewer, uses external
  `app.js` + `styles.css`.

## tracks/ and blacklists

`tracks/hg38.combined_blacklist.bed` is the one to feed `svaba run
--blacklist`. It's a **regeneratable artifact** — produced by
`scripts/combine_blacklists.sh` from the component files in the same dir. Don't
hand-edit it; edit a component file and re-run the script.

Components (as of last pass):
- `hg38.blacklist.sorted.bed` — ENCODE-style high-signal regions.
- `hg38.high_runtime.bed` — regions empirically slow to assemble.
- `hg38.manual.blacklist.bed` — ad-hoc bad-list, curated.
- `hg38.nonstd_chr.blacklist.bed` — full-contig entries for every
  chrUn/*_decoy/*_alt/*_random/chrEBV/HLA-* contig in the reference,
  generated from `tracks/chr` (a GRCh38 fasta-header dump).
- `hg38.rmsk.simple_repeat.bed` — UCSC RepeatMasker simple-repeat
  regions.

`scripts/combine_blacklists.sh` has three modes: plain concat (default),
`--merge` (sort + `bedtools merge` with distinct-label aggregation), and
`--clip GENOME` (clip interval ends to real contig length, because some
input BEDs use oversize-end sentinels like `end=250000000` that would
otherwise inflate covered-bp totals into the trillions — ask me how I
know). Preferred invocation for refreshing the combined blacklist:

```bash
./scripts/combine_blacklists.sh --merge --clip tracks/chr \
  -o tracks/hg38.combined_blacklist.bed \
  tracks/hg38.blacklist.sorted.bed \
  tracks/hg38.high_runtime.bed \
  tracks/hg38.manual.blacklist.bed \
  tracks/hg38.nonstd_chr.blacklist.bed \
  tracks/hg38.rmsk.simple_repeat.bed
```

## Perf notes (empirical, from profiling the HLA region on M3 Ultra)

From a `sample`-based CPU profile with `-p 16`:

- **96% of worker time is inside `SvabaRegionProcessor::process`**.
  Four buckets soak it up:
  - BWA alignment: **~38%** (`mem_align1_core`, `mem_chain`, seed lookups
    `bwt_sa`/`bwt_occ*`, `ksw_extend2`)
  - Fermi-lite assembly: **~27%** (`fml_assemble`, `fml_fmi2mag*`,
    `rld_rank2a`, `unitig_unidir`)
  - BFC error correction: **~17%** (`SeqLib::BFC::ErrorCorrect/Train`,
    `kmer_correct`, `bfc_*`). User-controllable via the ec flags.
  - BAM walking: **~9%** (`svabaBamWalker::readBam`)
- **IPC ≈ 3.0** on Apple Silicon P-cores (theoretical max ~8, typical
  real workloads 1.2–2.5). The code is already running near the silicon's
  sustained envelope; algorithmic rewrites have diminishing returns.
- **Voluntary context switches ≈ 42 across a 56s run** — i.e. essentially
  zero. Which means **there is no allocator contention, no mutex
  contention, no I/O blocking** to fix. jemalloc/mimalloc on macOS **lost
  to libmalloc by 5-10×** in an A/B test (DYLD interposition overhead +
  Apple's nanomalloc fast path for small allocs). Save the allocator swap
  for the Linux production box — on Linux it's often a 10–20% win; on
  macOS don't bother.
- **BWA and fermi-lite default to `n_threads = 1` internally**
  (`mem_opt_init` at `bwamem.c:100`, `fml_opt_init` at `bfc.c:21`) and
  svaba never mutates them. `kt_for` samples in the profile are the same
  worker thread running inline, not spawned threads. So svaba's `-p N`
  is 1:1 with worker threads; there's no multiplicative thread fan-out.

Thread-count guidance on this 20 P-core + 8 E-core box:
- `-p 20` is the default sweet spot. One worker per P-core; E-cores absorb
  macOS background work, jemalloc bg thread, etc.
- `-p 26` or `-p 28` is possible but the extra threads land on E-cores
  and become tail latency (E-core is ~50% the speed of a P-core). Usually
  a wash or slight loss vs `-p 20`.
- Load imbalance is the main remaining utilization gap. On a 3 Mb test
  region, utilization was ~63%; on a 10 Mb region it climbed to ~78%; on
  full chromosome / WGS expect 85–90%. Tail regions (HLA, centromere,
  high-coverage weird spots) dominate wall clock.

Build flag recommendation for production runs (see "Build system"):
`-O3 -mcpu=native` on bwa + fermi-lite + (via `CMAKE_CXX_FLAGS_RELWITHDEBINFO`)
the svaba C++. Measured 5–15% wall-time gain expected on top of the default
RelWithDebInfo build.

## Cloud scatter-gather (`scripts/svaba_cloud.sh`)

`svaba_cloud.sh` parallelizes a WGS run across multiple GCP VMs — one
per chromosome partition — sharing a single read-only persistent disk.
Each VM runs `svaba run -k <partition>` independently; outputs go to a
GCS bucket; an optional `--merge` step concatenates the per-partition
`bps.txt.gz` files and runs `svaba_postprocess.sh`.

Architecture rationale: svaba's bottleneck is BWA FM-index random
lookups, which are latency-bound and NUMA-hostile. Multi-socket servers
waste half their threads on cross-socket access. Small single-socket VMs
(e.g. `c3d-highcpu-30`, AMD Genoa, ~128 MB L3, no NUMA) give each
worker full-speed local memory. Horizontal scaling across VMs beats
vertical scaling to more threads on a big box.

On dual-socket Xeon servers (measured on 2×20-core Xeon @ 2.8 GHz),
jemalloc is the single biggest optimization — 37% wall-time reduction at
38 threads by eliminating glibc arena lock contention. NUMA pinning
(`numactl --cpunodebind --membind`) is second-order at ≤40 threads on
a 2-socket box but matters more on 4-socket or at higher thread counts.

Disk I/O is ~9% of wall time (sequential BAM streaming). NVMe is
unnecessary; standard persistent disk or even gcsfuse over a GCS bucket
is fine.

Interchromosomal SVs: both breakends get assembled independently by
whichever partition contains the discordant read pileup. The merge +
dedup step in `svaba_postprocess.sh` pairs them. No calls are lost.

## Conventions

- **File naming**: `src/svaba/Svaba*.{cpp,h}` for the svaba-specific files
  (PascalCase, first letter capitalized). Intentional lowercase / snake_case
  exceptions live at the top of `src/svaba/` (`refilter`, `run_svaba`,
  `svaba`, `threadpool`, `tovcf`, `vcf`). Don't introduce new lowercase-first
  names; if you rename something on a case-insensitive filesystem (macOS),
  do it via `git rm --cached <lower> && git add <Upper>` or git's `core.ignorecase`
  will hide the rename from the index.
- **C++ style**: snake_case methods inside ClassName, 2-space indent,
  header/impl split. Don't introduce new formatting unless asked.
- **Statistical code lives in `SvabaModels.*`**; breakpoint glue lives in
  `BreakPoint.*`. Keep statistical primitives in the models file, not
  inlined into BreakPoint.
- **LL/LOD values in this codebase are always log10**, not natural log.
- **`aN, dN, aT, dT`** = alt count / depth in normal and tumor; **`f`** =
  allele fraction; **`e_fwd`/`e_rev`** = forward/reverse error rates from
  the artifact model. These names are used consistently in the analysis
  HTML too.
- **Option codes** in `SvabaOptions.cpp::longOpts`: 1001-1099 = mode,
  1100s = assembly, 1200s = EC, 1300s = discordant, 1400s = filter,
  1500 = chunking, 1600s = bwa-mem tuning, 1700s = output/DBs, 1800 =
  dump-reads. Keep the ranges coherent when adding new options.

## Useful jump points

- Somatic LOD calc: `src/svaba/SvabaModels.cpp:86`
- Per-sample LO: `src/svaba/BreakPoint.cpp:1610`
- Somatic LOD entry: `src/svaba/BreakPoint.cpp:975`
- INDEL somatic gate: `src/svaba/BreakPoint.cpp:1072`
- Region-queue blacklist prune: `src/svaba/run_svaba.cpp` (right after
  `loader.countJobs(regionsToRun)`)
- Per-read blacklist filter: `src/svaba/SvabaBamWalker.cpp:181-182`
- r2c TSV emitter: `src/svaba/AlignedContig.cpp::printToR2CTsv` +
  `::r2cTsvHeader`
- Postprocess (C++): `src/svaba/SvabaPostprocess.cpp`
- Postprocess (shell orchestration): `scripts/svaba_postprocess.sh`
- `svaba tovcf` driver: `src/svaba/tovcf.cpp::runToVCF`
- VCF engine (parse + dedup + emit): `src/svaba/vcf.cpp` + `vcf.h`
- Symbolic SV classifier: `vcf.cpp::classify_symbolic_kind`
- Single-file VCF writers: `vcf.cpp::writeSvsSingleFile` + `writeIndelsSingleFile`
- Blacklist combiner: `scripts/combine_blacklists.sh`
- Runtime-file schema: `src/svaba/SvabaUtils.cpp::svabaTimer::header`
- Options parsing: `src/svaba/SvabaOptions.cpp::SvabaOptions::parse`
- Analysis writeup (somlod/maxlod): `somlod_maxlod_analysis.html`
