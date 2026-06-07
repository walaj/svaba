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
- `R/`, `docs/`, `tracks/` — downstream analysis/visualization helpers. The
  HTML viewer suite now lives in `docs/` (it used to be `viewer/`; see the
  "Viewer suite" section).
- `tests/`, `example_data/` — test fixtures.
- `scripts/` — post-processing and utility shell helpers, all kept here
  (not at the repo root). Current set: `svaba_postprocess.sh`,
  `combine_blacklists.sh`, `mosdepth_lowmapq_blacklist.sh`,
  `extract_discordants.sh`, `filter_contig_supporting_reads.sh`,
  `r2c_for_contig.sh`, `extract_by_qname.sh`, `search_sequence.sh`,
  `sort_bps.sh`, `sort_and_dedupe_bps_old.sh` (legacy standalone sorter,
  kept for reference; the live path is `svaba postprocess`),
  `svaba_cloud.sh`, `gcloud_teardown.sh`, `update_svaba_image.sh`,
  `plot_learn.sh`, `svaba_explore.sh` (one-shot launcher for the bps
  explorer), `svaba_local_function.sh` (sourceable bash helpers, incl.
  `svaba_bps_cols`). Profiling helpers (`memprof*.sh`) live under `opt/`
  (the user's ad-hoc tooling dir).

The somatic-model deep-dive writeup is documented inline in the
"Statistical model" and "somlod / maxlod investigation" sections below;
the standalone `somlod_maxlod_analysis.html` it used to live in is no
longer in the tree.

Scripts that used to live here and are gone, in case you're looking for
them: `sort_output.sh` and `sort_and_deduplicate_bps.sh` were subsumed into
the postprocess pipeline (now the `svaba postprocess` subcommand);
`extract_pairs_by_seq.sh` was subsumed
into the `svaba extract-pairs` subcommand.

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
  - `SomaticLOD(...)` (~lines 70-83) — public wrapper; sets the per-sample
    forward-error caps (`eN_fwd`, `eT_fwd`) and forwards to the split-error
    implementation.
  - `SomaticLOD_withSplitErrors(...)` (~lines 86-335) — the active somatic
    model. Returns `log10( P(somatic) / P(any non-somatic) )` from the max
    over these sub-hypotheses:
      - SOMATIC (H1) = max of `SOM_true` (clean normal, tumor at its MLE VAF)
        and `SOM_art` (shared low-VAF artifact).
      - GERMLINE (H0) = max of `GERM_het` (0.5), `GERM_hom` (~1.0),
        `GERM_art` (low-VAF artifact), `GERM_shared` (pooled MLE; see below),
        and `GERM_free` (independent MLE; see below).
    (This note previously listed only `GERM_shared` on the germline side and
    capped the function at line ~189 — the model has since gained `GERM_free`
    and the `GERM_shared` shaping terms, and the function now runs to ~335.)
  - Line ~79: `const double eN_fwd = std::min(e_art_fwd, 0.005);` — hard cap
    on the normal-sample forward error rate. This cap is important: it means
    even in regions where the artifact model infers a high error rate, the
    normal sample is assumed to be clean. That's the knob you'd relax if you
    want somlod to be "artifact-aware" on the normal side. (`eT_fwd` is capped
    to `[1e-4, 0.02]` on line ~80.)
  - `GERM_shared` (~lines 143-231): pooled-MLE allele fraction across
    tumor+normal, to catch LOH germline events where tumor VAF can be much
    higher than 0.5 while normal is still ~0.5. No longer a plain free-MLE
    branch — the old 0.15 floor is gone and three shaping terms ride on top
    of the raw pooled-MLE log-likelihood: (1) a "cleanliness" penalty,
    `kCleanPenalty·exp(-3·excess_alt_N)`, that holds the branch back when the
    normal looks clean so a low-VAF tumor can still win as somatic; (2) a
    VAF-similarity bonus gated by `(1 - cleanliness)`; (3) a both-samples-
    low-VAF shared-artifact bonus. Still a main reason `somlod` approaches its
    ceiling slowly at sub-clonal tumor VAFs (see below).
  - `GERM_free` (~lines 238-277): NEW relative to the original writeup. MLEs
    `f_N` and `f_T` *independently* (not tied like `GERM_shared`), filling the
    gap where the normal carries low-but-real evidence while the tumor is
    near-clonal (LOH, clonal hematopoiesis, tumor-in-normal contamination,
    mosaic germline). Charges a BIC penalty `0.5·log10(dN+dT)` for the extra
    free parameter, and only activates when `aN > expected_errors_N + 1`
    (where `expected_errors_N = dN·eN_fwd`) so it doesn't degenerate to
    "SOM_true minus BIC" on every clean-normal event.

- `src/svaba/BreakPoint.cpp` / `.h` — per-breakpoint scoring glue. Points of
  interest:
  - `BreakPoint::score_somatic()` at ~line 1341 is the entry point that sets
    `LO_s` (the somatic LOD). Calls
    `SvabaModels::SomaticLOD(scaled_alt_n, a_cov_n, scaled_alt_t, a_cov_t, error_fwd, error_rev)`
    around line ~1394.
  - `SampleInfo::modelSelection()` (starts ~line 1954) computes the per-sample
    `LO = ll_alt - ll_err` at line ~2002. These are unnormalized log10
    likelihoods — absolute value is not meaningful, the difference is.
  - `max_lod` is computed at lines ~247-249 and ~1567-1569 as the max of
    `al.LO` across samples. This is the "is this artifact or not" score; it
    does grow with additional supporting reads because it compares a variant
    hypothesis to a pure-error hypothesis with no germline branch.
  - Lines ~1434-1439: the current INDEL somatic gate only tests `somlod`
    (`LO_s > cutoff`) — it does not use `maxlod` as a co-gate. This is one of
    the levers in the proposed fixes. (The SV gate is the sibling branch at
    ~1449.)

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

**Proposed fixes — current status.** (The standalone `somlod_maxlod_analysis.html`
writeup this list used to point at is no longer in the tree; this section is
now the reference, and `SvabaModels.cpp` has moved past the original proposal.)
- Fix 1 — error-rate-aware gate on the free-MLE pooled branch. **Partly
  landed.** Instead of a disjunction gate on `GERM_shared`, the model now
  (a) shapes `GERM_shared` with a cleanliness penalty keyed on
  `expected_errors_N = dN·eN_fwd`, and (b) added a separate `GERM_free`
  branch whose activation gate is `aN > expected_errors_N + 1` — both
  error-rate aware in the spirit of this fix.
- Fix 2 — loosen `eN_fwd` cap in known high-artifact regions. **Not landed**
  (cap still hard-coded at 0.005, `SvabaModels.cpp:79`).
- Fix 3 — BIC penalty on the free-MLE branch. **Landed** as `GERM_free`'s
  `kGermFreeBIC = 0.5·log10(dN+dT)`.
- Fix 4 — joint `maxlod` + `somlod` gate for INDELs. **Not landed** — the
  INDEL gate still tests only `LO_s` (`BreakPoint.cpp` ~1434-1439).
- Fix 5 — debug dump of sub-hypothesis LLs. **Present but disabled** — a full
  sub-hypothesis dump exists in `SomaticLOD_withSplitErrors` guarded by
  `if (false)` (`SvabaModels.cpp` ~284); flip it to debug a specific call.

**Related fix landed already (SvABA2.0 v3 split-coverage gate):** the
old `both_split && homlen > 0` / `one_split && homlen == 0` branching
in `BreakPoint::splitCoverage` was removed. A read is now credited as
a split-supporter iff (a) its r2c alignment scores strictly higher
than its native alignment (`r2c_score > native_score`, no percentage
margin — see `src/svaba/SvabaOptions.h`), and (b) it spans at least
one breakend on the contig. Long junction homology → r2c and native
tie → read doesn't credit either sample, which is the correct
conservative behavior (rather than the old "homology=0 one_split is
fine, homology>0 you need both_split" which nuked normal support
specifically when homology was long, biasing toward somatic calls).
The repeat_seq-length padding on the buffers is also gone — same
rationale, subsumed by the comparative score gate. See the user-
facing bp-id (v3 schema) work for how to trace a specific read's
current support attribution end-to-end.

**v3.1 fix — remove T_R2C_MIN_MARGIN (set to 0):**
the 10% `T_R2C_MIN_MARGIN` was killing all tumor alt-supporting reads
for small indels. A 1bp deletion on a 150bp read gives r2c=150 vs
native=143, a 4.9% improvement — mathematically impossible to clear
the 10% threshold. Traced via the `SVABA_TRACE_CONTIG` system on
contig `c_fermi_chr2_215869501_215894501_13C` (CIGAR `392M1D530M`):
every tumor read hit TP8 with r2c=150 vs threshold=157.3 → SKIP →
0 split support → LOWLOD → hasMinimal fail → variant dropped.

Fix: set `T_R2C_MIN_MARGIN` to 0.0 in `SvabaOptions.h` — same as
normal, strict greater-than only. Any percentage margin is inherently
read-length-dependent: a 1bp del gives 4.9% improvement on 150bp
reads but only 2.9% on 250bp reads, so any fixed percentage either
blocks long reads or is too loose for short reads. There's no
percentage that works across all read lengths.

The margin was belt-and-suspenders on top of the LOD model. In the
junction-homology case it was designed for: if both tumor and normal
credit borderline reads equally (r2c barely > native by 1-2 points),
the downstream LOD model sees similar split support in both samples →
low somlod → correctly not called somatic. Normal already used
margin=0, so the asymmetry was the only thing preventing normal from
crediting those reads — and it wasn't, because N_R2C_MIN_MARGIN was
always 0. The somatic/germline distinction is the LOD model's job.

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

Everything post-`svaba run` is now done by the **`svaba postprocess` C++
subcommand** (`src/svaba/SvabaPostprocess.cpp::runPostprocess`). It used to be
orchestrated by `scripts/svaba_postprocess.sh` (itself the replacement for the
old `sort_output.sh` + `sort_and_deduplicate_bps.sh` pair); the subcommand has
since absorbed every step that wrapper did -- per-thread merge, BAM sort/dedup,
the `bps.txt.gz` sort/dedup (formerly GNU `gsort`), and the r2c stamping -- so
the shell wrapper is now legacy/superseded. One command:

```bash
svaba postprocess -i ${ID} -t 8 -m 4G
```

Six phases per invocation, all idempotent (each auto-skips when its work is
already done, so reruns are near-instant):

0. **Merge per-thread outputs** (`mergeThreadBams` + `R2CDatabase::merge_from`)
   -- `${ID}.thread*.${suffix}.bam` -> `${ID}.${suffix}.bam` for `discordant`
   / `corrected` (open all, concatenate records, header from the first input;
   not yet coordinate-sorted -- Phase 1 does that), and per-thread
   `${ID}.thread${N}.r2c.db` -> `${ID}.r2c.db` via SQLite ATTACH + INSERT.
   Single-file inputs are renamed; no-file inputs skipped.
1. **Parallel sort** + **2. serial dedup / reheader / index** -- the BAM half
   of the job, structured exactly as before but now Phases 1-2 of the
   subcommand. For each suffix (`corrected`, `discordant`, `contigs`):
   - `samtools sort -@ per_job_threads -m MEM` (shell out — htslib doesn't
     expose its sort as a library call). **Auto-skipped** when the BAM
     already declares `@HD SO:coordinate` — `isCoordinateSorted()`
     inspects the header via `readHeaderOnly()` and logs
     "already coordinate-sorted; skipping sort" so reruns are a no-op.
   - Native streaming dedup (only for `corrected`/`discordant`):
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
       the PG line, so a second `svaba postprocess` run on the same
       outputs no-ops almost entirely.
     - The Phase 0 merge is itself a no-op when no per-thread `.thread*`
       files are present (nothing to merge), so all phases compose naturally.
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
3. **Sort + dedup + PASS/somatic filter of `bps.txt.gz`** -- fully in-process
   now (replaces the GNU `sort` + awk pipeline the shell script used).
   Decompress `bps.txt.gz` into one contiguous slab, build a lightweight
   index, and `std::sort` by `chr1(V), pos1, strand1, chr2(V), pos2, strand2,
   maxlod(desc)` (chr key matches GNU `sort -V`: chr1..22 -> 1..22, X -> 23,
   Y -> 24, M -> 25, non-standard -> end). **Dedup is by canonical breakpoint
   pair**: each pair is canonicalized (lesser breakend first) so reciprocal
   `A/B` and `B/A` rows collapse -- the old shell adjacency-sort dedup missed
   that. Per canonical pair the surviving row is the one with the **LOWEST
   somlod** (tie-break: higher maxlod) -- deliberately conservative against
   promoting a germline event to somatic. Emits FOUR files:
   `.bps.sorted.txt.gz` (all rows, sorted), `.bps.sorted.dedup.txt.gz` (one
   winner per canonical pair), `.bps.sorted.dedup.pass.txt.gz` (winners with
   `confidence == PASS`), and `.bps.sorted.dedup.pass.somatic.txt.gz` (PASS
   winners with the somatic flag set). `svaba tovcf` consumes
   `.bps.sorted.dedup.txt.gz`.
   - Column indices are hard-coded from `BreakPoint::header()` (0-based in
     the C++): `0 chr1 .. 5 strand2`, `29 cname`, `31 confidence`, `35 somatic
     flag`, `36 somlod`, `37 maxlod` (= 1-based cols 30/32/36/37/38). Change
     `toFileString`/`header` and these break.
4. **Stamp `pass_cnames` into `${ID}.r2c.db`** -- from the Phase-3 winner sets,
   `R2CDatabase::stamp_pass_cnames` writes a small `pass_cnames(cname,
   somatic)` table: every PASS winner's cname, with `somatic = 1` when the
   bps **somatic flag** (col 36, 1-based) is set -- NOT when `somlod >= 1`.
   Consumers then do `SELECT r.* FROM reads r JOIN pass_cnames p USING(cname)`
   (add `WHERE p.somatic = 1` for tumor-specific). Replaces the v3-era
   `.r2c.pass.txt.gz` / `.r2c.pass.somatic.txt.gz` copies -- same info, one
   lookup table.
5. **Optional split-by-source** (`--split`) -- demuxes the dedup-eligible
   BAMs by the first 4 chars of each QNAME into `${ID}.${suffix}.${prefix}.bam`.

CLI: `svaba postprocess -i <ID> [options]`. Flags: `-t/--threads`,
`-m/--mem` (per samtools-sort thread, e.g. `4G`), `--sort-only` (sort + @PG +
index, skip dedup), `--dedup-only` (skip sort, assume already coordinate-
sorted), `--split`, `-v/--verbose`, `-h/--help`. With the per-phase
auto-skips, a rerun on already-postprocessed files is effectively instant --
handy for refreshing just the index or the bps subsets.

The legacy `scripts/svaba_postprocess.sh` wrapper still exists and works, but
it predates the subcommand absorbing the merge / bps-sort / r2c steps and is
no longer the recommended path. Its old gotchas -- needing GNU `gsort` on
macOS for the bps step, and `gzip -dc` rather than BSD `zcat` -- don't apply
to `svaba postprocess`, whose bps sort is in-process.

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
  uses the same namespace as `r2c.db`'s `split_bps`/`disc_bps` and
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

## `svaba extract-pairs` subcommand

BAM-native pair extractor. Pulls every read pair from a BAM where either
mate's SEQ contains any of the given query sequences (or, by default,
their reverse complements). Lives in `src/svaba/SvabaExtractPairs.cpp`;
dispatch wired in `src/svaba/svaba.cpp`. The old
`scripts/extract_pairs_by_seq.sh` is gone — this subcommand fully
subsumes it (BAM-native pipeline, no SAM-text round-trip).

Two-pass design, both passes BAM-native:

- **Pass 1 — QNAME collection.** Stream input via `SeqLib::BamReader`
  with a BGZF decompression thread pool (`SetThreads`). For each
  record we walk `bam_get_seq(rec.raw())` directly through an
  Aho-Corasick automaton (5-letter alphabet `{A,C,G,T,N}`) using a
  16-entry nibble→alphabet lookup table — no `std::string` allocation
  per record, no SAM text, no regex engine. Patterns containing
  non-ACGTN bases are silently rejected at insertion (they would never
  match an htslib-stored sequence anyway). On any pattern hit the
  read's QNAME goes into a hash set. In `--counts` mode (see below) the
  AC search is `searchNibblesCollect` instead of the early-exit
  `searchNibbles`, so every pattern hit in the SEQ gets attributed back
  to the bp_id(s) that emitted that kmer.
- **Pass 2 — pair extraction.** Re-stream input. Every record whose
  QNAME is in the set gets written to the output BAM. Both mates plus
  any supplementary/secondary alignments survive together because they
  share the QNAME. BGZF write pool on the writer.
- **Sort.** Skipped iff input declares `@HD SO:coordinate` (the common
  case). Output preserves input record order, so coord-sorted in →
  coord-sorted out. For unsorted input we shell out to `samtools sort`
  — same fallback as `svaba postprocess`. The legacy shell script
  unconditionally sorted, which on a large already-sorted BAM is a
  wasted full read+write of the output.
- **Index.** `sam_index_build(out, 0)` directly — no `samtools index`
  fork.
- **@PG stamp.** A bare `@PG ID:svaba_extract_pairs PN:svaba VN:<ver>
  CL:<argv>` line is appended to the output header. Unlike
  `svaba postprocess`, we do NOT uniquify the ID — extract-pairs
  produces a fresh output BAM, so a duplicate ID only appears if the
  user piped extract output through extract again, where a duplicate
  trace is informative rather than wrong.

CLI:
```
svaba extract-pairs -i IN.bam -o OUT.bam (-s SEQ ... | -f FILE) [options]
  -s, --seq SEQ           Query sequence; repeatable.
  -f, --seq-file FILE     File of query sequences. Two formats accepted,
                          auto-detected by content: a plain one-seq-per-line
                          list (# / blank lines ignored), or a svaba
                          bps.txt[.gz] dump (the `jxn_kmer` column,
                          col 53, is used as the query; rows with kmer
                          == "." are skipped).
  -t, --threads N         BGZF reader+writer threads. [4]
      --no-rc             Skip reverse-complement augmentation.
      --no-pairs          Single-pass mode: emit only records whose own
                          SEQ matched. Skips the mate / supplementary
                          pickup that the default two-pass mode does;
                          ~2x faster (one BAM pass instead of two, no
                          QNAME hash set). Use when you just want to
                          inspect which reads contain a motif.
      --counts FILE       Emit a per-bp_id TSV of unique-non-dup reads
                          carrying each bp's kmer. Header is
                          `bp_id<TAB>n_unique_reads`. "Unique" = primary
                          alignments only (excludes flag 256 secondary,
                          2048 supplementary, 1024 duplicate), dedup'd
                          by (bp_id, qname, mate). Requires `-f` to be
                          a bps.txt[.gz] file (the only source of the
                          kmer↔bp_id map); incompatible with --no-pairs.
                          The counts file is written even when zero
                          reads matched — you get a header-only TSV in
                          that case.
  -v, --verbose 0-3
  -h, --help
```

Why this is faster than the old shell script (now deleted):
- No `samtools view → awk → samtools view -b` SAM-text round-trip.
  The legacy script was bottlenecked on a single-threaded awk consumer
  of decompressed SAM text; `samtools view -@ N`'s decompression
  threads couldn't help past the pipe.
- No regex alternation per read. Aho-Corasick gives O(read_len) per
  read regardless of pattern count.
- No third pass for `samtools sort` when input is already coord-sorted
  (skipped via `@HD SO:coordinate` header check).
- BGZF thread pools on both reader and writer.

Useful jump points:
- Entry point: `src/svaba/SvabaExtractPairs.cpp::runExtractPairs`
- Aho-Corasick (5-letter alphabet, nibble-direct search; carries
  `pattern_id` + `output_link` per node so `searchNibblesCollect` can
  enumerate every hit for the per-bp_id counting path):
  `class AhoCorasick` in same file.
- bps reader (kmer + bp_id extraction): `readJxnKmersFromBps` in same
  file. Returns a `LoadedSeqs` with both the kmer list and the
  `bp_ids_by_kmer` map.
- Per-bp_id tally: lives inside `collectQnames` when the
  `pattern_to_bp_ids` / `bp_id_counts` args are non-null. Dedup key
  is `bp_id + qname + mate`. TSV is written out of
  `runExtractPairs` after pass 1.
- Coord-sort detection helper: `isCoordinateSorted()` in same file
  (kept local rather than pulling in `SvabaPostprocess.h`).

## Homology / insertion strand convention (v4)

For SVs (ASSMB / ASDIS), `BreakPoint::set_homologies_insertions` reads
the homology and inserted-sequence bytes out of `seq` (= the contig
`m_seq` in assembly-native orientation) by slicing between `b1.cpos`
and `b2.cpos`. m_seq has no inherent strand — it's whatever bytes the
assembler emitted — so when BWA reports the left fragment with
`ReverseFlag=true`, those bytes are the reverse-complement of side 1's
forward-strand reference.

Convention: **homology and insertion are always stored on side 1's
forward strand.** When `b1.gr.strand == '-'` (which the
`isleft=true` branch of `BreakEnd::transferContigAlignmentData` sets to
match `ReverseFlag` directly), the slice is reverse-complemented before
being assigned to `homology` / `insertion`. The fix kicks in only when
the contig happens to come off the assembler in the reverse-of-side-1
orientation — for the common cases (deletions, tandem dups, contigs
that align forward to side 1) it's a no-op.

Implications:
- `bps.txt.gz` cols 27 (homol) and 28 (insert) are unambiguous
  forward-strand spellings. Users can grep the reference + strand for
  the value directly and find the insert.
- `extract-pairs -f bps.txt.gz` (which already searches both forward
  and reverse-complement of every query by default) is unaffected —
  the canonicalization just changes which orientation gets searched
  first; the result set is identical.
- `tovcf`'s INFO/HOMSEQ and INFO/INSSEQ inherit the canonical form
  through `BreakPoint`'s fields.
- Indels (the indel BreakPoint ctor at line ~931) read insertion bytes
  from `m_align->Sequence()`, which is BAM SEQ — already
  reference-forward by SAM convention. So the indel path needs no fix
  and gets none.

## Junction kmer (v4 schema)

Col 53 of `bps.txt.gz` is `jxn_kmer` — a 20-bp slice of the contig
sequence that spans the breakend junction. Lives in
`BreakPoint::junctionKmer()` (BreakPoint.cpp). Window definition: 10 bp
ending at `cpos_on_m_seq().first` + 10 bp starting at `+1`. So for SVs
without inserted novel sequence (the common case where `c2 == c1+1`),
the kmer is exactly the contig spelling that split-supporting reads
carry verbatim. For indels with an insertion, the right half lands
inside the inserted bases — a deliberate compromise to keep the kmer
contig-contiguous (matches reads that cross the upstream junction).

Same forward-strand canonicalization as homology / insertion: when
`b1.gr.strand == '-'` (frag_left reverse-aligned), the kmer is
reverse-complemented before storage so cols 27, 28, and 53 are all
in side 1's forward-strand spelling. For `extract-pairs` queries this
is functionally a no-op (the AC matcher already runs both
orientations), but it keeps the on-disk strings consistent and
greppable against the reference + strand.

Emitted as "." when no clean kmer is definable: imprecise BPs,
DSCRD-only clusters, missing/unset cpos, empty contig, or contig too
short to fit the window.

The kmer round-trips through refilter via a parsed-cache field
(`BreakPoint::jxn_kmer`) so re-emission preserves the original value
even when `seq` is empty post-parse.

Primary downstream consumer: `svaba extract-pairs -f bps.txt.gz`.
That command auto-detects bps input by sniffing the `#chr1\t` header
prefix and pulls the kmer column out of each row, skipping "." rows.
Both forward and reverse-complement searches are run by default
(`--no-rc` to disable).

Readers updated for the v4 schema:
- `BreakPoint::header()` and the bps parser (`BreakPoint.cpp`,
  same colon-test heuristic used for col 52 → col 53).
- `docs/bps_explorer.html` (new `V4` constant + version-detection
  branch in `ingest()`).
- `scripts/svaba_local_function.sh::svaba_bps_cols` reference text.
- `scripts/svaba_postprocess.sh` and `SvabaPostprocess.cpp` need
  no changes — sort/dedup keys (cols 30/32/37/38/52) are unaffected
  and the dedup logic preserves whole lines via byte slabs.

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
51 cols, LEGACY had 41). It's also carried into `r2c.db` (see
next section) so a read's support attribution is unambiguously
linked to the exact BP row in bps.txt, eliminating the old "which BP
on this contig did this read actually support?" puzzle.

`svaba_bps_cols` (from `scripts/svaba_local_function.sh`) documents
the full layout; column 52 is the bp_id field.

**BAM aux tags `bi:Z` / `bz:Z` (v3).** The two aux tags svaba stamps
on weird/discordant/corrected BAM outputs now live in *different*
identifier namespaces — choose the right one for the join you want:

- `bi:Z` — comma-joined list of **bp_ids** this read supports as
  ALT. Matches the per-BP resolution of `r2c.db`'s `split_bps`
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

## r2c SQLite database (v4)

`${ID}.r2c.db` is the queryable r2c alignment dump. It replaces the v3
TSV (`${ID}.r2c.txt.gz`), which itself replaced the v2 ASCII output —
same information content, written directly into SQLite instead of going
through a TSV intermediate. No build-string-then-gzip round trip on the
hot path; arbitrary SQL on the back end via `sqlite3` CLI / `sql.js` /
`DuckDB.attach`.

The legacy `r2c_explorer.html` still loads old `.r2c.txt.gz` files for
historical outputs; the new viewer is `docs/r2c_db_explorer.html`
(sql.js, runs queries client-side).

**Per-thread emission + postprocess merge.** Each svaba worker writes
its own `${ID}.thread${N}.r2c.db` during the run via
`svabaThreadUnit::r2c_db_` (a `std::unique_ptr<R2CDatabase>`; opened in
the ctor, closed in the dtor, gated on `opts.dump_alignments`).
`SvabaOutputWriter::writeUnit` calls `alc.writeToR2cDb(...)`
**before** `writeMutex_` is acquired — each thread has its own SQLite
handle, so inserts run in parallel across all workers with no cross-
thread locking. Postprocess merges the per-thread `.db` files via
`R2CDatabase::merge_from()` (ATTACH the source as `src`, then
`INSERT OR IGNORE INTO main.contigs SELECT * FROM src.contigs;
INSERT INTO main.reads SELECT * FROM src.reads;`, then DETACH). The
first input is renamed in-place to the target so most of the data is
reused without copy.

Schema (built by `R2CDatabase::R2CDatabase()` in
`src/svaba/R2CDatabase.cpp`):

```
contigs(cname TEXT PK, contig_len INT, seq TEXT, frags TEXT, bps TEXT, n_reads INT)
reads  (cname TEXT FK, read_id TEXT, chrom TEXT, pos INT, flag INT,
        r2c_cigar TEXT, r2c_start INT, r2c_rc INT, r2c_nm INT,
        support TEXT, split_bps TEXT, disc_bps TEXT,
        r2c_score REAL, native_score REAL, seq TEXT)
        -- index: idx_reads_cname ON reads(cname)
pass_cnames(cname TEXT PK, somatic INT)   -- stamped by postprocess Phase 4
```

`pass_cnames` is the SQL-friendly equivalent of the v3 TSV-era
`.r2c.pass.txt.gz` / `.r2c.pass.somatic.txt.gz` filtered subsets:
`SELECT r.* FROM reads r JOIN pass_cnames p USING(cname);` for the
PASS subset, add `WHERE p.somatic = 1` for tumor-specific. One small
lookup table replaces gigabytes of duplicated rows.

`split_bps` / `disc_bps` are comma-joined `bp_id` lists — the
unambiguous per-BP attribution of each read. The categorical `support`
field (`split` / `disc` / `both` / `none`) is derived from these and
kept for query-friendliness. The `bps` column on the contig row also
carries each BP's id as the 2nd subfield, so viewers can join back
without a second table.

`frags` and `bps` retain their nested-string encoding from the TSV:
- `frags`: `|`-separated per fragment; within a frag, `:`-separated
  `chr:pos:strand:cigar:mapq:cpos_break1:cpos_break2:gpos_break1:gpos_break2:flipped`.
- `bps`: `|`-separated; within a row, `:`-separated 10 fields:
  `kind:bp_id:chr1:pos1:strand1:chr2:pos2:strand2:span:insertion`.
  `kind ∈ {global, multi, indel}`. Insertion is `.` when absent.

These could be promoted to proper relational tables (`contig_frags`,
`contig_bps`) in a follow-up — viewers parse them client-side today
because the TSV used the same encoding and there's value in keeping
the parser portable.

**Sentinels.** Empty `r2c_cigar` / `seq` (rather than NULL) when no
r2c alignment is available; score columns are `0.0` in that case.
SQL queries can filter via `WHERE r2c_cigar != ''` or
`WHERE r2c_score > 0`.

Performance: per-thread inserts run inside a single open transaction
(committed in `svabaThreadUnit::clear`) with `journal_mode=WAL` and
`synchronous=OFF`. This is faster than gzip-text emission was — the
TSV had to format a string for every row; SQLite binds primitives
directly. On WGS data we measured ~1.7× speedup on the dump path
end-to-end (build + write).

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
drop every read at `SvabaBamWalker.cpp:218`. Now those regions never
hit a thread.

Safe because `sc.blacklist` has had `MergeOverlappingIntervals()` +
`CreateTreeMap()` called, so `FindOverlapWidth` can't double-count and
wrongly drop a partially-callable region.

The per-read blacklist check (`SvabaBamWalker.cpp:218`) and the per-BP
check (`SvabaRegionProcessor.cpp:1262`, `checkBlacklist`) still run for
regions that **partially** overlap — this prune only short-circuits the
100%-covered case.

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
    - per-thread `${ID}.thread${N}.r2c.db` (merged into `${ID}.r2c.db` by
      the postprocess step; see "r2c SQLite database")
  The fields are kept separate so individual callsites can key off
  their own concern (e.g. `svabaThreadUnit` gates `r2c_db_` opening
  on `dump_alignments` only), but there's no runtime path to toggle
  them individually — all three flip as a unit under `--dump-reads`.
  (Heads-up: the `--help` text and a comment in `SvabaOptions.{h,cpp}`
  still say `r2c.txt.gz` — stale strings; the actual emission is the
  SQLite `.r2c.db`.)
- Without `--dump-reads`, svaba produces the lean output set only:
  `bps.txt.gz`, VCFs, `contigs.bam`, `runtime.txt`, `discordant.txt.gz`
  (cluster-level, tiny). No per-thread `r2c.db`, no `corrected.bam`,
  no `discordant.bam`. This is the production default because the gated
  outputs can run to tens of gigabytes on deep samples.
- **`alignments.txt.gz` is gone.** The pre-rendered ASCII viewer output
  was replaced by the structured r2c dump (first `r2c.txt.gz`, now the
  SQLite `r2c.db` — same information, not pre-formatted).
  `AlignedContig::printToAlignmentsFile` and
  `BreakPoint::printDeletionMarksForAlignmentsFile` were removed; the
  surviving `AlignmentFragment::printToAlignmentsFile` is kept only
  because one `std::cerr` debug print in `BreakPoint.cpp` still calls
  it. `docs/alignments_viewer.html` still exists and still works on
  old `.alignments.txt.gz` files from previous runs, but new runs don't
  produce that file — use `docs/r2c_db_explorer.html` instead.

## Viewer suite (`docs/`)

Entry point is `docs/index.html`, a card grid pointing at the
sub-viewers (the suite moved from `viewer/` to `docs/`). All
client-side, no server required.

- **`bps_explorer.html`** — primary viewer. Sortable table of bps rows,
  numeric filters (somlod/maxlod/qual/span/etc.), chip filters (counts
  are live — they reflect the current filter, not the full dataset),
  per-sample detail panel, small histograms for somlod/maxlod/span,
  IGV-navigation links in the IGV1/IGV2 columns (fires
  `fetch('http://localhost:60151/goto?locus=…', {mode:'no-cors'})`;
  requires IGV running with port 60151 enabled). r2c re-plot sub-panel
  was removed — that capability now lives in the standalone
  `r2c_explorer.html` below.
- **`r2c_db_explorer.html`** — current r2c viewer. Loads a `${ID}.r2c.db`
  (sql.js, runs SQL client-side) and renders the alignment plots for
  v4 SQLite-path runs. Reach for this one first.
- **`r2c_explorer.html`** — legacy re-plotter for the structured
  r2c TSV (`.r2c.txt.gz` from older runs, or filtered to
  PASS / PASS-somatic by `svaba postprocess`). Upload an
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
- **`learn_explorer.html`** — explorer for svaba's insert-size learning
  output (per-read-group insert-size / read-length distributions, learned
  inside `svaba run`; pairs with `scripts/plot_learn.sh`).
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
`bps.txt.gz` files and runs `svaba postprocess`.

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
dedup step in `svaba postprocess` pairs them. No calls are lost.

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
  the artifact model. These names are used consistently throughout this
  file.
- **Option codes** in `SvabaOptions.cpp::longOpts`: 1001-1099 = mode,
  1100s = assembly, 1200s = EC, 1300s = discordant, 1400s = filter,
  1500 = chunking, 1600s = bwa-mem tuning, 1700s = output/DBs, 1800 =
  dump-reads. Keep the ranges coherent when adding new options.

## Mate-region lookup pipeline

When svaba encounters discordant reads (insert size too large or wrong
orientation), it collects their mate loci and considers doing a secondary
"mate-region" assembly to catch the other breakend of an SV.

**Six gates** a mate candidate must pass (in order):

1. **Primary MAPQ** (`minMateMAPQ`, default -1 = no gate): the discordant
   read itself must have MAPQ ≥ this. Set to e.g. 10 to skip
   multi-mapped primaries.
2. **Chromosome ID** (`maxMateChrID`, default 23): mate must land on
   chr ≤ this ID (0-indexed: 0=chr1 .. 22=chrX, 23=chrY). Skips
   chrM/alt/decoy in human. Set to -1 (via `--non-human`) to disable
   entirely for non-human genomes.
3. **Blacklist**: mate locus checked against `sc.blacklist`.
4. **Min count** (`mateRegionMinCount`, default 2): merged region must
   have ≥ N supporting reads to survive the BamWalker filter.
5. **Somatic mateLookupMin** (default 3, `MATE_LOOKUP_MIN`): in
   `SvabaRegionProcessor`, only look up regions with ≥ this many
   somatic-only reads.
6. **Max regions** (6): cap at 6 mate regions per assembly window.

All constants live in `SvabaOptions.h` as `inline constexpr` with
runtime overrides in the `SvabaOptions` class:

```
--min-mate-mapq N     (default -1, no gate)
--max-mate-chr  N     (default 23, through chrY; set -1 for no limit)
--mate-min-count N    (default 2)
--non-human           (sets maxMateChrID = -1, removes human assumptions)
```

Code: `SvabaBamWalker.cpp::calculateMateRegions()`.

## Compile-time read & contig tracing

Two zero-cost compile-time trace systems for debugging why a specific
read was or wasn't credited / a contig was or wasn't called:

**`SVABA_TRACE_READ`** — traces a single read (by QNAME) through the
entire pipeline: BamWalker intake → BFC correction → r2c alignment →
native realignment → splitCoverage scoring → output tagging.

```bash
cmake .. -DCMAKE_CXX_FLAGS='-DSVABA_TRACE_READ="\"LH00306:129:227V5CLT4:6:1204:38807:7191\""'
```

Trace points (19 total across 3 files):
- `SvabaBamWalker.cpp`: initial read filter decisions (existing)
- `SvabaRegionProcessor.cpp`: BFC correction result, r2c alignment
  per-contig, native realignment reuse/done/miss, bi:Z tagging
- `BreakPoint.cpp`: splitCoverage entry, TP8 r2c-vs-native comparison,
  TP9 del/ins near break, TP10 span check, TP11 del covers,
  CREDITED/NOT CREDITED final decision

**`SVABA_TRACE_CONTIG`** — traces a single contig through assembly,
alignment, and scoring:

```bash
cmake .. -DCMAKE_CXX_FLAGS='-DSVABA_TRACE_CONTIG="\"c_fermi_chr2_215869501_215894501_13C\""'
```

Both can be combined. Both are `#ifdef`-guarded so they compile to
nothing when not defined. See `src/svaba/SvabaDebug.h` for the macro
definitions and `README.md` for full recipes.

## Useful jump points

- Somatic LOD calc: `src/svaba/SvabaModels.cpp:86`
- Per-sample LO: `src/svaba/BreakPoint.cpp:2002`
- Somatic LOD entry: `src/svaba/BreakPoint.cpp:1341`
- INDEL somatic gate: `src/svaba/BreakPoint.cpp:1434`
- Region-queue blacklist prune: `src/svaba/run_svaba.cpp` (right after
  `loader.countJobs(regionsToRun)`)
- Per-read blacklist filter: `src/svaba/SvabaBamWalker.cpp:218`
- r2c SQLite writer: `src/svaba/AlignedContig.cpp::writeToR2cDb` +
  `src/svaba/R2CDatabase.cpp` (old `printToR2CTsv`/`r2cTsvHeader` TSV
  emitter is gone)
- Postprocess driver: `src/svaba/SvabaPostprocess.cpp::runPostprocess`
- Postprocess (legacy shell wrapper): `scripts/svaba_postprocess.sh`
- `svaba tovcf` driver: `src/svaba/tovcf.cpp::runToVCF`
- VCF engine (parse + dedup + emit): `src/svaba/vcf.cpp` + `vcf.h`
- Symbolic SV classifier: `vcf.cpp::classify_symbolic_kind`
- Single-file VCF writers: `vcf.cpp::writeSvsSingleFile` + `writeIndelsSingleFile`
- Blacklist combiner: `scripts/combine_blacklists.sh`
- Runtime-file schema: `src/svaba/SvabaUtils.cpp::svabaTimer::header`
- Options parsing: `src/svaba/SvabaOptions.cpp::SvabaOptions::parse`
- Somlod/maxlod analysis: the "Statistical model" + "somlod / maxlod
  investigation" sections of this file (the standalone HTML writeup is gone)
- Mate-region lookup: `src/svaba/SvabaBamWalker.cpp::calculateMateRegions`
- Mate-region constants: `src/svaba/SvabaOptions.h` (lines 113, 122-143)
- Read trace macro: `src/svaba/SvabaDebug.h`
- Read trace (BFC/r2c/native): `src/svaba/SvabaRegionProcessor.cpp`
- Read trace (splitCoverage): `src/svaba/BreakPoint.cpp`
- Debugging recipes: `README.md` (in svaba_opt root)
