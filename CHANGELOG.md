# Changelog

All notable changes to svaba are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and svaba follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
(`MAJOR.MINOR.PATCH`):

- **MAJOR** — incompatible changes to output schemas (e.g. `bps.txt.gz`
  columns), the VCF representation, or CLI behavior that breaks existing
  pipelines.
- **MINOR** — new functionality (subcommands, flags, outputs) that is
  backward compatible.
- **PATCH** — bug fixes and internal changes with no interface impact.

## How versioning works in this repo

- The human-edited version and date are the two constants
  `SVABA_VERSION` / `SVABA_DATE` in
  [`src/svaba/SvabaOptions.h`](src/svaba/SvabaOptions.h). To cut a release,
  bump those, add a section here, and tag the commit (`git tag vX.Y.Z`).
- The **build commit** is stamped automatically by CMake at build time into
  a generated header (`SvabaGitVersion.h`) — no manual step. `svaba
  --version` prints the semantic version plus that commit (and `describe`,
  and a `dirty` flag if the working tree had uncommitted changes), so a
  binary's exact provenance is recoverable even from a Docker image whose
  `.git` was removed. The commit also appears in the run startup banner and
  in the `@PG` lines of output BAMs.
- If svaba is built from a non-git source tree, the commit shows as
  `unknown` and `--version` prints only the semantic version.

## [Unreleased]

_Nothing yet._

## [2.5.0] - 2026-08-18

### Added
- **`svaba extract-pairs` accepts a multi-record FASTA** as `-f`. The file
  format is auto-detected from the first non-empty line, so `-f` now takes a
  FASTA, a plain one-sequence-per-line list, or a `bps.txt[.gz]` dump (plain
  or gzip in every case). A record's name is the first whitespace-delimited
  token of its `>` header; wrapped sequence lines are concatenated; duplicate
  names are uniquified (`name`, `name.2`, …) so each record keeps its own
  row; records with an empty sequence are skipped with a warning.
- **`--counts` gained an `n_unique_qnames` column** — the number of distinct
  read pairs (QNAMEs) carrying each query, alongside the existing
  `n_unique_reads` (distinct mates). A pair whose two mates both match counts
  once in the new column and twice in the old one.

### Changed
- **`--counts` works for every query source, not just `bps.txt[.gz]`.** The
  row id is the `bp_id` (bps), the FASTA record name, or the sequence itself
  (`-s` / plain list); the header labels follow (`bp_id`/`jxn_kmer` vs
  `seq_name`/`sequence`). For FASTA / plain / `-s` queries every query gets a
  row, including ones that matched nothing (zero counts are a real answer for
  a small probe set); `bps` input still lists only ids with at least one
  match. The previous "`--counts` requires `-f` to be a bps file" error is
  gone; `--counts` remains incompatible with `--no-pairs`.
- **`--counts` writes CSV when the output filename ends in `.csv`** (RFC-4180
  quoted), TSV otherwise. Existing TSV output is unchanged apart from the
  appended fifth column.

## [2.4.0] - 2026-06-25

### Added
- **Recurrent 5′ soft-clip "tag" detection + auto-trim.** An untrimmed 5′
  molecular tag (UMI / inline barcode / adapter / library spacer) mismatches the
  reference at the read's 5′ end on nearly every read. That mismatch lands at the
  leading edge of each read's overlap with its neighbor, which breaks the overlap
  (string-graph) assembly genome-wide — so svaba can silently lose calls on an
  otherwise-fine BAM (measured: 133→94 calls on a 2 Mb region with a synthetic
  5 bp tag).
  - **Detect during learning:** `BamReadGroup::addRead` accumulates a
    5′-soft-clip-length histogram (strand-aware) and `detectTag()` flags a sharp
    short-clip spike (≥20% of reads at the modal length, ±1 bp) as a tag,
    classifying *fixed-sequence* (adapter/spacer) vs *variable* (UMI).
  - **Warn:** a prominent per-read-group `WARNING` naming the length, fraction,
    and fixed sequence (or "likely UMI"), with a pointer to trim upstream.
  - **Record + propagate:** `tag_trim_5p` is written to the `.learn.tsv.gz` and
    to the `--bam-params` cache, so scatter-gather shards reuse the detected
    value with no re-detection.
  - **Auto-trim:** svaba trims `min(5′-softclip, L)` bases off each read's 5′ end
    before assembly (`svabaRead::TrimTag5p`). Self-targeting (only clipped,
    non-genomic bases up to `L`; a long split-read clip keeps everything past
    `L`) and assembly-only (the original BAM record/CIGAR is untouched, so r2c,
    scoring, and coordinates are unaffected).
  - **Flags:** `--tag-trim N` forces N bp; `--no-tag-trim` disables the action
    (detection + warning still run).
  - Verified: auto-trim fully recovers the broken assembly (94→133, exact match
    to the untagged control); inert and byte-identical on clean BAMs (no false
    positives).

## [2.3.3] - 2026-06-25

### Fixed
- **`--max-cov` no longer changes the error model of clean reads in the same
  window.** The weird-coverage subsample demotes pileup reads out of fermi
  assembly to cut the (super-linear) assembly cost, but it had *also* removed
  them from BFC k-mer **training**. Because BFC auto-learns `k` and its k-mer
  spectrum from the pooled window reads, capping a pileup changed how *every*
  read in that ~25 kb assembly window was error-corrected — silently re-calling
  (and sometimes dropping) clean variants nowhere near the pileup.
  - Diagnosed on a clean somatic Alu–Alu deletion at `chr1:44147919/44149082`
    (somlod 2.14): instrumentation showed **0 of its breakend reads were
    demoted** — every demoted read was in a pileup ~9 kb upstream — yet the call
    vanished at `--max-cov 100`, and only with error-correction on (`--ec-type 0`
    never produced it). The coupling was the window-global BFC spectrum, not the
    coverage measure (which is already per-base).
  - **Fix:** demoted reads now carry `svabaRead::cov_demoted` and remain in the
    BFC **training** pool (Phase 1) while staying excluded from the correction
    pool and from fermi. Training is ~linear and cheap, so the k/spectrum match
    the no-cap run → reads elsewhere correct identically → the call returns
    **byte-identical** (same somlod/maxlod/FILTER; only the contig serial
    differs). Measured chr1 runtime give-back: ~+27 s CPU (~1.5%), i.e. the
    `--max-cov` speedup is preserved.
  - Remaining cap-vs-no-cap differences are now only the unavoidable fermi-graph
    effect of not assembling a pileup, confined to the marginal (somlod ≤ 3)
    tier. An env-gated probe `SVABA_DEBUG_SUBSAMPLE=1` reports per-region
    demotion counts and the max weird-coverage seen.

## [2.3.2] - 2026-06-25

### Fixed
- **Byte-identical output across runs (determinism, part 2).** After 2.3.1 made
  the *calls* deterministic, two remaining sources of run-to-run variation
  remained, both cosmetic but enough to make a plain `diff` of the output look
  different:
  - **`bp_id` (bps col 52) was thread-derived.** The old format
    `bpTTTNNNNNNNN` encoded the worker-thread ID + a per-thread counter, so the
    same variant got a different id depending on which thread (and in what
    order) processed its region. The `bi:Z` BAM tag, `r2c.db`
    `split_bps`/`disc_bps`, and the VCF `EVENT` field all key on it, so those
    joins weren't reproducible across runs either. It is now derived from the
    BP's stable content (contig name + both breakends + their contig positions +
    SV kind) via FNV-1a → `bp` + 16 hex digits. The same BP yields the same id
    across runs, threads, and `-p` settings; the tuple uniquely identifies a BP,
    so the hash can never merge two distinct variants onto one id.
  - **`svaba postprocess` bps sort could reorder equal-key rows.** The sort
    (`std::sort`, unstable) keyed on coords + maxlod, so rows with identical
    keys were left in worker-thread emission order. Added a final total-order
    tie-break on the raw line bytes (now deterministic), so the sorted output is
    stable.
  - **Result:** a plain `diff` of `.bps.sorted.txt.gz` (and the
    `.dedup` / `.pass` / `.pass.somatic` subsets) is now `0` across runs.

### Changed
- **`bp_id` format** changed from `bpTTTNNNNNNNN` (thread+counter) to `bp` + 16
  hex digits (content hash). Tools that treat bp_id opaquely (all of svaba's
  own — `extract-pairs`, the viewers, the r2c/VCF joins) are unaffected;
  anything that parsed the old thread+counter layout must adapt.

## [2.3.1] - 2026-06-24

### Fixed
- **Nondeterminism under multithreading (`-p > 1`).** Identical inputs and
  command produced different calls run-to-run: ~0.3% of all breakpoints, and
  ~8% of the somatic-PASS subset (somatic classification is a threshold on the
  sensitive `somlod`, which amplifies any upstream wobble). Confirmed by
  bisection: deterministic at `-p 1`, nondeterministic at `-p 8`, and
  deterministic again with error-correction disabled (`--ec-type 0`).
  - **Root cause:** the BFC error-corrector auto-learns its k-mer size from the
    reads, but that value was **cached in the per-thread BFC object, which is
    reused across regions** (an alloc-reuse perf optimization). So `k` was
    frozen from whichever region a worker thread processed *first*. Because
    region→thread assignment is nondeterministic with multiple threads, the EC
    k-mer — and therefore error-correction, assembly, and the entire call set —
    varied between runs.
  - **Fix:** re-learn `k` from each region's own reads on every `BFC::Train()`
    (SeqLib), rather than caching it; an explicit `SetKmer()` still pins it
    (`kmer_fixed`). The count-hash allocation is still reused when consecutive
    regions resolve to the same `k`, so the perf optimization is preserved.
    **Verified:** `-p 8` now produces byte-identical calls to `-p 1`, and
    repeated `-p 8` runs are identical.
  - **Behavior change:** results differ slightly from 2.3.0 for regions whose
    own learned `k` differed from the frozen one — this is the corrected,
    reproducible output. (Not an RNG; there was no missing seed.)

## [2.3.0] - 2026-06-23

### Added
- **`--max-normal-weird-cov N` (opt-in, 0 = off) — somatic-safe artifact skip.**
  In any sub-region where the **normal** weird-read coverage exceeds `N`, the
  reads are dropped from assembly + error-correction + r2c + scoring entirely
  (DP coverage and discordant clustering already happened, so both are
  unaffected). Rationale: a pathological *normal* weird-read pileup is a shared
  artifact/repeat, and a somatic event requires a **clean** normal — so skipping
  these regions can never drop a somatic call. It is keyed on the **normal** only
  (a tumor-only pileup, which could be real somatic, is never skipped), which is
  what makes it somatic-safe by construction. Targets the few high-pileup regions
  that dominate runtime even after `--max-cov`: measured ~40%+ region-time
  reduction once it fires. The threshold must sit above germline weird-coverage
  (which scales with normal depth); ~200 is reasonable for deep WGS normals.
  Validate recall/precision on the sim panel at your depth before relying on it.
  Implemented cross-walker in `SvabaRegionProcessor::process` (after discordant
  clustering, before BFC), querying each normal walker's `weird_cov`.

## [2.2.1] - 2026-06-23

### Fixed
- **`--max-cov` was dead code; now enforced.** The per-position weird-coverage
  subsampler (`subSampleToWeirdCoverage`) had its only call site commented out in
  `svabaBamWalker::readBam`, and `walker->max_cov` was never wired to
  `opts.maxCov` — so the cap was parsed, stored, and even logged ("Max cov to
  assemble: N") but never applied. High-weird-coverage pileups therefore flooded
  the BFC + fermi-assembly + r2c pipeline (the dominant cost on deep regions),
  which is why high-depth windows had outsized runtime despite `--max-cov 100`.
  Fix: wire `walker->max_cov = opts.maxCov` in `SvabaRegionProcessor`, restore
  the call guarded by `if (get_coverage && max_cov > 0)`.
  - **Demote, don't delete.** Where per-position **weird-read** coverage exceeds
    the cap, the excess reads are not discarded — they're **demoted to
    second-class** (`to_assemble=false`): excluded from fermi assembly + BFC
    error-correction (the read-count-scaling costs, ~43% of compute on deep
    regions = the speedup) but **kept** in the read set so they are still
    r2c-aligned, discordant-clustered, and scored. This preserves support for
    very-high-coverage real events (amplicons / double-minutes) and — critical
    for somatic safety — never samples away a normal read that could refute a
    somatic call. Deterministic hash-by-qname (mate pairs demoted together);
    symmetric across tumor/normal. `--max-cov 0` disables it.
  - **Behavior change:** results differ from 2.2.0 only in regions whose
    weird-read coverage exceeds the cap; validate recall/precision on the sim
    panel before relying on it for calling.

## [2.2.0] - 2026-06-22

### Added
- **`--bam-params FILE` / `--learn-only` — cache insert-size learning across a
  scatter-gather run.** `learnParams()` samples up to ~1 M reads/window across
  every chromosome midpoint (≈ many millions of reads, ×2 for tumor+normal) to
  estimate per-read-group insert size — and it ran on **every** `svaba run`
  regardless of `-k`. In a 322-shard scatter that genome-wide sweep (+ 5.4 GB
  index load) was paid 322×, ≈ 5 CPU-min/shard of pure redundant overhead (≈ a
  quarter of the total CPU-h) plus billions of random reads against the shared
  BAM disk. Now: `svaba run --learn-only --bam-params p.tsv` learns once and
  writes a tiny TSV; every shard runs `--bam-params p.tsv` to LOAD it and skip
  the sweep entirely. The cache stores per-RG `isize_median`/`sd_isize` + per-BAM
  `readlen/mapq/isize` maxima, matched to the run's BAMs by path; if the file is
  missing/incomplete svaba transparently learns (and writes it). No effect on
  results — only on how insert size is obtained.

## [2.1.2] - 2026-06-22

### Changed
- **Contig-alignment confidence now penalizes alignment non-uniqueness.**
  `scoreContigAlignment` gained Rule E: when BWA's suboptimal score `XS`
  approaches the primary `AS` (`XS/AS > 0.80`, scaling to a full penalty at a
  dead tie), the contig has another near-equal placement — a repeat / paralog /
  segmental-dup locus that BLATs to many sites — so `contig_conf` is lowered and
  the existing `WEAKCONTIG` gate catches it. The `as_xs_gap` signal was already
  computed but unused; the prior confidence measured only absolute quality
  (mismatch rate, AS-density), not multi-mapping. Targets the dominant remaining
  false-positive mode: weak/ambiguous contig alignments svaba treated as
  confident. Caller-side → needs a fresh run. The per-breakend `contig_conf`
  (bps cols 40/41) is unchanged in meaning, just better calibrated; `benchmark.py
  --events-out` now also records it (min over breakends) so it can be tuned as a
  slider in `docs/sim_commit_compare.html`.

## [2.1.1] - 2026-06-22

### Fixed
- **Span-0 FR discordant false somatics.** `score_dscrd`'s `LOWSPANDSCRD` gate
  guarded on `getSpan() > 0`, but an intrachromosomal same-position discordant
  cluster has `getSpan() == 0` — a degenerate zero-length "deletion" formed when
  mildly-over-insert FR pairs (tagged discordant by a low insert cutoff) collapse
  to a single base. These leaked straight to `PASS` as false somatic DSCRD calls
  (e.g. 312 of 319 somatic-PASS DSCRD at 10× in the sim panel). Changed the guard
  to `getSpan() >= 0`, so span-0 FR clusters are correctly labeled `LOWSPANDSCRD`;
  interchromosomal events (`getSpan()` = −1) remain excluded. Caller-side → needs
  a fresh run.

## [2.1.0] - 2026-06-21

Backward-compatible feature additions plus a batch of false-somatic / coordinate
correctness fixes accumulated since 2.0.0. The `bps.txt.gz` and VCF additions are
trailing/additive (older parsers that read by name or early columns are
unaffected), so this is a MINOR release. See [`CLAUDE.md`](CLAUDE.md) for the
full engineering notes. Versioning is now Claude-maintained (see CLAUDE.md
"Versioning").

### Added
- **`--annotation <BED>`** (repeatable): overlap each breakend with a labeled
  BED (RepeatMasker / segmental-dups / custom). New trailing bps columns
  **`repeat_anno`** (col 55, `b1labels|b2labels`) and **`poly_a`** (col 56,
  longest poly-A/T run in the insertion = MEI signal); VCF `INFO/REPEAT_ANNO`
  and `INFO/POLYA`. Non-filtering. (bps schema v5 → **v6**.)
- **Somatic-safety junction-kmer net** → new per-sample FORMAT subfield
  **`KC`**. Reads excluded from assembly/r2c (adapter read-through,
  blacklist-self) that carry a breakpoint's junction kmer (≥19/20 bp, either
  strand) are counted and folded into the normal alt count for the somatic LOD,
  so a junction-spanning normal read can never be silently lost.
- **`sim/`** — self-contained SV simulator: non-overlapping SVs across the
  genome → diploid donor → ART → BWA-MEM, in-silico tumor impurity admixture,
  ground truth (VCF/BEDPE/TSV), `benchmark.py`, and one-shot panel/run drivers
  (`build_sim_panel.sh`, `run_svaba_on_panel.sh`).
- **Viewers**: `docs/sim_benchmark.html` (ROC / PR vs truth, dual somlod+maxlod
  gate, IGV-linked TP/FP/FN tables, blacklist masking) and `docs/truth_store.js`
  (in-browser Valid/Artifact/Germline/Ambiguous labeling with auto-save to a
  committed `ground_truth.json`), wired into `bps_explorer.html` /
  `comparison.html`. New `tracks/hg38.rmsk.bed.gz`, `tracks/hg38.segdups.bed.gz`.

### Changed
- **DUPREADS** now keys on the number of UNIQUE split-read genomic start
  positions (`nsplit_starts`), not the contig-footprint span — fixes
  false-flagging of real low-coverage / short-contig / duplicate-free events.
- Distant-read recovery for assembly-only large SVs; blacklist-aware
  region-queue pruning; coverage (DP) now counts blacklist-overlapping reads;
  `addCovs` averages only populated breakends.

### Fixed
- **`hasAdapter`**: a 5′-end soft-clip is split-read signal, never adapter
  read-through — adapter geometry is judged on the 3′ clip only, so
  junction-spanning reads across small SVs are no longer dropped.
- **SV breakend off-by-one**: `gbreak` stored 0-based; every SV breakend
  coordinate shifts −1 vs prior output (indels unaffected). Symbolic-SV REF base
  in `tovcf` now taken at the min breakend. VCF POS/END now consistently 1-based.
- **Insert-size robust SD (MAD)**: tames heavy-tailed isize → fixes inflated
  discordant cutoffs that blinded the normal and produced false somatic calls;
  interchromosomal RF discordant reads no longer dropped.
- Confidence read-count vetoes removed for SV paths (evidence sufficiency is the
  LOD's job); `LOWICSUPPORT` intra-chrom guard restored.

## [2.0.0] - 2026-06-07

The **SvABA 2.0** overhaul. This is a major release: output schemas, the
VCF representation, and several internal pipelines changed relative to the
1.x line. See [`CLAUDE.md`](CLAUDE.md) for the exhaustive engineering notes;
the highlights:

### Added
- **`svaba --version`** — prints the semantic version, the git commit the
  binary was built from, `git describe`, a dirty flag, and the build date.
- **`svaba postprocess`** subcommand — in-process merge, BAM sort/dedup
  (with BGZF thread pools), `bps.txt.gz` sort/dedup, and r2c stamping;
  replaces the old `sort_output.sh` / `sort_and_deduplicate_bps.sh`
  shell pipeline.
- **`svaba tovcf`** subcommand — standalone converter from a deduped
  `bps.txt.gz` to VCFv4.5 (`*.sv.vcf.gz` + `*.indel.vcf.gz`), with
  symbolic `<DEL>/<DUP>/<INV>` alleles where orientation is unambiguous
  and paired BND records otherwise.
- **`svaba extract-pairs`** subcommand — BAM-native read-pair extraction
  by sequence match (Aho-Corasick, both strands), replacing
  `extract_pairs_by_seq.sh`.
- **Per-breakpoint stable IDs** (`bp_id`, `bps.txt.gz` col 52) threaded
  through the BAM `bi:Z` tag and the r2c database for end-to-end read↔
  variant traceability.
- **Junction k-mer** (`jxn_kmer`, `bps.txt.gz` col 53) — a 20 bp
  contig-native sequence spanning each breakend junction.
- **r2c SQLite database** (`${ID}.r2c.db`) — queryable per-read
  read-to-contig alignment dump, written per-thread and merged in
  postprocess; replaces the pre-rendered `alignments.txt.gz`.
- **Blacklist-aware region pruning** — fully blacklisted regions are
  dropped before reaching a worker thread.
- **jemalloc support** via `-DUSE_JEMALLOC=ON` (recommended on Linux at
  `-p 16+`), plus the `svaba_jemalloc` LD_PRELOAD wrapper.
- **HTML viewer suite** under `docs/` (bps / r2c / runtime / learn
  explorers and a two-run comparison view).

### Changed
- **Somatic statistical model** reworked (`SvabaModels.cpp`): split-error
  somatic LOD with additional germline sub-hypotheses (`GERM_shared`
  shaping terms and a new independent-MLE `GERM_free` branch with a BIC
  penalty).
- **Split-coverage gate** is now a direct comparison — a read supports a
  breakend iff its r2c alignment scores strictly higher than its native
  alignment and it spans the breakend. The old homology-conditioned
  `both_split`/`one_split` branching and the `T_R2C_MIN_MARGIN` percentage
  margin were removed (the margin was read-length-dependent and silently
  dropped tumor support for small indels).
- **fermi-lite** is the default local-assembly engine; its k-mer hash is
  pooled per worker thread and reused across regions.
- **VCF output** declares `VCFv4.5`; somatic and germline records share one
  SV file and one indel file, distinguished by the `SOMATIC` INFO flag.
- **Homology / inserted sequence / junction k-mer** are canonicalized to
  side-1 forward-strand spelling in `bps.txt.gz`.
- Under jemalloc, the per-flush `malloc_trim(0)` calls are now compiled out
  (gated on `SVABA_USE_JEMALLOC`): `malloc_trim` is a glibc-only call that
  can't reclaim jemalloc-held memory and only exercised glibc's arena
  machinery. They remain active for glibc builds, where they bound RSS on
  high-coverage samples.

### Removed
- `alignments.txt.gz` and its ASCII emitter (superseded by the r2c
  database).
- The standalone post-processing and pair-extraction shell scripts
  (`sort_output.sh`, `sort_and_deduplicate_bps.sh`,
  `extract_pairs_by_seq.sh`), now built-in subcommands.

[Unreleased]: https://github.com/walaj/svaba/compare/v2.0.0...HEAD
[2.0.0]: https://github.com/walaj/svaba/releases/tag/v2.0.0
