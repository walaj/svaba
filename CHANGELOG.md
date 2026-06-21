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
