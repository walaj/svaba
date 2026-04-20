## *SvABA* — Structural variation and indel analysis by assembly

SvABA (formerly *Snowman*) is an SV and indel caller for short-read BAMs.
It performs genome-wide local assembly with a vendored SGA, realigns
contigs with BWA-MEM, and scores variants by reassembled read support.
Tumor/normal, trios, and single-sample modes are supported; variants
are emitted as VCF plus a verbose tab-delimited companion
(`bps.txt.gz`) that carries the full per-sample evidence.

**License:** [GNU GPLv3](./LICENSE). Uses the [SeqLib][seqlib] API for
BAM I/O, BWA-MEM alignment, interval trees, and the assembly front-end.

## Install

CMake is required; htslib is an external dependency (no bundled
copy). If htslib is installed system-wide, svaba auto-detects it
via pkg-config or the default include/lib search paths:

```bash
git clone --recursive https://github.com/walaj/svaba
cd svaba && mkdir build && cd build
cmake .. && make -j
```

For a non-standard htslib install location, point at it explicitly:

```bash
cmake .. -DHTSLIB_DIR=/path/to/htslib-1.xx
make -j
```

The binary lands at `build/svaba`. To install it system-wide (or
into this repo's `bin/`), use CMake's install target — the default
prefix is `/usr/local` and needs sudo, or override it with
`-DCMAKE_INSTALL_PREFIX`:

```bash
make install                                           # /usr/local/bin/svaba
cmake --install build --prefix $(pwd)/..               # repo-local bin/svaba
```

Build type defaults to `RelWithDebInfo` (`-O2 -g -DNDEBUG`). The
vendored `SeqLib/bwa` and `SeqLib/fermi-lite` hardcode `-O2` in their
Makefiles; see `CLAUDE.md` for how to push them to `-O3 -mcpu=native`
(typically a 5–15% wall-time win).

### Assembler selection

By default svaba assembles contigs with **fermi-lite** (ropebwt2-based,
the faster path). The vendored SGA (String Graph Assembler) is still
fully supported but off by default — switch it on at compile time by
passing `-DSVABA_ASSEMBLER_FERMI=0` to CMake:

```bash
cmake .. -DSVABA_ASSEMBLER_FERMI=0
make -j
```

The choice is a single compile-time symbol in
`src/svaba/SvabaAssemblerConfig.h`; you can also flip it by editing that
file directly and rebuilding. Runtime logs report the active
assembler under `svaba::kAssemblerName`, so you can always confirm
which engine a build was compiled with.

## Quick start

Three steps. Run the caller, merge + dedup + index the per-thread
outputs, then convert the deduped `bps.txt.gz` to VCF. The bundled
combined blacklist is strongly recommended — it keeps svaba out of
well-known pileup / high-complexity regions that would otherwise
dominate wall-clock time with no real calls to show for it. See
`tracks/README.md` for customizing or extending the blacklist.

```bash
# 1. Call: tumor/normal on chr22 with 4 threads, bundled blacklist
svaba run -t tumor.bam -n normal.bam -G ref.fa -a my_run -p 4 \
          -k chr22 \
          --blacklist tracks/hg38.combined_blacklist.bed

# 2. Post-process: merge per-thread BAMs, sort+dedup+index,
#    sort+dedup bps.txt.gz, filter r2c to PASS.
scripts/svaba_postprocess.sh -t 8 -m 4G my_run

# 3. Convert the deduped bps.txt.gz to VCFv4.5 (SV + indel)
svaba tovcf -i my_run.bps.sorted.dedup.txt.gz -b tumor.bam -a my_run
```

A single-sample call drops `-n`. Any number of cases/controls can be
jointly assembled; prefix on the sample ID drives case/control routing
(`t*` = case, `n*` = control).

## Subcommands

SvABA is a multi-tool binary. `svaba help` lists everything. The main
ones:

`svaba run` performs the whole assembly + variant-calling pipeline.
Takes a BAM (or many) plus a reference, a blacklist, and an output ID.
Emits `bps.txt.gz`, per-sample VCFs, `contigs.bam`, `runtime.txt`, and
(with `--dump-reads`) per-thread `*.discordant.bam`,
`*.corrected.bam`, `*.weird.bam`, and `*.r2c.txt.gz`.

`svaba postprocess` sorts, deduplicates, @PG-stamps, and indexes the
per-thread BAMs and the `bps.txt.gz` emitted by `svaba run`. Typically
invoked via `scripts/svaba_postprocess.sh` which also merges per-thread
files and builds the PASS / PASS-somatic r2c subsets.

`svaba tovcf` converts a deduplicated `bps.txt.gz` into VCFv4.5 output
(one SV VCF, one indel VCF; somatic distinguished by the `SOMATIC`
INFO flag). Clean intrachrom events emit as symbolic `<DEL>`/`<DUP>`/
`<INV>`; everything else stays paired BND. Input is assumed already
sorted/deduped by `svaba_postprocess.sh` — use `--dedup` to opt back
into the legacy internal dedup.

The `SOMATIC` flag is stamped when a record's somatic LOD
(`INFO/SOMLOD`) meets or exceeds a configurable cutoff; the default
is **1.0**, which is a reasonable "confident somatic" gate and keeps
marginal calls out of the somatic set. Tune it with `--somlod`:

```bash
svaba tovcf -i deduped.bps.txt.gz -b tumor.bam -a my_run            # default somlod >= 1
svaba tovcf ... --somlod 2.0                                         # stricter: only strong somatic calls
svaba tovcf ... --somlod 0.0                                         # permissive: flag anything with LO_s >= 0
```

The raw score is always written to `INFO/SOMLOD` regardless of the
cutoff, so downstream `bcftools view -i 'INFO/SOMLOD >= 3'` still
works if you want to re-threshold after the fact without regenerating
the VCF.

`svaba refilter` re-runs the LOD cutoffs / PASS logic on an existing
`bps.txt.gz` with new thresholds, regenerating VCFs without rerunning
assembly. Use it when you want to tune sensitivity/specificity after
the fact.

## Output files

`${ID}.bps.txt.gz` is the primary output — one row per breakpoint,
with a v3 schema of 52 core columns + per-sample blocks (see
`BreakPoint::header` for column names). The 52nd column is a unique
`bp_id` of the form `bpTTTNNNNNNNN` that joins back to the BAM aux
tags and the VCF `EVENT=` field. `${ID}.contigs.bam` holds every
assembled contig, `${ID}.runtime.txt` holds per-region timing, and
`${ID}.log` carries the run log.

The VCF files (`${ID}.sv.vcf.gz`, `${ID}.indel.vcf.gz`) declare
`VCFv4.5`, use symbolic alleles where unambiguous, and carry the
canonical scoring in INFO: `MAXLOD` (variant-vs-error, per-sample
max), `SOMLOD` (somatic LLR), `SOMATIC` (flag), and `SVCLAIM`
(evidence class). VCF `QUAL` defaults to `.` — filter on `FILTER=PASS`
or the two LOD fields, not QUAL. See `CLAUDE.md` for the full scoring
model.

Opt-in outputs (gated behind `--dump-reads`): `${ID}.r2c.txt.gz` is a
structured, re-plottable dump of every contig + its r2c-aligned reads;
`${ID}.corrected.bam` / `${ID}.weird.bam` / `${ID}.discordant.bam`
carry per-read evidence streams. These can run to tens of GB on deep
samples, so they're off by default.

## Post-processing pipeline

`svaba run` emits per-thread, unsorted BAMs and a raw per-thread
`r2c.txt.gz`. Merge + sort + dedup + filter them with one command:

```bash
scripts/svaba_postprocess.sh -t 8 -m 4G my_run
```

Five idempotent steps: merge per-thread BAMs and r2c files, run
`svaba postprocess` (sort + stream-dedup + @PG-stamp + index), sort
and dedup `bps.txt.gz` with PASS filter, emit PASS / PASS-somatic
subsets of `r2c.txt.gz`, and (optional) demux the BAMs by source
prefix. See `CLAUDE.md` for the full flag surface.

## Viewers

All-HTML, no server, drop files in from `file://`. Entry point:
`viewer/index.html`. The primary viewer is `bps_explorer.html` —
sortable `bps.txt.gz` table with chip filters, per-sample detail
panel, log10 histograms for somlod/maxlod/span, and click-to-IGV
navigation. `r2c_explorer.html` re-plots the structured r2c TSV in
browser (contig ruler, fragment rows, indel `||` marker rows with
labels, per-read gap-expanded CIGAR, bp_id filter dropdown).
`runtime_explorer.html` visualizes `runtime.txt`; `comparison.html`
does side-by-side diffs of two runs.

`viewer/alignments_viewer.html` still renders the legacy
`alignments.txt.gz` ASCII output, but new runs don't produce that
file — `r2c_explorer.html` is the replacement.

## Blacklists

`tracks/hg38.combined_blacklist.bed` is the ready-made blacklist;
feed it to `svaba run --blacklist`. It is a regeneratable union of
component BED files in `tracks/` (ENCODE, high-runtime regions, manual
curations, simple repeats, non-standard contigs, and a
low-mappability blacklist derived from paired mosdepth runs). See
`tracks/README.md` for the recipe.

## Recipes

Germline-only. Raise the mate-lookup threshold so only larger
discordant clusters trigger a cross-region lookup — more conservative
and usually appropriate for a single-sample germline run where we
don't have a built-in control to guard against mapping artifacts
(a tumor would typically see smaller, subclonal clusters we want
to chase down). Add `--single-end` to disable mate lookups entirely
if you want to be even more conservative:

```bash
svaba run -t germline.bam -p 8 --mate-min 6 -a germline_run \
          -G ref.fa \
          --blacklist tracks/hg38.combined_blacklist.bed
```

Targeted assembly over a list of regions (BED, chr, or IGV-style
string):

```bash
svaba run -t sample.bam -k targets.bed -a exome_cap -G ref.fa
svaba run -t sample.bam -k chr17:7,541,145-7,621,399 -a TP53 -G ref.fa
```

Dump all per-read evidence (large output, only for debugging a
specific call):

```bash
svaba run -t sample.bam -G ref.fa -a debug_run --dump-reads
scripts/svaba_postprocess.sh -t 8 -m 4G debug_run
```

Snapshot where a long-running job currently is:

```bash
tail my_run.log
```

## Further reading

`CLAUDE.md` at the repo root is the crash-safety-net doc — conventions,
file landmarks, build-system quirks (the `-O2` hardcoding in
submodules), the somatic LOD model, the postprocess pipeline details,
performance notes, and open investigations. Always update `CLAUDE.md`
when you change something non-obvious.

`scripts/svaba_local_function.sh` is a sourceable bash helper library
with `svaba_*` utilities for grepping contigs, following a bp_id
through the output set, opening locations in IGV, etc. Source it from
your shell rc file to use.

## Issues and contact

Please file bug reports, feature requests, and questions on the GitHub
issues tracker: https://github.com/walaj/svaba/issues.

## Attributions

SvABA is developed by Jeremiah Wala in the Rameen Berkoukhim lab at
Dana-Farber Cancer Institute, in collaboration with the Cancer Genome
Analysis team at the Broad Institute. Particular thanks to Cheng-Zhong
Zhang (HMS DBMI) and Marcin Imielinski (Weill Cornell / NYGC).

Additional thanks to Jared Simpson for SGA, Heng Li for htslib and
BWA, and the SeqLib contributors.

The SvABA 2.0 release and its documentation were built with the
assistance of OpenAI Codex and Anthropic Claude, with extensive
human-in-the-loop review, testing, and decision-making throughout.

[seqlib]: https://github.com/walaj/SeqLib
