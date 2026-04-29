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

To link jemalloc at compile time (recommended on Linux for `-p 16+`
runs — eliminates allocator contention, no `LD_PRELOAD` needed):

```bash
cmake .. -DUSE_JEMALLOC=ON
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

### jemalloc (Linux, high thread count)

On Linux, running svaba with **jemalloc** as the allocator typically
wins 10–20 % wall-time at `-p 16+` (large WGS, tumor/normal). svaba's
per-thread assembly loop generates heavy malloc/free traffic and
glibc's default allocator serializes on its arena locks under that
contention pattern; jemalloc's thread-local arenas sidestep the
problem. The repo ships a drop-in wrapper at `./svaba_jemalloc` that
`LD_PRELOAD`s the system `libjemalloc.so.2` and then exec's svaba:

```bash
# one-liner replacement for `svaba run ...`:
./svaba_jemalloc run -t tumor.bam -n normal.bam -G ref.fa -a my_run -p 16 \
                     --blacklist tracks/hg38.combined_blacklist.bed
```

Install jemalloc first (`apt install libjemalloc2`, `yum install
jemalloc`, `dnf install jemalloc`). The wrapper auto-detects the
library under the standard distro paths; if yours is elsewhere, set
`JEMALLOC_LIB=/path/to/libjemalloc.so.2`. If the svaba binary isn't on
your `$PATH`, set `SVABA=/path/to/svaba`.

For very high concurrency (`-p 24+`), also pass jemalloc's own tuning
knobs via `MALLOC_CONF`:

```bash
MALLOC_CONF=background_thread:true,narenas:24,dirty_decay_ms:10000 \
    ./svaba_jemalloc run -p 24 ...
```

`narenas` should match or modestly exceed the thread count;
`background_thread:true` asks jemalloc to reclaim dirty pages on a
background thread instead of in the hot path.

**macOS users should not use jemalloc.** Apple's native libmalloc
(with its nanomalloc fast path for small allocations) and the
DYLD_INSERT_LIBRARIES mechanism combine to run 5–10× slower than
system malloc on this exact workload in our A/B tests. The wrapper
refuses to run on Darwin for that reason. Stick with the system
allocator on macOS.

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

## Debugging recipes

### Trace a specific read through the entire pipeline

Compile with `-DSVABA_TRACE_READ` to get per-decision-point stderr output
for a single read (by QNAME) from BAM ingestion through to output tagging.
Zero runtime cost when not compiled in (all trace macros compile to no-ops).

```bash
cmake -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_CXX_FLAGS='-DSVABA_TRACE_READ="\"LH00306:129:227V5CLT4:6:1204:38807:7191\""'
cmake --build build -j$(nproc)

# Run on just the region containing the read (faster iteration)
svaba run -t tumor.bam -n normal.bam -G ref.fa -p 4 \
    -k chr2:215000000-216000000 --dump-reads -a trace_run 2> read_trace.log

grep READ_TRACE read_trace.log
```

The trace covers every gate the read hits:

| Stage | Log prefix | What it tells you |
|-------|-----------|-------------------|
| BAM read & filter | (various, from SvabaBamWalker) | Dup/QC/blacklist gates, rule_pass, NM salvage, adapter filter, quality trim, buffer admission |
| BFC error correction | `BFC corrected:` | Whether BFC changed the sequence and by how much |
| R2C alignment | `R2C SKIP:` / `R2C HIT:` | Perfect-ref-match skip, or which contigs the read aligned to (AS score, CIGAR) |
| Native realignment | `NATIVE_REALIGN:` | Whether the corrected seq was re-aligned to ref or reused the BAM CIGAR |
| Split coverage scoring | `SPLIT_COV enter` | Entry into per-BP scoring with r2c coords and breakpoint positions |
| R2C vs native gate | `SPLIT_COV TP8` | The critical score comparison: r2c_score, native_score, both CIGARs, NM values |
| Indel-at-break check | `SPLIT_COV TP9` | Whether an r2c CIGAR indel lands at the breakpoint |
| Span check | `SPLIT_COV TP10` | Whether the read spans the breakpoint(s), with exact coord bounds |
| DEL-covers-break | `SPLIT_COV TP11` | Whether an r2c deletion masks the breakpoint |
| Final verdict | `SPLIT_COV CREDITED` / `NOT CREDITED` / `SKIPPED` | Did this read count as a variant supporter? |
| Output tagging | `OUTPUT TAG bi:Z` | BP id stamped on the read, confidence, somatic status |

Example trace (abridged) for a read supporting a somatic 1bp deletion:

```
[READ_TRACE:SvabaBamWalker.cpp:203] read=LH00306:129:... flag=163 mapq=60 pos=chr2:215869800 ...
[READ_TRACE:SvabaBamWalker.cpp:236] PASS mr.isValid rule_pass=1
[READ_TRACE:SvabaBamWalker.cpp:340] ADDED to read buffer (n=847)
[READ_TRACE:SvabaRegionProcessor.cpp:377] BFC corrected: changed=NO ...
[READ_TRACE:SvabaRegionProcessor.cpp:808] R2C HIT: contig=c_fermi_chr2_... AS=150 rc=0 CIGAR=75M1D74M
[READ_TRACE:SvabaRegionProcessor.cpp:978] NATIVE_REALIGN: reusing BAM CIGAR (seq unchanged by BFC)
[READ_TRACE:BreakPoint.cpp:379] SPLIT_COV enter contig=c_fermi_chr2_... r2c_start=340 r2c_end=489 ...
[READ_TRACE:BreakPoint.cpp:476] SPLIT_COV TP8 PASS r2c>native r2c_score=150 native_score=143 ...
[READ_TRACE:BreakPoint.cpp:601] SPLIT_COV TP10 span check issplit1=1 issplit2=1 one_split=1
[READ_TRACE:BreakPoint.cpp:697] SPLIT_COV CREDITED as variant supporter sample=t000 tumor=1
[READ_TRACE:SvabaRegionProcessor.cpp:1346] OUTPUT TAG bi:Z bp_id=bp00100000042 confidence=PASS somatic=1
```

### Trace a specific contig

```bash
cmake -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_CXX_FLAGS='-DSVABA_TRACE_CONTIG="\"c_fermi_chr2_215869501_215894501_13C\""'
cmake --build build -j$(nproc)
```

Shows assembly filtering (TP1–TP5), variant identification (TP6),
per-read split-coverage scoring (TP8–TP11), breakpoint confidence
(TP13–TP16, TP23).

### Trace both a read AND a contig simultaneously

```bash
cmake -B build \
      -DCMAKE_CXX_FLAGS='-DSVABA_TRACE_READ="\"LH00306:129:227V5CLT4:6:1204:38807:7191\"" \
                          -DSVABA_TRACE_CONTIG="\"c_fermi_chr2_215869501_215894501_13C\""'
cmake --build build -j$(nproc)
```

Independent prefixes (`[READ_TRACE:...]` vs `[TRACE:...]`) — grep for either.

### Trace ALL reads or ALL contigs (very noisy)

```bash
cmake ... -DCMAKE_CXX_FLAGS='-DSVABA_TRACE_ALL_READS=1'   # every read
cmake ... -DCMAKE_CXX_FLAGS='-DSVABA_TRACE_ALL=1'          # every contig
```

Best combined with a small `-k` region.

### Finding a read's QNAME to trace

```bash
# From bps.txt.gz — get contig name (col 30) and bp_id (col 52)
zcat results.bps.txt.gz | awk -F'\t' '$1=="chr2" && $2 > 215869000 && $2 < 215870000'

# From the corrected BAM — reads tagged with a specific bp_id
samtools view results.corrected.bam chr2:215869000-215870000 | grep "bi:Z:.*bp00100000042"

# From the r2c TSV — all reads for a contig
zcat results.r2c.txt.gz | awk -F'\t' '$2 == "c_fermi_chr2_215869501_215894501_13C" && $1 == "read"' | cut -f8
```

### Restrict assembly to reads containing a specific kmer

```bash
cmake -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_CXX_FLAGS='-DSVABA_KMER_RESTRICT="\"CCATGCAGAGTGTTGAAGAAAAGGC\""'
cmake --build build -j$(nproc)
```

After BFC error correction, only reads whose corrected sequence
contains the specified kmer (or its reverse complement) survive.
All other reads get `to_assemble = false` — they won't enter fermi
assembly, r2c alignment, corrected.bam, or scoring. The log prints
kept/dropped counts per region.

Use case: you see a suspicious contig and want to know whether
assembly still produces it when restricted to reads from a particular
locus. Pick a 25-mer unique to that locus, compile with
`SVABA_KMER_RESTRICT`, and run on the same region:

```bash
svaba run -t tumor.bam -n normal.bam -G ref.fa -k chr11:5000000-5100000 \
          -a kmer_test --dump-reads -p 1
```

If the chimeric contig still assembles, the kmer-carrying reads are
sufficient to produce it. If it disappears, reads from elsewhere
were required.

Can be combined with read/contig tracing:

```bash
cmake -B build \
      -DCMAKE_CXX_FLAGS='-DSVABA_KMER_RESTRICT="\"CCATGCAGAGTGTTGAAGAAAAGGC\"" \
                          -DSVABA_TRACE_CONTIG="\"c_fermi_chr11_5000001_5100001_3C\""'
```

### Disable the r2c-vs-native gate (for debugging only)

```bash
cmake ... -DCMAKE_CXX_FLAGS='-DSVABA_R2C_NATIVE_GATE=0'
```

Restores old behavior where any r2c-spanning read credits as variant supporter.
Will reintroduce false-positive somatic calls from homology-trap reads.

---

## Alt-contig demotion

BWA-MEM sometimes places a contig fragment on an alt contig (e.g.
`chr11_JH159136v1_alt`) when chr11 proper has an equally good
alignment. If the alt later gets blacklisted, the real breakpoint
is lost.

svaba now requests secondary alignments for contigs and runs a
post-alignment step (`preferStandardChromosomes`): for each
primary/supplementary fragment on a non-standard chromosome
(ChrID > `maxMateChrID`, default 23), it looks for a secondary
alignment on a standard chromosome that:

- Covers ≥ 80% of the same query range (reciprocal overlap)
- Has AS ≥ 95% of the non-standard alignment's AS

If found, the standard-chr alignment is promoted (gets the
primary/supplementary flag) and the non-standard one is demoted
to secondary. The contig trace (`SVABA_TRACE_CONTIG`) logs each
swap with both chromosome IDs and alignment scores.

This is controlled by `--max-mate-chr` (same constant used for
mate-region lookup). `--non-human` sets it to -1, which disables
alt-demotion entirely (no chromosome is "non-standard").

---

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

Non-human genome (e.g. mouse, zebrafish, _C. elegans_). By default
svaba assumes a human reference: mate-region lookups skip chromosomes
past chrY (ChrID > 23), and insert-size learning samples only the
primary assembly (chr1–chrY). `--non-human` removes these gates so
every contig in the reference is treated equally:

```bash
svaba run -t mouse.bam -G mm39.fa -a mouse_run --non-human -p 16
```

The flag sets `maxMateChrID = -1` internally. You can also fine-tune
individually: `--max-mate-chr 19` (mouse has 19 autosomes + X + Y =
ChrIDs 0–20 in a typical reference), `--min-mate-mapq 10` (require
MAPQ ≥ 10 on the primary read before chasing its mate).

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
