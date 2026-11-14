# SvABA — developer notes

Build tuning, debugging harnesses, and internals. User-facing docs are in
[README.md](./README.md); the authoritative internals reference (conventions,
file landmarks, the somatic LOD model, the postprocess pipeline, perf notes,
open investigations) is `CLAUDE.md`.

## Build tuning

Build type defaults to `RelWithDebInfo` (`-O2 -g -DNDEBUG
-fno-omit-frame-pointer`). The single biggest free-perf knob is pushing the
vendored assemblers to `-O3 -mcpu=native` — they hardcode `-O2` in their own
Makefiles and CMake's build type doesn't reach them (`libbwa.a` ≈ 38% of
wall-time, `libfml.a` ≈ 27%):

```bash
make -C SeqLib/bwa        clean
make -C SeqLib/fermi-lite clean
make -C SeqLib/bwa        -j CFLAGS="-g -O3 -mcpu=native -fno-omit-frame-pointer -Wall -Wno-unused-function"
make -C SeqLib/fermi-lite -j CFLAGS="-g -O3 -mcpu=native -fno-omit-frame-pointer -Wall -Wno-unused-function"
cd build && make -j
```

Sanity-check the flags that actually reached the compiler:
`make VERBOSE=1 2>&1 | grep -oE -- "-O[0-9sg]" | sort | uniq -c`. See the
"Build system" section of `CLAUDE.md` for the full story on the `-O2`
hardcoding.

### Assembler selection

By default svaba assembles with **fermi-lite** (ropebwt2-based, the faster
path). The vendored SGA (String Graph Assembler) is still fully supported but
off by default — switch it on at compile time:

```bash
cmake .. -DSVABA_ASSEMBLER_FERMI=0
make -j
```

The choice is a single compile-time symbol in
`src/svaba/SvabaAssemblerConfig.h`. Runtime logs report the active assembler
under `svaba::kAssemblerName`.

## jemalloc (Linux, high thread count)

On Linux, running with **jemalloc** typically wins 10–20% wall-time at `-p 16+`
(large WGS, tumor/normal): svaba's per-thread assembly loop generates heavy
malloc/free traffic and glibc's default allocator serializes on its arena locks
under that contention pattern. Two ways to enable it.

Compile-time link (no `LD_PRELOAD` needed):

```bash
cmake .. -DUSE_JEMALLOC=ON
make -j
```

Or the drop-in wrapper `./svaba_jemalloc`, which `LD_PRELOAD`s the system
`libjemalloc.so.2` and exec's svaba:

```bash
./svaba_jemalloc run -t tumor.bam -n normal.bam -G ref.fa -a my_run -p 16 \
                     --blacklist tracks/hg38.combined_blacklist.bed
```

Install jemalloc first (`apt install libjemalloc2` / `yum install jemalloc` /
`dnf install jemalloc`). The wrapper auto-detects the library under standard
distro paths; override with `JEMALLOC_LIB=/path/to/libjemalloc.so.2` and, if
svaba isn't on `$PATH`, `SVABA=/path/to/svaba`.

For very high concurrency (`-p 24+`), also pass jemalloc's tuning knobs:

```bash
MALLOC_CONF=background_thread:true,narenas:24,dirty_decay_ms:10000 \
    ./svaba_jemalloc run -p 24 ...
```

`narenas` should match or modestly exceed the thread count;
`background_thread:true` reclaims dirty pages off the hot path.

**macOS users should not use jemalloc.** Apple's native libmalloc (with its
nanomalloc fast path) plus the `DYLD_INSERT_LIBRARIES` mechanism run 5–10×
slower than system malloc on this workload in our A/B tests. The wrapper refuses
to run on Darwin. Stick with the system allocator on macOS.

## Debugging recipes

Three zero-cost compile-time trace/restrict systems. All are `#ifdef`-guarded
and compile to no-ops when not defined — no runtime cost in a normal build. See
`src/svaba/SvabaDebug.h` for the macro definitions.

### Trace a specific read through the entire pipeline

`-DSVABA_TRACE_READ` emits per-decision-point stderr output for a single read
(by QNAME) from BAM ingestion through output tagging:

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

Shows assembly filtering (TP1–TP5), variant identification (TP6), per-read
split-coverage scoring (TP8–TP11), breakpoint confidence (TP13–TP16, TP23).

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

# From the r2c.db — all reads for a contig
sqlite3 results.r2c.db "SELECT read_id FROM reads WHERE cname='c_fermi_chr2_215869501_215894501_13C';"
```

### Restrict assembly to reads containing a specific kmer

```bash
cmake -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_CXX_FLAGS='-DSVABA_KMER_RESTRICT="\"CCATGCAGAGTGTTGAAGAAAAGGC\""'
cmake --build build -j$(nproc)
```

After BFC error correction, only reads whose corrected sequence contains the
specified kmer (or its reverse complement) survive. All other reads get
`to_assemble = false` — they won't enter fermi assembly, r2c alignment,
`corrected.bam`, or scoring. The log prints kept/dropped counts per region.

Use case: you see a suspicious contig and want to know whether assembly still
produces it when restricted to reads from a particular locus. Pick a 25-mer
unique to that locus, compile with `SVABA_KMER_RESTRICT`, and run on the same
region:

```bash
svaba run -t tumor.bam -n normal.bam -G ref.fa -k chr11:5000000-5100000 \
          -a kmer_test --dump-reads -p 1
```

If the chimeric contig still assembles, the kmer-carrying reads are sufficient
to produce it. If it disappears, reads from elsewhere were required. Can be
combined with read/contig tracing:

```bash
cmake -B build \
      -DCMAKE_CXX_FLAGS='-DSVABA_KMER_RESTRICT="\"CCATGCAGAGTGTTGAAGAAAAGGC\"" \
                          -DSVABA_TRACE_CONTIG="\"c_fermi_chr11_5000001_5100001_3C\""'
```

### Disable the r2c-vs-native gate (debugging only)

```bash
cmake ... -DCMAKE_CXX_FLAGS='-DSVABA_R2C_NATIVE_GATE=0'
```

Restores old behavior where any r2c-spanning read credits as a variant
supporter. Will reintroduce false-positive somatic calls from homology-trap
reads.

## Alt-contig demotion

BWA-MEM sometimes places a contig fragment on an alt contig (e.g.
`chr11_JH159136v1_alt`) when chr11 proper has an equally good alignment. If the
alt later gets blacklisted, the real breakpoint is lost.

svaba requests secondary alignments for contigs and runs a post-alignment step
(`preferStandardChromosomes`): for each primary/supplementary fragment on a
non-standard chromosome (ChrID > `maxMateChrID`, default 23), it looks for a
secondary alignment on a standard chromosome that:

- covers ≥ 80% of the same query range (reciprocal overlap), and
- has AS ≥ 95% of the non-standard alignment's AS.

If found, the standard-chr alignment is promoted (gets the
primary/supplementary flag) and the non-standard one is demoted to secondary.
The contig trace (`SVABA_TRACE_CONTIG`) logs each swap with both chromosome IDs
and alignment scores.

Controlled by `--max-mate-chr` (same constant as mate-region lookup).
`--non-human` sets it to `-1`, disabling alt-demotion entirely (no chromosome is
"non-standard").
