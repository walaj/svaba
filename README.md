## *SvABA* — Structural variation and indel analysis by assembly

SvABA (formerly *Snowman*) is an SV and indel caller for short-read BAMs.
It performs genome-wide local assembly, realigns contigs with BWA-MEM, and
scores variants by reassembled read support. Tumor/normal, trios, and
single-sample modes are supported; variants are emitted as VCF plus a verbose
tab-delimited companion (`bps.txt.gz`) carrying the full per-sample evidence.

**License:** [GNU GPLv3](./LICENSE). Uses the [SeqLib][seqlib] API for BAM I/O,
BWA-MEM alignment, interval trees, and the assembly front-end.

For debugging recipes, build tuning, and internals, see **[README.dev.md](./README.dev.md)**
and `CLAUDE.md`.

## Install

| Package      | Required? | Purpose                                     |
|--------------|-----------|---------------------------------------------|
| CMake ≥ 3.14 | yes       | build system                                |
| htslib       | yes       | BAM/CRAM/VCF I/O                            |
| zlib         | yes       | gzip in/out                                 |
| pthread      | yes       | worker threads                              |
| BZip2, lzma  | yes       | htslib compression deps                     |
| sqlite3      | optional  | enables `--dump-reads` `r2c.db` output      |
| jemalloc     | optional  | recommended on Linux at `-p 16+`            |

```bash
git clone --recursive https://github.com/walaj/svaba
cd svaba && mkdir build && cd build
cmake .. && make -j
```

The binary lands at `build/svaba`. For a non-standard htslib, pass
`-DHTSLIB_DIR=/path/to/htslib-1.xx`. Install system-wide with `make install`
(prefix `/usr/local`, override with `-DCMAKE_INSTALL_PREFIX`).

Build type defaults to `RelWithDebInfo` (`-O2 -g -DNDEBUG`). See
[README.dev.md](./README.dev.md) for `-O3 -mcpu=native` tuning (5–15% faster),
the fermi-lite/SGA assembler switch, and jemalloc (10–20% faster on Linux at
high thread counts).

When sqlite3 isn't found, svaba builds fine but `--dump-reads` skips the
`${ID}.r2c.db` file (the other `--dump-reads` outputs are unaffected). To
enable it: `brew install sqlite` / `apt install libsqlite3-dev` /
`dnf install sqlite-devel`, then re-run cmake.

## Quick start

Three steps: call, post-process, convert to VCF. The bundled combined
blacklist is strongly recommended — it keeps svaba out of pileup /
high-complexity regions that otherwise dominate wall-clock with no real calls.

```bash
# 1. Call: tumor/normal on chr22 with 4 threads, bundled blacklist
svaba run -t tumor.bam -n normal.bam -G ref.fa -a my_run -p 4 -k chr22 \
          --blacklist tracks/hg38.combined_blacklist.bed

# 2. Post-process (one command): merge per-thread BAMs, sort+dedup+index,
#    sort+dedup bps.txt.gz, stamp PASS reads into r2c.db.
svaba postprocess -i my_run -t 8 -m 4G

# 3. Convert the deduped bps.txt.gz to VCFv4.5 (SV + indel)
svaba tovcf -i my_run.bps.sorted.dedup.txt.gz -b tumor.bam -a my_run
```

A single-sample call drops `-n`. Any number of cases/controls can be jointly
assembled; the sample-ID prefix drives routing (`t*` = case, `n*` = control).

## Subcommands

SvABA is a multi-tool binary; `svaba help` lists everything.

- **`svaba run`** — the whole assembly + variant-calling pipeline. Emits
  `bps.txt.gz`, per-sample VCFs, `contigs.bam`, `runtime.txt`, and (with
  `--dump-reads`) per-thread `*.discordant.bam`, `*.corrected.bam`, `*.r2c.db`.
- **`svaba postprocess`** — one-command post-processing: merges per-thread
  BAMs and `r2c.db` files, coordinate-sorts + stream-dedups + `@PG`-stamps +
  indexes the BAMs, sorts/dedups `bps.txt.gz` (writing PASS / PASS-somatic
  subsets), and stamps PASS reads into `r2c.db`. Six idempotent phases; reruns
  are near-instant.
- **`svaba tovcf`** — converts a deduplicated `bps.txt.gz` into VCFv4.5 (one SV
  VCF, one indel VCF; somatic marked by the `SOMATIC` INFO flag). Clean
  intrachrom events emit as symbolic `<DEL>`/`<DUP>`/`<INV>`; everything else
  stays paired BND. The `SOMATIC` flag is stamped when `INFO/SOMLOD` ≥ a cutoff
  (default **1.0**, tune with `--somlod`). The raw score is always written to
  `INFO/SOMLOD`, so `bcftools view -i 'INFO/SOMLOD >= 3'` re-thresholds without
  regenerating the VCF.
- **`svaba refilter`** — re-runs LOD cutoffs / PASS logic on an existing
  `bps.txt.gz` with new thresholds, regenerating VCFs without re-assembling.
- **`svaba extract-pairs`** — pulls every read pair whose SEQ contains a query
  sequence (or, with `-f bps.txt.gz`, every read carrying a variant's junction
  kmer). See `svaba extract-pairs -h`.

## Output files

`${ID}.bps.txt.gz` is the primary output — one row per breakpoint, with a v4
schema of 53 core columns plus one FORMAT-style block per sample. Full column
reference below.

The VCF files (`${ID}.sv.vcf.gz`, `${ID}.indel.vcf.gz`) declare `VCFv4.5`, use
symbolic alleles where unambiguous, and carry the canonical scoring in INFO:
`MAXLOD` (variant-vs-error, per-sample max), `SOMLOD` (somatic LLR), `SOMATIC`
(flag), `SVCLAIM` (evidence class). VCF `QUAL` defaults to `.` — filter on
`FILTER=PASS` or the two LOD fields, not QUAL.

`${ID}.contigs.bam` holds every assembled contig; `${ID}.runtime.txt` holds
per-region timing; `${ID}.log` is the run log.

Opt-in outputs (behind `--dump-reads`, can run to tens of GB on deep samples):
`${ID}.r2c.db` (queryable SQLite of every contig + its r2c-aligned reads),
`${ID}.corrected.bam` / `${ID}.discordant.bam` (per-read evidence streams).

### `bps.txt.gz` column reference

Columns are 1-indexed (awk `$1` == `chr1`). This is the v4 schema emitted by
`BreakPoint::toFileString`; positions are 1-based (SAM/VCF convention).
`x` marks an empty string field; `.` marks an empty numeric/id field.

| # | Name | Description |
|---|------|-------------|
| 1 | `chr1` | Chromosome of breakend 1 |
| 2 | `pos1` | Position of breakend 1 (1-based) |
| 3 | `strand1` | Orientation of breakend 1 (`+`/`-`) |
| 4 | `chr2` | Chromosome of breakend 2 |
| 5 | `pos2` | Position of breakend 2 (1-based) |
| 6 | `strand2` | Orientation of breakend 2 (`+`/`-`) |
| 7 | `ref` | Reference allele at the breakpoint |
| 8 | `alt` | Alternate allele |
| 9 | `span` | Event span in bp (`-1` for interchromosomal) |
| 10 | `split` | Total split-read support across all samples |
| 11 | `alt_count` | Total alt-read count |
| 12 | `cov` | Read coverage at the breakpoint |
| 13 | `cigar` | Indel support from read CIGAR strings |
| 14 | `cigar_near` | CIGAR support near (but not exactly at) the breakpoint |
| 15 | `dmq1` | Discordant-cluster MAPQ, side 1 |
| 16 | `dmq2` | Discordant-cluster MAPQ, side 2 |
| 17 | `dcn` | Discordant read count, normal |
| 18 | `dct` | Discordant read count, tumor |
| 19 | `mapq1` | Contig-alignment MAPQ, side 1 |
| 20 | `mapq2` | Contig-alignment MAPQ, side 2 |
| 21 | `nm1` | Contig-alignment edit distance (NM), side 1 |
| 22 | `nm2` | Contig-alignment edit distance (NM), side 2 |
| 23 | `as1` | Contig-alignment score (BWA AS), side 1 |
| 24 | `as2` | Contig-alignment score (BWA AS), side 2 |
| 25 | `sub1` | Number of sub-optimal contig alignments, side 1 |
| 26 | `sub2` | Number of sub-optimal contig alignments, side 2 |
| 27 | `homol` | Microhomology at the junction (side-1 forward strand); `x` if none |
| 28 | `insert` | Inserted novel (non-template) sequence at the junction; `x` if none |
| 29 | `repeat` | Repeat-sequence context at the breakpoint; `x` if none |
| 30 | `contig_and_region` | Contig name (`cname`) — the r2c / `r2c.db` join key |
| 31 | `naln` | Number of contig-alignment fragments |
| 32 | `conf` | Confidence / FILTER (`PASS`, `LOWLOD`, `LOWSUPPORT`, …) |
| 33 | `type` | Evidence type (`ASSMB`, `ASDIS`, `DSCRD`, `INDEL`, …) |
| 34 | `qual` | Variant quality |
| 35 | `2ndary` | Secondary-alignment flag |
| 36 | `somatic` | Somatic call flag (`NA` when no normal sample present) |
| 37 | `somlod` | Somatic LOD (`LO_s`), capped at 99 |
| 38 | `maxlod` | Max per-sample `LO` (variant-vs-error) |
| 39 | `dbsnp` | dbSNP overlap (rs id, or `x`) |
| 40 | `contig_conf1` | Contig-alignment confidence, side 1 |
| 41 | `contig_conf2` | Contig-alignment confidence, side 2 |
| 42 | `cpos1` | Contig-relative breakend position, side 1 (`-1` = unset) |
| 43 | `cpos2` | Contig-relative breakend position, side 2 (`-1` = unset) |
| 44 | `lmatch` | Left flanking match length |
| 45 | `rmatch` | Right flanking match length |
| 46 | `scov1` | Split-coverage lower bound on the contig |
| 47 | `scov2` | Split-coverage upper bound on the contig |
| 48 | `local1` | Per-end LocalAlignment enum (0–3), side 1 |
| 49 | `local2` | Per-end LocalAlignment enum (0–3), side 2 |
| 50 | `ctglen` | Contig length |
| 51 | `flipped` | Contig orientation flag (1 = contig flipped relative to side 1) |
| 52 | `bp_id` | Unique per-BP identifier `bpTTTNNNNNNNN` (thread `TTT`, counter `NNNNNNNN`). Joins to `bps.txt` col 52, the BAM `bi:Z` tag, `r2c.db`'s `split_bps`/`disc_bps`, and the VCF `EVENT=` field. `.` if unset |
| 53 | `jxn_kmer` | 20 bp contig sequence spanning the breakend junction — a read-search query for `svaba extract-pairs -f bps.txt.gz`. `.` when no precise junction exists |
| 54+ | *per-sample block* | One `FORMAT`-style block per BAM (order matches the header row) — see below |

**Per-sample block** (columns 54+, one per BAM), colon-delimited
`GT:AD:DP:SR:DR:GQ:PL:LO:LO_n`:

| Sub | Name | Description |
|-----|------|-------------|
| GT | genotype | Genotype call |
| AD | alt depth | Alt-supporting read count (for indels, `max(alt, cigar)`) |
| DP | depth | Total read coverage |
| SR | split reads | Split-read support count |
| DR | discordant reads | Discordant-read support (**indels:** this field holds the CIGAR-indel count instead) |
| GQ | genotype quality | Phred genotype quality |
| PL | phred likelihoods | Comma-separated genotype likelihoods |
| LO | log-odds | Variant vs error (per-sample); `maxlod` is the max of these |
| LO_n | log-odds ref | Confidence the site is REF-only in this sample (the germline-vs-somatic discriminant; the normal's `LO_n` drives `somlod`) |

`scripts/svaba_local_function.sh::svaba_bps_cols` prints this same reference
from the shell.

## Viewers

All-HTML, no server — drop files in from `file://`. Entry point:
`docs/index.html`. Primary viewer is `bps_explorer.html` (sortable
`bps.txt.gz` table, chip filters, per-sample detail, histograms, click-to-IGV).
`r2c_db_explorer.html` loads a `${ID}.r2c.db` and plots each contig with its
r2c-aligned reads. `runtime_explorer.html` visualizes `runtime.txt`;
`comparison.html` diffs two runs; `learn_explorer.html` plots insert-size
distributions.

## Blacklists

`tracks/hg38.combined_blacklist.bed` is the ready-made blacklist; feed it to
`svaba run --blacklist`. It's a regeneratable union of component BEDs in
`tracks/` (ENCODE, high-runtime regions, manual curations, simple repeats,
non-standard contigs, low-mappability). See `tracks/README.md` for the recipe.

## Recipes

**Germline-only.** Raise the mate-lookup threshold so only larger discordant
clusters trigger a cross-region lookup (add `--single-end` to disable mate
lookups entirely):

```bash
svaba run -t germline.bam -p 8 --mate-min 6 -a germline_run -G ref.fa \
          --blacklist tracks/hg38.combined_blacklist.bed
```

**Targeted assembly** over a region list (BED, chr, or IGV-style string):

```bash
svaba run -t sample.bam -k targets.bed -a exome_cap -G ref.fa
svaba run -t sample.bam -k chr17:7,541,145-7,621,399 -a TP53 -G ref.fa
```

**Non-human genome** (mouse, zebrafish, …). By default svaba assumes a human
reference (mate lookups skip past chrY, insert-size learning samples chr1–chrY).
`--non-human` removes those gates:

```bash
svaba run -t mouse.bam -G mm39.fa -a mouse_run --non-human -p 16
```

Fine-tune individually with `--max-mate-chr N` and `--min-mate-mapq N`.

## Issues and contact

File bug reports, feature requests, and questions on GitHub:
https://github.com/walaj/svaba/issues.

## Attributions

SvABA is developed by Jeremiah Wala in the Rameen Berkoukhim lab at Dana-Farber
Cancer Institute, in collaboration with the Cancer Genome Analysis team at the
Broad Institute. Particular thanks to Cheng-Zhong Zhang (HMS DBMI) and Marcin
Imielinski (Weill Cornell / NYGC).

Additional thanks to Jared Simpson for SGA, Heng Li for htslib and BWA, and the
SeqLib contributors.

The SvABA 2.0 release and its documentation were built with the assistance of
OpenAI Codex and Anthropic Claude, with extensive human-in-the-loop review,
testing, and decision-making throughout.

[seqlib]: https://github.com/walaj/SeqLib
