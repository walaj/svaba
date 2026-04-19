# svaba tracks/

This directory holds the BED files consumed by `svaba run --blacklist`
(and by the scripts that build the combined blacklist from component
files). Everything here is a regeneratable artifact — don't hand-edit
`hg38.combined_blacklist.bed` or `lowmap*perc.bed`; edit the upstream
component or re-run the generator script.

## `hg38.combined_blacklist.bed` — the master blacklist

This is the file you point svaba at. It is the union of several
component blacklists, produced by `scripts/combine_blacklists.sh` from
the components in this same directory. The three primary content
sources are the manual blacklist (`hg38.manual.blacklist.bed`),
empirically-slow-to-assemble regions (`hg38.high_runtime.bed`), and
UCSC RepeatMasker simple repeats (`hg38.rmsk.simple_repeat.bed`). The
non-standard-contig file (`hg38.nonstd_chr.blacklist.bed`) and the
ENCODE-style high-signal blacklist (`hg38.blacklist.sorted.bed`) are
typically folded in as well.

To refresh the combined file after editing any component:

```bash
./scripts/combine_blacklists.sh --merge --clip tracks/chr \
    -o tracks/hg38.combined_blacklist.bed \
    tracks/hg38.manual.blacklist.bed \
    tracks/hg38.high_runtime.bed \
    tracks/hg38.rmsk.simple_repeat.bed
```

`--merge` runs `bedtools merge` to collapse overlapping intervals,
with a comma-separated list of contributing source tags in col 4.
`--clip` clips any oversize-end sentinels (e.g. synthetic
`end=250000000` rows in some BED files) to real contig lengths.
See the header comment in `combine_blacklists.sh` for the full flag
surface.

## `lowmap{30,50}perc.bed` — low-mappability blacklists

These are blacklists of regions where the vast majority of reads at a
locus fail a mapping-quality threshold — i.e. the locus is effectively
"unmapped" for any caller that only trusts uniquely-placed reads.
They're derived from paired `mosdepth` runs over known-non-tumor BAMs
using `scripts/mosdepth_lowmapq_blacklist.sh`, then intersected across
two independent non-tumor samples so only regions confirmed as low-
mappability in both samples survive.

The number in the filename is the ratio cutoff used: `lowmap30perc.bed`
keeps bins where fewer than 30% of reads pass MAPQ (stricter), and
`lowmap50perc.bed` keeps bins below 50% (the default). The generation
recipe has three steps per sample, plus a final intersection.

### Step 1 — mosdepth coverage per sample, unfiltered vs MAPQ-filtered

Run mosdepth twice per BAM, differing only in `--mapq`. The
unfiltered pass (`--mapq 0`) counts every alignment; the filtered pass
(`--mapq 20`) counts only uniquely-placed reads. `--by 300` bins the
coverage into 300 bp windows, which is finer than the svaba
blacklist's natural resolution and gives us room to merge later.

```bash
mosdepth --by 300 --fast-mode --mapq  0 unf_T0  NP_WGS_Normal_DNA.recal.bam
mosdepth --by 300 --fast-mode --mapq 20 filt_T0 NP_WGS_Normal_DNA.recal.bam
mosdepth --by 300 --fast-mode --mapq  0 unf_T2  <second_non_tumor>.recal.bam
mosdepth --by 300 --fast-mode --mapq 20 filt_T2 <second_non_tumor>.recal.bam
```

Each mosdepth invocation produces a `<prefix>.regions.bed.gz` with
`chrom, start, end, mean_depth` columns.

### Step 2 — per-sample ratio + flagging

`mosdepth_lowmapq_blacklist.sh` joins the two mosdepth outputs bin by
bin, computes `filt/unf` (the fraction of reads that pass MAPQ), and
flags bins below a threshold where `unf` is high enough to be
statistically meaningful. `--max-ratio` controls strictness: `0.5`
means "flag bins where fewer than half of reads pass MAPQ" (default),
`0.3` is stricter, `0.7` is looser.

```bash
~/git/svaba/scripts/mosdepth_lowmapq_blacklist.sh --max-ratio 0.5 --out lowmap0.5_T0 unf_T0.regions.bed.gz filt_T0.regions.bed.gz
~/git/svaba/scripts/mosdepth_lowmapq_blacklist.sh --max-ratio 0.5 --out lowmap0.5_T2 unf_T2.regions.bed.gz filt_T2.regions.bed.gz
~/git/svaba/scripts/mosdepth_lowmapq_blacklist.sh --max-ratio 0.7 --out lowmap0.7_T2 unf_T2.regions.bed.gz filt_T2.regions.bed.gz
```

Each run emits `<prefix>.flagged.bed` (raw per-bin flags),
`<prefix>.flagged.merged.bed` (intervals merged within 1 kb), and a
`<prefix>.summary.txt` with totals.

### Step 3 — cross-sample intersection

The per-sample `lowmap0.Xp_T*.flagged.merged.bed` files are noisier
than we want as a blacklist input. Taking the `bedtools intersect` of
two independent non-tumor samples gives a more conservative result:
only regions flagged as low-mappability in BOTH samples survive.

```bash
bedtools intersect -a lowmap0.5_T0.flagged.merged.bed -b lowmap0.5_T2.flagged.merged.bed > lowmap50perc.bed
bedtools intersect -a lowmap0.3_T0.flagged.merged.bed -b lowmap0.3_T2.flagged.merged.bed > lowmap30perc.bed
```

Both samples must be non-tumor — tumor BAMs have structural variation
that can pile up multi-mappers at real break sites, which we don't
want mistaken for low-mappability and blacklisted.

Intersecting two independent non-tumor samples gets you a much smaller
(and more trustworthy) blacklist than either sample's raw output. Per
the v3 mosdepth helper's summary: typical single-sample flagged volume
on hg38 WGS is 200-400 Mb; after intersection with a second non-tumor
sample, the `lowmap50perc.bed` settles in the ~30-50 Mb range, and
`lowmap30perc.bed` shrinks further.

## Other component files

`hg38.blacklist.sorted.bed` is the ENCODE-style high-signal blacklist.
`hg38.high_runtime.bed` is a list of regions empirically slow to
assemble; curated from svaba runtime.txt tails. `hg38.manual.blacklist.bed`
is the curated ad-hoc bad-list. `hg38.nonstd_chr.blacklist.bed` is
full-contig entries for every chrUn/*_decoy/*_alt/*_random/chrEBV/HLA-*
contig in the reference (regenerated from `tracks/chr` which is a
GRCh38 fasta-header dump). `hg38.rmsk.simple_repeat.bed` is UCSC
RepeatMasker simple-repeat regions.

`genome.hg38.sorted.bed`, `hg38.bed`, and
`hg38_arms_excl_centromeres.bed` are reference-coordinate helpers used
by downstream analysis scripts, not by `svaba run` directly.
`region_generator.R` is an R script used to regenerate some of the
above.

## Adding a new component to the combined blacklist

Put the new `.bed` file in this directory, then re-run
`combine_blacklists.sh` with it added to the argument list. The script
tolerates any well-formed 3-column BED. If the new file has a 4th
column with per-interval labels, those get preserved in the merged
output's label column.
