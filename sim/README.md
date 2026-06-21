# svaba SV simulator (`sim/`)

A self-contained pipeline to generate a **reproducible structural-variant truth
set** and matching simulated WGS BAMs for validating svaba (or any SV caller).
It places thousands of **non-overlapping** SVs spanning the full size range
(small indels → inter-chromosomal translocations) across the genome, builds a
diploid donor genome, simulates reads with **ART**, aligns with **BWA-MEM**, and
admixes a normal for **tumor-impurity** series.

```
sim/
├── simulate_sv_genome.py   # SVs + diploid donor genome + truth (VCF/BEDPE/TSV)  [stdlib + samtools]
├── run_sim.sh              # genome → ART → BWA-MEM → tumor BAMs at N coverages
├── mix_tumor_normal.sh     # in-silico tumor impurity (purity admixture)
├── benchmark.py            # score a caller's output vs truth (recall/precision)
├── sim_config.sh           # all paths/parameters (edit or export overrides)
├── environment.yml         # conda env: art, bwa, samtools
└── README.md
```

## Requirements
`art_illumina`, `bwa`, `samtools`, `python3` (standard library only — no pysam/numpy).
The reference must be BWA-indexed (`.bwt`) and faidx'd (`.fai`).
```bash
mamba env create -f sim/environment.yml && mamba activate svaba-sim
```

## Quickstart

**Smoke test** (one small chromosome, low coverage — minutes):
```bash
cd sim
CHROMS=chr22 MAXCOV=10 COVERAGES="5 10" \
GEN_ARGS="--n-del-small 60 --n-ins-small 60 --n-del 60 --n-dup 40 --n-inv 40 --n-tra 10 --min-gap 50000" \
OUTDIR=$PWD/smoke ./run_sim.sh
```

**Full genome** (edit `sim_config.sh` for paths, then):
```bash
cd sim && ./run_sim.sh          # writes $OUTDIR/tumor.<cov>x.bam + truth.*
```

**Run svaba against it** and score:
```bash
svaba run -t sim_out/tumor.30x.bam -n ~/Downloads/svaba_compare/blood.recal.bam \
          -G $REF -a sim_eval -p 8
svaba postprocess -i sim_eval
python3 benchmark.py --truth sim_out/truth.bedpe \
          --calls sim_eval.bps.sorted.dedup.txt.gz --pass-only --somatic-only
```

**Interactive ROC / PR viewer:** open `docs/sim_benchmark.html`, drop in
`truth.bedpe` + a `bps.txt[.gz]`, and sweep a score threshold (somlod / maxlod /
qual) to get Precision-Recall and Recall-vs-FalsePositives (FROC) curves, an
operating-point slider, and recall by SV type / size bin. It streams large
gzip and auto-restricts to the truth's chromosomes (so a genome-wide bps from a
real-normal pairing isn't drowned in off-target germline). Tip: for a small,
fast bps, run svaba with `-k <chroms>` matching the simulated `CHROMS`.

## How it works

1. **`simulate_sv_genome.py`** places SVs by rejection sampling: each event
   reserves a padded footprint (`--min-gap`, default 50 kb) so no two SVs
   overlap/stack, none crosses an N-gap, and none falls in `--avoid-bed`
   (pass the svaba blacklist here to skip centromeres etc.). Chromosomes are
   weighted by callable length for even ("thick") coverage. Each SV is het/hom
   (`--hom-frac`) and somatic/germline (`--germline-frac`, default 0 = all
   somatic). It then builds **diploid donor haplotypes** (`tumor.hap1.fa`,
   `tumor.hap2.fa`) by splicing the edits into the reference; reciprocal
   translocations become derivative chromosomes. Truth is written in reference
   coordinates as `truth.vcf`, `truth.bedpe`, `truth.sv_table.tsv`.

2. **`run_sim.sh`** runs ART on each haplotype at half the target coverage,
   concatenates, aligns with BWA-MEM to the **full reference**, sorts/indexes,
   and stamps a tumor read group. It simulates once at `MAXCOV` and produces the
   other `COVERAGES` by `samtools view -s` downsampling (fast). Every step is
   idempotent (skips finished outputs).

3. **`mix_tumor_normal.sh`** builds an impure tumor: `purity·tumor +
   (1−purity)·normal` at a target depth (estimating source depths from
   `idxstats` if not given). A het somatic SV at VAF 0.5 in pure tumor lands at
   ~`0.5·purity` in the mix. The unmixed normal stays the matched normal.

## SV types & sizes (defaults)
| type | count | size | model |
|---|---:|---|---|
| small DEL | 800 | 1–50 bp | deletion |
| small INS | 800 | 1–50 bp | novel-sequence insertion |
| DEL | 600 | 50 bp–3 Mb (log-unif) | deletion |
| DUP | 400 | 200 bp–1 Mb | tandem duplication |
| INV | 400 | 200 bp–1 Mb | inversion |
| TRA | 200 | inter-chromosomal | balanced reciprocal translocation |

≈3,200 SVs by default; scale with `--n-*`. All sizes/counts are CLI-tunable
(`python3 simulate_sv_genome.py --help`).

## Tumor impurity series
```bash
for p in 1.0 0.8 0.6 0.4 0.2; do
  ./mix_tumor_normal.sh -t sim_out/tumor.30x.bam -n $NORMAL_BAM \
      -p $p -d 30 -o sim_out/tumor.p${p}.30x.bam --tumor-cov 30
done
```

## One-shot panel builder (`build_sim_panel.sh`)
End-to-end driver that writes a **self-contained, reproducible panel** to a new
timestamped+seeded subfolder (default `/Volumes/wala24T/sim/<tag>_<date>_seedN/`)
containing the indexed impure tumor BAMs, all ground-truth files, and a README
recording exactly how it was made (incl. the seed and the regenerate command).
It generates many SVs (default chr22), simulates the pure tumor (ART→BWA), once
downsamples a **real whole-genome normal** (a 2nd normal from the same person)
to a cached contamination pool, then admixes it at each coverage×purity:
```bash
# preview the plan (no work):
DRYRUN=1 ./build_sim_panel.sh
# build it (uses up to THREADS=10):
./build_sim_panel.sh
# fast end-to-end validation with tiny inputs:
SMOKE=1 OUTROOT=/tmp/p CONTAM_NORMAL_BAM=some_small.bam ./build_sim_panel.sh
```
Key knobs (env): `COVERAGES="5 10 15 20"`, `PURITIES="0.5"`, `NORMAL_DOWNSAMPLE=0.5`,
`CONTAM_NORMAL_BAM=…`, `MATCHED_NORMAL_BAM=…` (the OTHER normal, for the svaba
T/N pair; README-only), `SEED`, `CHROMS`, `GEN_ARGS`, `OUTROOT`. Coverage is
matched on the SV chromosome. Pair the resulting tumor BAMs with your *other*
normal and score in `docs/sim_benchmark.html`.

## Run svaba over a whole panel (`run_svaba_on_panel.sh`)
Runs svaba T/N on **every** `tumor.*x*.bam` in a panel folder, then postprocess
+ benchmark, using the preferred invocation (`--blacklist`, `-p 20`, `-a`,
`--annotation rmsk/segdups`, `time`):
```bash
./run_svaba_on_panel.sh /Volumes/wala24T/sim/wgs_..._seed42   # or omit -> newest wgs_*
DRYRUN=1 ./run_svaba_on_panel.sh    # print the commands only
```
Outputs land in `<panel>/svaba_runs/<tag>/` and a recall/precision
`benchmark_summary.tsv` is written. **`-n` (matched normal) is the SECOND WGS
normal from the SAME individual** as the tumor's contamination (`NORMAL=` env;
default `~/Downloads/svaba_compare/blood.recal.bam`). Same person ⇒ shared
germline (germline SVs aren't called somatic); different reads from the
contamination ⇒ no identical-read artifact; only the simulated SVs are somatic.

## Coordinate conventions
All internal coordinates are 0-based half-open. Truth breakends (`truth.bedpe`):
- **DEL** affected `[a,b)` → breakends `(a,b)`; **INS** at `p` → `(p,p)`;
- **DUP** tandem `[a,b)` → `(a,b)`; **INV** `[a,b)` → `(a,b)`;
- **TRA** → `(chrA,posA)` and `(chrB,posB)`.

`benchmark.py` matches with a position **tolerance** (default 100 bp,
orientation-agnostic), so exact convention differences between callers don't
matter.

## Scaling & caveats
- Full hg38 at high coverage is large: ~10M reads/×-coverage, BAMs of tens of
  GB at 30×. Restrict `CHROMS` and/or lower `MAXCOV` to size the job. Lower
  coverages are free (downsampled).
- ART produces no PCR duplicates and a clean error model; this isolates SV
  *detectability*, it is not a model of library/artifact noise. Mixing with a
  **real** blood normal (default) reintroduces realistic germline background and
  makes the simulated SVs unambiguously somatic.
- Donor reads carry no germline SNPs unless `--germline-frac > 0` (which also
  writes a matched `normal.hap*.fa`). For a fully-controlled T/N pair use
  `--germline-frac` plus `SIM_NORMAL=1`.
- Microhomology/inserted bases at junctions are not yet modeled (blunt
  breakpoints; INS uses random novel sequence) — a planned knob.
