# CLAUDE.md — svaba working notes

This file captures conventions, file landmarks, and open investigations for the
svaba SV/indel caller project so future sessions can pick up quickly. Update it
as understanding changes.

## Project at a glance

svaba is a structural variant (SV) and indel caller that uses local assembly +
read realignment to call variants from short-read BAMs. The canonical use case
is tumor/normal somatic calling, but it also supports germline and multi-sample
modes.

Top-level layout:

- `src/svaba/` — main C++ sources. Entry point is `run_svaba.cpp`; assembly,
  realignment, breakpoint scoring, and VCF output live here.
- `src/SGA/` — String Graph Assembler sources (vendored).
- `src/svabautils/` — small utility lib shared across components.
- `SeqLib/` — vendored htslib/bwa wrapper used for BAM I/O and alignment.
- `bin/`, `build/` — build artifacts; don't edit by hand.
- `R/`, `viewer/`, `tracks/` — downstream analysis/visualization helpers.
- `tests/` — test fixtures.
- `sort_output.sh`, `extract_discordants.sh`, `memprof*.sh` — post-processing
  and profiling shell helpers that live at the repo root.
- `somlod_maxlod_analysis.html` — deep-dive writeup of the somatic log-odds
  scoring model, what goes wrong, and a menu of proposed fixes. See below.

## Statistical model — the files that matter

For anything related to how variants are scored (LOD, somatic vs germline,
error model), the two files to read first are:

- `src/svaba/svabaModels.cpp` / `.h` — self-contained statistical primitives.
  - `LogLikelihood(d, a, f, e_fwd, e_rev)` (lines ~11-63) is the per-sample
    two-state error model: `p_ref = (1-f)(1-e_fwd) + f*e_rev` and
    `p_alt = f(1-e_rev) + (1-f)*e_fwd`. Returns log10 likelihood of observing
    `a` alt reads out of `d`. This is the primitive every higher-level score
    is built from.
  - `SomaticLOD(...)` (lines ~70-83) — public wrapper; forwards to the
    split-error implementation.
  - `SomaticLOD_withSplitErrors(...)` (lines ~86-189) — the active somatic
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
    `al.LO` across samples. This is the "is this artifact or not" score and
    it does grow with additional supporting reads because it compares a
    variant hypothesis to a pure-error hypothesis with no germline branch.
  - Lines ~1072-1073: the current INDEL somatic gate only tests `somlod` —
    it does not use `maxlod` as a co-gate. This is one of the levers in the
    proposed fixes.

## The somlod / maxlod investigation (current open thread)

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
  positive. I verified this numerically while writing the analysis.
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

**Important correctness notes from the investigation (earned the hard way):**
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

## sort_output.sh — recent edit

`sort_output.sh` now has an optional source-splitting module at the bottom
gated by a `SPLIT_BY_SOURCE` flag (default `0`). When enabled it splits the
`discordant`, `weird`, and `corrected` BAMs by the first 4 chars of each
read's QNAME (the "source tag" like `t001`) into per-source BAMs named
`${ID}.${suffix}.${prefix}.bam`. The split is a single-pass `awk` route into
per-prefix SAM files in a tempdir, then reassembly via
`samtools view -bS | sort | index`. Invoke with
`SPLIT_BY_SOURCE=1 ./sort_output.sh ID` or flip the default at the top.

## Conventions

- C++ style matches the existing codebase (snake_case methods inside
  ClassName, 2-space indent, header/impl split). Don't introduce new
  formatting unless asked.
- Statistical code lives in `svabaModels.*`; breakpoint glue lives in
  `BreakPoint.*`. Keep statistical primitives in the models file, not
  inlined into BreakPoint.
- LL/LOD values in this codebase are always **log10**, not natural log.
- `aN, dN, aT, dT` = alt count / depth in normal and tumor; `f` = allele
  fraction; `e_fwd`/`e_rev` = forward/reverse error rates from the artifact
  model. These names are used consistently in the analysis HTML too.

## Useful jump points

- Somatic LOD calc: `src/svaba/svabaModels.cpp:86`
- Per-sample LO: `src/svaba/BreakPoint.cpp:1610`
- Somatic LOD entry: `src/svaba/BreakPoint.cpp:975`
- INDEL somatic gate: `src/svaba/BreakPoint.cpp:1072`
- Analysis writeup: `somlod_maxlod_analysis.html`
