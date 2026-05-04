#pragma once

// Post-process svaba outputs: merge per-thread files, coordinate-sort BAMs,
// streaming-dedup, sort/dedup/filter bps.txt.gz, filter r2c.txt.gz, and
// optionally split BAMs by QNAME source prefix.
//
// Fully replaces svaba_postprocess.sh. Everything runs natively in C++
// except `samtools sort` (shelled out — htslib doesn't expose its sort
// as a library call).
//
// Phases (all auto-detect inputs; no-op if absent):
//
//   Phase 0: Merge per-thread outputs.
//     - BAMs: ${ID}.thread*.${suffix}.bam → ${ID}.${suffix}.bam via native
//       htslib record streaming. No samtools dependency.
//     - r2c: ${ID}.thread*.r2c.txt.gz → ${ID}.r2c.txt.gz via binary gzip
//       concatenation (RFC 1952 concat-safe).
//
//   Phase 1: Parallel coordinate-sort of BAMs (samtools sort).
//     One thread-pool slice per suffix, all suffixes concurrently.
//
//   Phase 2: Serial dedup + reheader + index per BAM.
//     - Streaming dedup: (qname, flag) keyed per-locus, with bi:Z/bz:Z
//       comma-token union (preserves contig support evidence from
//       overlapping assembly windows). Full BGZF thread pool per BAM.
//     - @PG stamp: chained onto existing PG chain.
//     - .bai index: native sam_index_build.
//     Idempotent: auto-skips when header already has svaba_postprocess PG.
//
//   Phase 3: Sort + dedup + filter bps.txt.gz.
//     In-memory: contiguous slab + lightweight sort index. Canonical AB==BA
//     dedup with lowest-somlod winner selection (conservative against false
//     somatic calls). Outputs: .sorted, .sorted.dedup, .sorted.dedup.pass,
//     .sorted.dedup.pass.somatic.
//
//   Phase 4: Filter r2c.txt.gz to PASS / PASS+somatic contigs.
//     Single streaming pass writing two gzipped outputs, keyed on the PASS
//     cname sets collected during Phase 3.
//
//   Phase 5: Split-by-source (optional, --split).
//     For each dedup-eligible suffix (weird/corrected/discordant), reads
//     the BAM and routes records by the first 4 chars of QNAME into
//     per-prefix BAMs. Each prefix BAM is then sorted + indexed.
//     Output: ${ID}.${suffix}.${PREFIX}.bam(.bai).
//
// CLI:  svaba postprocess -i <ID> [-t THREADS] [-m MEM] [-v V]
//                         [--sort-only | --dedup-only] [--split]
void runPostprocess(int argc, char** argv);
