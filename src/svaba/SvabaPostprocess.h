#pragma once

// Post-process svaba output BAMs: coordinate-sort via samtools, then run a
// native streaming dedup on the suffixes that can accumulate exact duplicate
// records from overlapping assembly windows (weird, corrected, discordant).
//
// Replaces the slow path of sort_output.sh, which went through
// `samtools view | awk | samtools view` — that pipeline decompresses every
// record to SAM text, maintains a large awk hash, and recompresses on the
// way out. This module keeps everything in native htslib / SeqLib records
// and resets the dedup set at every locus boundary, so memory is O(reads
// per locus) and throughput is close to what the disk+decoder can sustain.
//
// Scope:
//   - Sort      → shell out to `samtools sort` (already highly tuned; htslib
//                 doesn't expose its sort as a library call).
//   - Dedup     → native streaming, (qname, flag) keyed, per-locus reset.
//                 Stamps an @PG line onto the output header as a free
//                 side-effect of rewriting the file.
//   - Reheader  → for suffixes that don't go through dedup (e.g. `contigs`,
//                 or anything under --sort-only), shell out to
//                 `samtools reheader` to stamp the same @PG line. Streams the
//                 BGZF body as opaque blocks — cheap.
//   - Index     → native `sam_index_build` (htslib) producing a .bai
//                 alongside every final BAM.
//   - Per-suffix parallelism via std::thread; suffixes run independently.
//
// End state: ${ID}.${suffix}.bam is the single authoritative output per
// suffix — sorted, deduped (where applicable), @PG-stamped, and indexed.
// Intermediate files use a .postprocess.*.tmp.bam suffix and are cleaned up
// on both success (rename-over) and failure (best-effort unlink).
//
// Explicitly NOT in scope here (still in sort_output.sh):
//   - thread-BAM merging (pre-svaba-postprocess step)
//   - SPLIT_BY_SOURCE qname-prefix demultiplexing
//
// CLI:  svaba postprocess -i <ID> [-t THREADS] [-m MEM] [-v V]
//                         [--sort-only | --dedup-only]
void runPostprocess(int argc, char** argv);
