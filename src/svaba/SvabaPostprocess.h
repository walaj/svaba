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
// Scope (deliberately small, see conversation with user):
//   - Sort  → shell out to `samtools sort` (already highly tuned; htslib
//             doesn't expose its sort as a library call).
//   - Dedup → native streaming, (qname, flag) keyed, per-locus reset.
//   - Per-suffix parallelism via std::thread; suffixes run independently.
//
// Explicitly NOT in scope here (still in sort_output.sh):
//   - thread-BAM merging (pre-svaba-postprocess step)
//   - BAI indexing (post-step)
//   - SPLIT_BY_SOURCE qname-prefix demultiplexing
//
// CLI:  svaba postprocess -i <ID> [-t THREADS] [-m MEM] [-v V]
//                         [--sort-only | --dedup-only]
void runPostprocess(int argc, char** argv);
