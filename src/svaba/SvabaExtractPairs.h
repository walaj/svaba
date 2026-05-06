#pragma once

// Extract every read pair from a BAM where either mate's SEQ matches any of
// a small set of query sequences (or their reverse complements). Drop-in
// replacement for scripts/extract_pairs_by_seq.sh, but BAM-native end-to-end
// — no SAM text round-trip, no gratuitous final sort.
//
// Algorithm (two passes, neither leaves BAM space):
//
//   Pass 1 — QNAME collection.
//     Stream input via SeqLib::BamReader with a BGZF decompression thread
//     pool. For each record, walk the 4-bit-packed SEQ directly (no
//     std::string alloc) through an Aho-Corasick automaton built once over
//     {query, revcomp(query)} for every query. On any pattern hit, insert
//     QNAME into a global hash set.
//
//   Pass 2 — pair extraction.
//     Re-stream input. For every record whose QNAME is in the set, write to
//     the output BAM. Both mates plus any supplementary/secondary
//     alignments survive together because they all share QNAME. Output
//     header is the input header + an appended @PG line stamped
//     "svaba_extract_pairs".
//
//   Sort.
//     Skipped iff the input declares @HD SO:coordinate (the common case).
//     Output preserves input record order, so coord-sorted in → coord-sorted
//     out, no sort needed. For unsorted input we shell out to
//     `samtools sort` — same fallback as `svaba postprocess`. The legacy
//     shell script unconditionally sorted, which on a large already-sorted
//     BAM is a wasted full read+write of the output; this skips that.
//
//   Index.
//     sam_index_build(out, 0) → out.bai. No samtools index fork.
//
// Why this is fast vs. the shell script:
//   * No samtools-view → awk → samtools-view-b round-trip. The shell
//     pipeline decompresses BAM → ASCII SAM → re-encodes ASCII SAM → BAM,
//     bottlenecked on a single-threaded awk consumer of the SAM text. The
//     C++ path stays in BAM throughout and uses BGZF thread pools on both
//     reader and writer.
//   * No regex alternation on the SEQ field. Aho-Corasick gives one linear
//     scan of each read SEQ regardless of how many query sequences (and
//     their reverse complements) are in flight.
//   * No third pass for `samtools sort` when the input is already
//     coord-sorted — dropped on the basis of the @HD SO: tag.
//
// CLI:
//   svaba extract-pairs -i IN.bam -o OUT.bam (-s SEQ ... | -f FILE) [-t N]
//                                            [--no-rc] [-v V]
//
// See SvabaExtractPairs.cpp for the implementation.
void runExtractPairs(int argc, char** argv);
