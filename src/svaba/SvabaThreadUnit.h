#pragma once

#include <sstream>
#include <map>
#include <memory>
#include <unordered_map>
#include <vector>
#include <mutex>
#include <cstdio>   // std::snprintf in next_bp_id()
#include <string>

#include "SvabaBamWalker.h"
#include "SvabaLogger.h"
#include "SvabaOptions.h"
#include "AlignedContig.h"
#include "DiscordantCluster.h"
#include "BreakPoint.h"
#include "SvabaUtils.h"
#include "gzstream.h"  // ogzstream for the per-thread r2c.txt.gz

#include "SeqLib/RefGenome.h"
#include "SvabaSharedConfig.h"
#include "SeqLib/BWAAligner.h"

using SeqLib::BamRecordPtrVector;

class svabaBamWalker;
namespace SeqLib {
  class BamWriter;
}
using WalkerMap = std::map<std::string, std::shared_ptr<svabaBamWalker>>;
using WriterMap = std::map<std::string, std::shared_ptr<SeqLib::BamWriter>>;

class svabaThreadUnit {
  
public:
  
  //svabaThreadUnit() = default;
  ~svabaThreadUnit();
  
  svabaThreadUnit(SvabaSharedConfig& sc_,
		  int thread);

  void flush();

  // local version of aligner class, but will hold shared memory index
  std::shared_ptr<SeqLib::BWAAligner> bwa_aligner; //(sc.bwa_idx);

  size_t processed_count = 0;
  size_t total_count = 0; // total to process
  size_t processed_since_memory_dump = 0;

  // results
  std::vector<AlignedContig>                 master_alc;
  BamRecordPtrVector                         master_contigs;
  BreakPointPtrVector                        m_bps;
  DiscordantClusterMap                       m_disc;
  //size_t                                     m_bamreads_count = 0;
  //size_t                                     m_disc_reads     = 0;
  //SeqLib::GRC                                badd; // bad regions
  int                                        threadId;

  // SvABA2.0 (v3): map from read UniqueName -> comma-separated list
  // of bp_ids this read supports as ALT. Populated at BP
  // finalization so `SvabaOutputWriter::writeUnit` can stamp the
  // `bi:Z:<bp_id,...>` BAM aux tag on records in
  // `all_corrected_reads` (which are *newly aligned* BamRecords
  // produced by bwa, so they don't share identity with the original
  // svabaRead pointers and couldn't have been tagged via pointer).
  // The same map is used to (re)tag weird and discordant records at
  // write time for consistency.
  //
  // Prior to v3 this carried cnames (deterministic contig-window
  // identifiers). Switched to bp_ids so the tag matches the
  // granularity of the r2c.txt.gz `split_bps` / `disc_bps` columns:
  // a single contig can carry several BPs (global + multi + indel),
  // and a read may support one, some, or all of them. Keying off the
  // BP, not the contig, removes the "which variant on this contig
  // does this read support?" ambiguity and makes
  // `samtools view | grep bi:Z:bp00100000042` return exactly the
  // ALT-supporting reads for that specific variant row in
  // bps.txt.gz. cname-level grouping is still available on the
  // companion `bz:Z` tag (below).
  std::unordered_map<std::string, std::string> alt_bp_ids_by_name;

  // SvABA2.0: map from read UniqueName -> comma-separated list of
  // contig cnames this read aligned to (r2c), regardless of whether
  // the alignment ended up supporting a variant. Populated inside the
  // r2c alignment loop in SvabaRegionProcessor right alongside
  // `svabaRead::AddR2C(...)`, and consumed at write time by
  // `SvabaOutputWriter::writeUnit` to stamp the `bz:Z:<cname>` aux
  // tag. This stays as cnames because "which contigs did this read
  // align to" is a contig-level concept — unlike `bi:Z` which is the
  // per-variant ALT-supporter attribution. The set of reads with a
  // `bz:Z` tag is a superset of the set with `bi:Z` (every
  // ALT-supporter aligned to the contig; not every r2c alignment
  // yields ALT support).
  std::unordered_map<std::string, std::string> all_cnames_by_name;

  // very verbose outpout
  svabaReadPtrVector                            all_weird_reads;
  BamRecordPtrVector                            all_corrected_reads;    

  // time and temp log dump until goes to file
  svabaUtils::svabaTimer                      st;
  std::stringstream                           ss;
  
  // store the BAM .bai indicies for for this thread
  WalkerMap                            walkers;
  WriterMap                            writers;

  // SvABA2.0: per-thread counter for generating unique BreakPoint IDs
  // without cross-thread coordination. Format emitted by next_bp_id()
  // is "bpTTTNNNNNNNN" (3-digit thread + 8-digit counter). Up to
  // 1000 threads and ~1e8 BPs per thread, far beyond any realistic
  // svaba run. If either ever overflows, the resulting ID is still a
  // unique string — just one digit wider — so downstream tools that
  // treat it opaquely are unaffected.
  size_t                               bp_id_counter = 0;
  std::string next_bp_id() {
    char buf[32];
    // %03d pads thread ID; %08zu pads the counter. post-increment so
    // the first BP for thread 1 is bp00100000000 (the worker pool
    // numbers threads 1..N — see threadpool.h).
    std::snprintf(buf, sizeof(buf), "bp%03d%08zu",
                  threadId, bp_id_counter++);
    return std::string(buf);
  }

  // SvABA2.0: per-thread r2c TSV stream. Mirrors the per-thread BAM
  // writers above — when --dump-reads is set (opts.dump_alignments),
  // each thread writes its variant-bearing contigs + r2c-aligned reads
  // into ${ID}.thread${N}.r2c.txt.gz directly, no shared handle, no
  // mutex. Postprocess merges the per-thread files via
  //   cat ${ID}.thread*.r2c.txt.gz > ${ID}.r2c.txt.gz
  // which produces a valid gzip (gzip is concatenation-safe per
  // RFC 1952). The first worker (threadId == 1; the pool numbers
  // workers 1..N) writes the column-header line once at open time
  // so the merged file has exactly one header at the top.
  //
  // Held via unique_ptr because ogzstream inherits from std::ios (which
  // has deleted copy/move). Embedding it by value would implicitly
  // delete svabaThreadUnit's move constructor and produce
  // -Wdefaulted-function-deleted. unique_ptr is trivially movable, so
  // the `svabaThreadUnit(svabaThreadUnit&&) noexcept = default;`
  // below stays valid. Same reason the BAM writers live in
  // shared_ptr<SeqLib::BamWriter> inside the `writers` map above.
  std::unique_ptr<ogzstream>           r2c_out_;
  
  // non-copyable, movable
  svabaThreadUnit(const svabaThreadUnit&) = delete;
  svabaThreadUnit& operator=(const svabaThreadUnit&) = delete;
  svabaThreadUnit(svabaThreadUnit&&) noexcept = default;
  svabaThreadUnit& operator=(svabaThreadUnit&&) noexcept = delete;

  void clear();

  bool MemoryLimit(size_t readLimit, size_t contLimit) const;

  // store the faidx index for this thread
  std::unique_ptr<SeqLib::RefGenome>   ref_genome;

  SvabaSharedConfig& sc;

private:

};
