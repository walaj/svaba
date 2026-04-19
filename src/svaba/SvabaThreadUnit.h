#pragma once

#include <sstream>
#include <map>
#include <memory>
#include <unordered_map>
#include <vector>
#include <mutex>

#include "SvabaBamWalker.h"
#include "SvabaLogger.h"
#include "SvabaOptions.h"
#include "AlignedContig.h"
#include "DiscordantCluster.h"
#include "BreakPoint.h"
#include "SvabaUtils.h"

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

  // SvABA2.0: map from read UniqueName -> comma-separated list of
  // contig cnames this read supports as ALT. Populated at BP
  // finalization so `SvabaOutputWriter::writeUnit` can stamp the
  // `bi:Z:<cname>` BAM aux tag on records in `all_corrected_reads`
  // (which are *newly aligned* BamRecords produced by bwa, so they
  // don't share identity with the original svabaRead pointers and
  // couldn't have been tagged via pointer). The same map is used to
  // (re)tag weird and discordant records at write time for
  // consistency. cname is a deterministic identifier coming from the
  // assembly window and is stable across runs (unlike a per-thread
  // minted counter), so `samtools view | grep bi:Z:<cname>` pulls all
  // ALT support for a given contig's variants.
  std::unordered_map<std::string, std::string> alt_cnames_by_name;

  // SvABA2.0: map from read UniqueName -> comma-separated list of
  // contig cnames this read aligned to (r2c), regardless of whether
  // the alignment ended up supporting a variant. Populated inside the
  // r2c alignment loop in SvabaRegionProcessor right alongside
  // `svabaRead::AddR2C(...)`, and consumed at write time by
  // `SvabaOutputWriter::writeUnit` to stamp the `bz:Z:<cname>` aux
  // tag. This is strictly a superset of `alt_cnames_by_name` — any
  // read with `bi:Z` will also have `bz:Z` (not necessarily the same
  // cname list, since a read can align to a contig yet not support a
  // variant on it).
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
