//run_svaba.cpp

#include <thread>
#include <mutex>
#include "threadpool.h"
#include <memory>

#include <algorithm>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <map>
#include <vector>
#include <cassert>

#include "SeqLib/GenomicRegion.h"
#include "SeqLib/UnalignedSequence.h"
#include "SeqLib/ReadFilter.h"
#include "SeqLib/BamHeader.h"
#include "SeqLib/BWAAligner.h"
#include "SeqLib/ReadFilter.h"

#include "SvabaRegionProcessor.h"
#include "SvabaLogger.h"
#include "KmerFilter.h"
#include "vcf.h"
#include "DBSnpFilter.h"
#include "SvabaUtils.h"
#include "SeqLib/BFC.h"

#include "SvabaFileLoader.h"
#include "SvabaOutputWriter.h"
#include "SvabaAssemblerEngine.h"
#include "SvabaAssemblerConfig.h"

using SeqLib::BWAIndex;
using SeqLib::BWAAligner;
using SeqLib::BamHeader;
using SeqLib::Filter::ReadFilterCollection;
using SeqLib::GenomicRegion;
using SeqLib::GRC;
using SeqLib::BWAIndexPtr;

using std::string;
using std::cerr;
using std::endl;
using std::make_shared;

// 1 is normal logging
// 2 is heavy
constexpr inline int log_level = 1;

// forward declaration, need non-static and not in
// anonymous namespace to keep external linkage
void  runsvaba       (int argc, char** argv);

struct SvabaBatchWorkItem {
  // batch of region + unique region ID pairs
  std::vector<std::pair<SeqLib::GenomicRegion,int>> items;
  SvabaRegionProcessor& processor;

  // constructor takes a reference to the shared processor
  SvabaBatchWorkItem(const std::vector<std::pair<SeqLib::GenomicRegion,int>>& items_,
                     SvabaRegionProcessor& proc)
    : items(items_), processor(proc) {}

  // call operator matching SvabaWorkItem signature
  bool operator()(svabaThreadUnit& unit, size_t threadId) const {
    for (auto& [region, id] : items) {
      processor.process(region, unit, threadId);
    }
    return true;
  }
};

// Modified sendThreads() to submit batches of regions instead of one-by-one jobs
void sendThreads(const SeqLib::GRC& regionsToRun,
                 SvabaSharedConfig& sc) {
  
  // shared processor instance
  SvabaRegionProcessor proc(sc);

  // thread pool now handles SvabaBatchWorkItem jobs
  ThreadPool<SvabaBatchWorkItem> pool(sc);
  
  size_t batch_size = 1; 
  if (sc.opts.numThreads > 0) {
    size_t per_thread = sc.total_regions_to_process / sc.opts.numThreads;
    if (per_thread < batch_size) {
      batch_size = std::max<size_t>(1, per_thread); // prevent batch_size = 0
    }
  }
  
  std::vector<std::pair<SeqLib::GenomicRegion,int>> batch;
  batch.reserve(batch_size);

  int count = 0;
  for (auto& region : regionsToRun) {
    batch.emplace_back(region, ++count);
    if (batch.size() == batch_size) {
      pool.submit(std::make_unique<SvabaBatchWorkItem>(batch, proc));
      batch.clear();
    }
  }
  // submit any leftover regions as one final batch
  if (!batch.empty()) {
    pool.submit(std::make_unique<SvabaBatchWorkItem>(batch, proc));
  }

  pool.shutdown();
}
  
void makeVCFs(SvabaSharedConfig& sc) {

//   // make the VCF file
   string file = sc.opts.analysisId + ".bps.txt.gz";  
   sc.logger.log(true, true, "...loading the bps file ", file," for conversion to VCF");

//   // make the header
   VCFHeader header;
   header.filedate = svabaUtils::fileDateString();
   header.source = sc.args;
   header.reference = sc.opts.refGenome;
  
   for (int i = 0; i < sc.header.NumSequences(); ++i)
     header.addContigField(sc.header.IDtoName(i),sc.header.GetSequenceLength(i));

   for (auto& b : sc.opts.bams) {
     string fname = b.second; //bpf.filename();
     header.addSampleField(fname);
     header.colnames += "\t" + fname; 
   }

   // check if it has a matched control. If so, output "somatic / germline" vcfs
   bool case_control_run = false;
   for (auto& b : sc.opts.bams)
     if (b.first.at(0) == 'n')
       case_control_run = true;

   // primary VCFs
   sc.logger.log(true, true, "...making the primary VCFs (unfiltered and filtered) from file ",file);
   VCFFile snowvcf(file, sc.opts.analysisId, sc, header, false, sc.opts.verbose > 0);
  
   string basename = sc.opts.analysisId + ".svaba.unfiltered.";
   snowvcf.include_nonpass = true;
   sc.logger.log(true, true, "...writing unfiltered VCFs");
   snowvcf.writeIndels(basename, false, !case_control_run, sc.header);
   snowvcf.writeSVs(basename, false,    !case_control_run, sc.header);
  
   sc.logger.log(true, true, "...writing filtered VCFs");
   basename = sc.opts.analysisId + ".svaba.";
   snowvcf.include_nonpass = false;
   snowvcf.writeIndels(basename, false, !case_control_run, sc.header);
   snowvcf.writeSVs(basename, false,    !case_control_run, sc.header);

}

void runsvaba(int argc, char** argv) {

  // parse command line
  SvabaOptions opts;
  try {
    opts = SvabaOptions::parse(argc, argv);
    if (opts.help) {
      SvabaOptions::printUsage();
      return;
    }
  }
  catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return;
  }

  // instantiate logger, options, writer here.
  // its just easier then to store refernce in SvabaSharedConfig
  // to avoid header dependencies
  SvabaLogger logger;
  logger.init(opts.analysisId + ".log");
  //logger.welcome(opts); // initial message

  // open the human reference
  logger.log(true, true, "...loading the human reference sequence for BWA");  
  BWAIndexPtr bwa_idx = make_shared<SeqLib::BWAIndex>();
  bwa_idx->LoadIndex(opts.refGenome);
  
  // get the dictionary from reference
  BamHeader bwa_header = bwa_idx->HeaderFromIndex();
  
  // open the writer
  SvabaOutputWriter writer(logger, opts);
  writer.init(opts.analysisId, bwa_header);

  // shared information to be passed around thread units
  SvabaSharedConfig sc(logger, opts, writer);
  sc.bwa_idx = bwa_idx;
  sc.header = bwa_header;
  
  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &sc.start);

  // report which local-assembly engine was compiled in
  std::cerr << "...local-assembly engine: " << svaba::kAssemblerName
#if SVABA_ASSEMBLER_FERMI
            << " (FERMI assembly enabled)"
#else
            << " (SGA assembly enabled)"
#endif
            << std::endl;

  // check that the two headers are equivalant
  // open the main bam to get header info
  SeqLib::BamReader first_tumor_bam_reader;  
  if (!first_tumor_bam_reader.Open(opts.main_bam)) {
    cerr << "ERROR: Cannot open main bam file: " << opts.main_bam << endl;
    exit(EXIT_FAILURE);
  }
  SeqLib::BamHeader first_header = first_tumor_bam_reader.Header();
  
  svabaUtils::checkHeaderCompatibility(first_header, bwa_header, logger);

  // Check input BAM @PG tags for BWA. When --always-realign-corrected is
  // NOT set, svaba reuses the input BAM's native CIGAR/NM for unchanged
  // reads in the r2c-vs-native gate. That comparison is only valid if the
  // input was aligned with BWA (same scoring model as svaba's internal
  // aligner). If the PG chain doesn't mention "bwa", warn loudly.
  if (!opts.alwaysRealignCorrected) {
    std::string hdr_text = first_header.AsString();
    // case-insensitive search: look for "bwa" anywhere in @PG lines
    bool found_bwa = false;
    std::istringstream hdr_stream(hdr_text);
    std::string line;
    while (std::getline(hdr_stream, line)) {
      if (line.substr(0, 3) != "@PG") continue;
      // case-insensitive: lowercase the line for matching
      std::string lower_line = line;
      std::transform(lower_line.begin(), lower_line.end(),
                     lower_line.begin(), ::tolower);
      if (lower_line.find("bwa") != std::string::npos) {
        found_bwa = true;
        break;
      }
    }
    if (!found_bwa) {
      logger.log(true, true,
        "");
      logger.log(true, true,
        "************************************************************");
      logger.log(true, true,
        "* WARNING: No BWA @PG tag found in input BAM header.      *");
      logger.log(true, true,
        "* svaba reuses the input BAM's CIGAR/NM for reads that    *");
      logger.log(true, true,
        "* BFC did not modify. This is only valid when the input   *");
      logger.log(true, true,
        "* was aligned with BWA-MEM (same scoring model as svaba's *");
      logger.log(true, true,
        "* internal aligner).                                      *");
      logger.log(true, true,
        "*                                                         *");
      logger.log(true, true,
        "* If this BAM was NOT aligned with BWA, re-run with:      *");
      logger.log(true, true,
        "*   --always-realign-corrected                            *");
      logger.log(true, true,
        "* to force re-alignment of every read to the reference.   *");
      logger.log(true, true,
        "************************************************************");
      logger.log(true, true,
        "");
    }
  }

  // open file loader
  SvabaFileLoader loader(sc); 
  
  // open the blacklist, load into sc.blacklist
  for (const auto& path : sc.opts.blacklistFile) {
    SeqLib::GRC this_blacklist;
    loader.loadBedRegions(path, this_blacklist);
    sc.blacklist.Concat(this_blacklist);
  }
  sc.blacklist.MergeOverlappingIntervals();
  sc.blacklist.CreateTreeMap();
  
  logger.log(true, true, "...loaded ", SeqLib::AddCommas(sc.blacklist.size()),
	     " blacklist regions across ", SeqLib::AddCommas(sc.opts.blacklistFile.size()),
	     " files");
  
  // open the germline sv database
  loader.loadBedRegions(sc.opts.germlineSvFile, sc.germline_svs);
  
  // open the DBSnpFilter
  if (sc.opts.dbsnpVcf.length()) {
    sc.dbsnp_filter = std::make_shared<DBSnpFilter>(opts.dbsnpVcf, bwa_header, logger);
  }
  
  // needed for aligned contig
  for (auto& b : opts.bams)
    sc.prefixes.insert(b.first);

  // parse the region file, count number of jobs
  SeqLib::GRC regionsToRun;
  loader.countJobs(regionsToRun);

  // SvABA2.0: drop any queued region that is 100% covered by the blacklist.
  // svaba already filters reads by blacklist *inside* each region (see
  // SvabaRegionProcessor.cpp:74), so a chunk that sits entirely in a
  // blacklisted contig (e.g. the chrUn / *_decoy / *_alt contigs in the
  // nonstd_chr blacklist) would otherwise still get its full region
  // pipeline run — open the BAM walker, learn params, attempt local
  // assembly on a stream of reads that all get dropped — just to produce
  // zero callable bases. Net behavior is unchanged; this is pure
  // optimization of the queue.
  //
  // The blacklist has had MergeOverlappingIntervals() + CreateTreeMap()
  // called already (a few lines above), so FindOverlapWidth is an O(k log n)
  // tree walk over the blacklist intervals that touch each region. Runs
  // once, pre-queue, so it's off the critical path.
  if (!sc.blacklist.empty() && regionsToRun.size() > 0) {
    SeqLib::GRC kept;
    size_t dropped_n  = 0;
    size_t dropped_bp = 0;
    for (const auto& r : regionsToRun) {
      const int w = r.Width();
      if (w > 0 &&
          static_cast<int>(sc.blacklist.FindOverlapWidth(r, /*ignore_strand=*/true)) >= w) {
        ++dropped_n;
        dropped_bp += static_cast<size_t>(w);
        continue;
      }
      kept.add(r);
    }
    if (dropped_n > 0) {
      logger.log(true, true,
        "...pruned ", SeqLib::AddCommas(dropped_n), " of ",
        SeqLib::AddCommas(regionsToRun.size()),
        " queued regions (", SeqLib::AddCommas(dropped_bp),
        " bp) that were 100% covered by the blacklist");
      regionsToRun = std::move(kept);
      // No runtime.txt line will be emitted for a pruned region; the
      // accounting in sc.total_regions_to_process below reflects what
      // actually gets queued.
    }
  }

  sc.total_regions_to_process = regionsToRun.size();
  if (sc.total_regions_to_process < opts.numThreads) {
    opts.numThreads = sc.total_regions_to_process;
  }

  // debug
  /*  regionsToRun.clear();
  for (int i = 0; i < 10; i++)
    regionsToRun.add(GenomicRegion(2,89000000,91500000));
  */
  
  // --- learn the insert-sizes ---
  logger.log(true, true,"...learning insert size distribution across all BAMs; this may take a while");
  
  // learn from the BAM files
  for (const auto& b : opts.bams) {
    auto [it, inserted] = sc.bamStats.emplace(b.first, LearnBamParams(sc, b.second));
    logger.log(true,true,"......learning BAM: ", b.second);
    it->second.learnParams();
  }

  // --- report what we learned ---
  for (const auto& bs : sc.bamStats) {
    logger.log(opts.verbose > 1, true, "=====", bs.first);
    for (const auto& brg : bs.second.bam_read_groups)
      logger.log(opts.verbose > 1, true, brg.first, " - ", brg.second);
  }

  int globalReadLen = 0;
  int globalMaxMapQ = 0;
  double globalInsertSize = 0;
  for (const auto& [_, learn_bam_param] : sc.bamStats) {
    globalReadLen = std::max(globalReadLen, learn_bam_param.readlen_max);
    globalMaxMapQ = std::max(globalMaxMapQ, learn_bam_param.mapq_max);
    globalInsertSize = std::max(globalInsertSize, learn_bam_param.isize_max);    
  }
  sc.readlen = globalReadLen;
  sc.insertsize = globalInsertSize;
  
  // --- set the SGA min overlap if user didn't ---
  if (opts.sgaMinOverlap == 0) {
    opts.sgaMinOverlap = std::max(30, int(0.6 * globalReadLen));
  }
  
  logger.log(true, true,
    "... found max read length = ", globalReadLen,
    "; SGA minOverlap = ", opts.sgaMinOverlap,
    "; max MAPQ = ", globalMaxMapQ
  );

  // --- compute seed parameters ---
  int seedLength, seedStride;
  {
    svabaAssemblerEngine tester(
      "test", opts.sgaErrorRate,
      opts.sgaMinOverlap,
      globalReadLen
    );
    tester.calculateSeedParameters(
      globalReadLen,
      opts.sgaMinOverlap, seedLength, seedStride
    );
  }
  logger.log(opts.verbose > 0, true, "...seedLength = ", seedLength,", readlen=", globalReadLen, ")" );

  // --- build per-RG discordant size rules ---
  for (auto const& [sample, pm] : sc.bamStats) {
    for (auto const& [rg, bp] : pm.bam_read_groups) {
      double cutoffd = bp.isize_mean + bp.sd_isize * sc.opts.sdDiscCutoff * 3.0;
      int cutoff = int(std::floor(cutoffd));
      opts.addFRRule(rg, cutoff);
    }
  }

  // set the ReadFilterCollection to be applied to each region
  logger.log(opts.verbose > 1, true, opts.rulesJson);
  sc.mr = ReadFilterCollection(opts.rulesJson, sc.header);
  logger.log(opts.verbose > 0, true, sc.mr);

  // put args into string for VCF later
  sc.args += "(v" + string(SVABA_VERSION) + ") ";
  for (int i = 0; i < argc; ++i)
    sc.args += string(argv[i]) + " ";

  if (regionsToRun.size()) {
    logger.log(true, true, "...running on ", SeqLib::AddCommas(regionsToRun.size()),
		" chunks"); 
  } else {
    logger.log(true, true, "Chunk was <= 0: Reading in whole-genome at once");
  }

  // send the jobs to the queue
  logger.log(true, true, "Starting detection pipeline");
  sendThreads(regionsToRun, sc); 

  // close the writer
  writer.close();
  
  // make the VCF file
  //makeVCFs(sc);
  
  cerr << SeqLib::displayRuntime(sc.start) << endl;
}


