//run_svaba.cpp

#include <thread>
#include <mutex>
#include "threadpool.h"
#include <memory>

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
#include "svabaLogger.h"
#include "KmerFilter.h"
#include "vcf.h"
#include "DBSnpFilter.h"
#include "svabaUtils.h"
#include "SeqLib/BFC.h"
#include "svaba_params.h"

#include "svabaFileLoader.h"
#include "svabaOutputWriter.h"
#include "svabaAssemblerEngine.h"

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

// time
static struct timespec start;

// forward declaration, need non-static and not in
// anonymous namespace to keep external linkage
void  runsvaba       (int argc, char** argv);

// --- one small helper functor ---
struct SvabaWorkItem {
  SeqLib::GenomicRegion    region;
  int                      id;
  SvabaRegionProcessor&    processor;

  SvabaWorkItem(const SeqLib::GenomicRegion& r,
                int n,
                SvabaRegionProcessor& proc)
    : region(r), id(n), processor(proc)
  {}

  bool operator()(svabaThreadUnit& unit, size_t threadId) const {
    return processor.process(region, unit, threadId);
  }
};

void sendThreads(const SeqLib::GRC& regionsToRun,
	 SvabaSharedConfig& sc) {
  
  // single shared processor:
  SvabaRegionProcessor proc(sc); 
  
  // our threadpool of work items:
  ThreadPool<SvabaWorkItem> pool(sc);
  
  int count=0;
  // submit one job per region  
  for (auto& r : regionsToRun)
    pool.submit(std::make_unique<SvabaWorkItem>(r, ++count, proc));
  
  // if no intervals, submit a single wholegenome job  
  if (regionsToRun.IsEmpty())
    pool.submit(std::make_unique<SvabaWorkItem>(SeqLib::GenomicRegion(), ++count, proc));
  pool.shutdown();
}
  
void makeVCFs(SvabaSharedConfig& sc) {

  // make the VCF file
  string file = sc.opts.analysisId + ".bps.txt.gz";  
  sc.logger.log(true, true, "...loading the bps file ", file," for conversion to VCF");

  // make the header
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
  VCFFile snowvcf(file, sc.opts.analysisId, sc.header, header, !false, sc.opts.verbose > 0);
  
  string basename = sc.opts.analysisId + ".svaba.unfiltered.";
  snowvcf.include_nonpass = true;
  sc.logger.log(true, true, "...writing unfiltered VCFs");
  snowvcf.writeIndels(basename, false, !case_control_run);
  snowvcf.writeSVs(basename, false,    !case_control_run);
  
  sc.logger.log(true, true, "...writing filtered VCFs");
  basename = sc.opts.analysisId + ".svaba.";
  snowvcf.include_nonpass = false;
  snowvcf.writeIndels(basename, false, !case_control_run);
  snowvcf.writeSVs(basename, false,    !case_control_run);

} //anonymous namespace

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
  
  // check that the two headers are equivalant
  // open the main bam to get header info
  SeqLib::BamReader first_tumor_bam_reader;  
  if (!first_tumor_bam_reader.Open(opts.main_bam)) {
    cerr << "ERROR: Cannot open main bam file: " << opts.main_bam << endl;
    exit(EXIT_FAILURE);
  }
  SeqLib::BamHeader first_header = first_tumor_bam_reader.Header();
  
  svabaUtils::checkHeaderCompatibility(first_header, bwa_header, logger);
  
  // open file loader
  SvabaFileLoader loader(sc); 
  
  // open the blacklist, load into sc.blacklist
  loader.loadBedRegions(sc.opts.blacklistFile, sc.blacklist);
  
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
  SeqLib::GRC regions_torun;
  loader.countJobs(regions_torun);  

  // --- learn the insert-sizes ---
  logger.log(true, true,"... learning insert size distribution across all BAMs; this may take a while");
  
  // learn from the BAM files
  for (const auto& b : opts.bams) {
    auto [it, inserted] = sc.bamStats.emplace(b.first, LearnBamParams(sc, b.second));
    it->second.learnParams();
  }
  logger.log(true, true,"... done learning insert size distribution");

  // --- report what we learned ---
  int globalReadLen = 0;
  int globalMaxMapQ = 0;
  for (const auto& ll : sc.bamStats) {
    globalReadLen = std::max(globalReadLen, ll.second.readlen_max);
    globalMaxMapQ = std::max(globalMaxMapQ, ll.second.mapq_max);
  }
  sc.readlen = globalReadLen;
  
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
  logger.log(false, true,
    "... calculated seedLength = ", seedLength,
    " (error_rate=", opts.sgaErrorRate,
    ", readlen=", globalReadLen, ")"
  );

  // --- build per-RG discordant size rules ---
  std::stringstream ssRules;
  std::unordered_set<std::string> seenRG;
  std::unordered_map<std::string,int> minIsizeForDisc;
  for (auto const& [sample, pm] : sc.bamStats) {
    for (auto const& [rg, bp] : pm.bam_read_groups) {
      double cutoffd = bp.isize_mean + bp.sd_isize * sc.opts.sdDiscCutoff * 3.0;
      int cutoff = int(std::floor(cutoffd));
      ssRules << R"({"isize":[)" << cutoff << R"(,0],"rg":")" << rg << R"("})" << ",";
    }
  }
  if (!ssRules.str().empty())
    ssRules.seekp(-1, ssRules.cur);  // drop trailing comma
  sc.logger.log(true, true, "[INFO] Learned discordant size cutoffs by RG:");
  
  // plug into your JSON rules template
  sc.opts.rulesJson = svabaUtils::myreplace(
    sc.opts.rulesJson, "FRRULES", ssRules.str().empty()
      ? "{}"
      : ssRules.str()
  );

  // similarly replace any READLENLIM token
  if (sc.opts.rulesJson.find("READLENLIM") != std::string::npos) {
    int clipLen = std::max(30, int(0.3 * globalReadLen));
    sc.opts.rulesJson = svabaUtils::myreplace(
      sc.opts.rulesJson, "READLENLIM", std::to_string(clipLen)
    );
  }

  // set the ReadFilterCollection to be applied to each region
  logger.log(opts.verbose > 1, true, opts.rulesJson);
  sc.mr = ReadFilterCollection(opts.rulesJson, sc.header);
  logger.log(opts.verbose > 1, true, sc.mr);

  // put args into string for VCF later
  sc.args += "(v" + string(SVABA_VERSION) + ") ";
  for (int i = 0; i < argc; ++i)
    sc.args += string(argv[i]) + " ";

  // start the timer
#ifndef __APPLE__
  clock_gettime(CLOCK_MONOTONIC, &start);
#endif

  // send the jobs to the queue
  logger.log(true, true, "Starting detection pipeline");
  sendThreads(regions_torun, sc); 

  // make the VCF file
  makeVCFs(sc);
  
#ifndef __APPLE__
  cerr << SeqLib::displayRuntime(start) << endl;
#endif
}


