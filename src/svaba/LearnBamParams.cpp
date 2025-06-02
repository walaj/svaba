// LearnBamParams.cpp

#include "LearnBamParams.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/BamReader.h"
#include "SvabaSharedConfig.h"
#include "svabaOptions.h"
#include "svabaLogger.h"

#include <stdexcept>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace {
  
  // helper to extract read group ids from the header
  std::vector<std::string> extractReadGroups(const SeqLib::BamHeader& hdr) {
    std::vector<std::string> groups;
    std::istringstream lines(hdr.AsString());
    std::string line;
    while (std::getline(lines, line)) {
      if (line.rfind("@RG", 0) != 0) continue;
      std::istringstream fields(line);
      std::string field;
      while (std::getline(fields, field, '\t')) {
	if (field.rfind("ID:", 0) == 0) {
	  groups.push_back(field.substr(3));
	  break;
	}
      }
    }
    return groups;
  }
  
} // anonymous namespace

LearnBamParams::LearnBamParams(SvabaSharedConfig& sc_,
			       const std::string& bamPath) 
  : sc(sc_),
    bam_(bamPath)
{

  // make the new BamReader
  reader_ = std::make_shared<SeqLib::BamReader>();
  
  // open the BAM for learning
  if (!reader_->Open(bam_)) {
    throw std::runtime_error("Could not open bam file " + bam_ + " in LearnBamParams");
  }
}

void LearnBamParams::learnParams() {

  // gather read group ids from the header
  auto hdr = reader_->Header();
  auto groups = extractReadGroups(hdr);

  // read-group, count
  std::unordered_map<std::string, size_t> rg_count;

  // track how many read groups are satisfied. Break when all are done
  size_t satisfied = 0;
  
  // set the BAM reader to the beginning of the BAM
  reader_->Reset();

  // scan through records
  SeqLib::BamRecord r;  
  while (reader_->GetNextRecord(r)) {
    
    // get read group tag
    std::string rg;
    if (!r.GetZTag("RG", rg)) {
      rg = "NA";
    }

    // already seen too many of these reads
    if (rg_count[rg] == sc.opts.perRgLearnLimit) {
      satisfied++;

      // compute isize stats
      bam_read_groups[rg].computeStats();

      // check if we satisifed all the read groups
      if (satisfied == groups.size()) {

	// store the max readlen and mapq for this entire bam across RGs
	for (const auto& br : bam_read_groups) {
	  readlen_max = std::max(readlen_max, br.second.readlen_max);
	  mapq_max = std::max(mapq_max, br.second.mapq_max);	  
	}
	break;
      }
    }
    
    // refer to existing BamReadGroup or make new
    auto& bstats = bam_read_groups[rg];
    bstats.addRead(r); // add the read for learning

    // update the counter for this group
    rg_count[rg]++;
    
  }
}

void BamReadGroup::addRead(const SeqLib::BamRecord &r)
{

  // track the flags
  ++reads;
  if (r.SecondaryFlag())
    ++supp;
  if (r.QCFailFlag())
    ++qcfail;
  if (r.DuplicateFlag())
    ++duplicate;
  if (!r.MappedFlag())
    ++unmap;
  if (!r.MateMappedFlag())
    ++mate_unmap;

  // track the mapq 
  mapq_max = std::max(mapq_max, r.MapQuality());
  
  // track the insert size
  int isizer = -1;
  if (!r.PairMappedFlag())
    isizer = -2;
  else if (!r.Interchromosomal())
    isizer = std::abs(r.InsertSize());
  isize_vec.push_back(isizer);

  // track the read length
  readlen_max = std::max(readlen_max, r.Length());

}

 void BamReadGroup::computeStats() {

   isize_mean = std::accumulate(isize_vec.begin(), isize_vec.end(), 0.0);

   // Calculate variance
   double sq_sum = 0.0;
   for (int val : isize_vec) {
     sq_sum += (val - isize_mean) * (val - isize_mean);
   }
   sd_isize = std::sqrt(sq_sum / isize_vec.size());  // population SD
   
 }
 
 
/*BamLearningResult LearnBamParams::learnAll(
    const std::map<std::string,std::string>& bamFiles,
    size_t perRGLimit,
    int sdCutoff,
    SvabaLogger& logger)
{

  // this is the master map of bams : (RG : params)
  BamLearningResult R;

  // loop through the BAM files and do the learning
  for (auto const& entry : sc.opts.bams) {
    
    auto const& name = entry.first;  // t001 etc
    auto const& path = entry.second; // filename

    sc.logger.log(opts.verbose > 1,true,"learning insert size for bam ", name, " at ", path);

    // the learn the RG params for a single BAM file
    LearnBamParams learner(path, sc.opts.perRGLimit, sc.logger);
    R.perFile[name] = learner.learnParams();

    // update global metrics
    for (auto const& grp : R.perFile[name]) {
      auto const& params = grp.second;
      R.globalReadLen = std::max(R.globalReadLen, params.readlen);
      R.globalMaxMapQ = std::max(R.globalMaxMapQ, params.max_mapq);
      int cutoff = int(std::floor(params.mean_isize 
                         + params.sd_isize * sdCutoff));
      R.globalMinDiscordantSize = 
        std::max(R.globalMinDiscordantSize, cutoff);
    }
  }

  return R;
}
*/
