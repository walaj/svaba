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

#include <ostream>

std::ostream& operator<<(std::ostream& os, const BamReadGroup& bg) {
  os << "r="  << bg.reads
     << " s="  << bg.supp
     << " u="  << bg.unmap
     << " qc=" << bg.qcfail
     << " d="  << bg.duplicate
     << " mu=" << bg.mate_unmap
     << " MQ=" << bg.mapq_max
     << " RL=" << bg.readlen_max
     << " mIS=" << bg.isize_mean
     << " sdIS=" << bg.sd_isize;
  return os;
}


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
  size_t countr = 0;

  // track how many read groups are satisfied. Break when all are done
  size_t satisfied = 0;
  
  // set the BAM reader to the beginning of the BAM
  reader_->Reset();

  // scan through records
  while (auto r = reader_->Next()) {

    // get read group tag
    std::string rg;
    if (!r->GetZTag("RG", rg)) {
      rg = "NA";
    }

    // print update
    ++countr;
    if (countr % 10000000==0) {
      sc.logger.log(true,true,"......learning for read ",
		    SeqLib::AddCommas(countr), " and learned ",
		    satisfied, " read groups of ", groups.size());
      
      for (const auto& rgc : rg_count)
	if (rgc.second < 100000)
	  std::cerr << rgc.first << ":" << rgc.second << std::endl;
      
    }
    
    // already seen too many of these reads
    if (rg_count[rg] == sc.opts.perRgLearnLimit) {
      satisfied++;

      // check if we satisifed all the read groups
      if (satisfied == groups.size()) {
	break;
      }
    }
    
    // refer to existing BamReadGroup or make new
    auto& bstats = bam_read_groups[rg];
    bstats.addRead(*r); // add the read for learning

    // update the counter for this group
    rg_count[rg]++;
    
  }
  
  // compute isize stats
  // store the max readlen and mapq for this entire bam across RGs
  for (auto& br : bam_read_groups) {
    readlen_max = std::max(readlen_max, br.second.readlen_max);
    mapq_max = std::max(mapq_max, br.second.mapq_max);	  
    br.second.computeStats();
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
  if (!r.PairMappedFlag() || r.Interchromosomal() || r.PairOrientation() != FRORIENTATION)
    ;
  else if (!r.Interchromosomal())
    isize_vec.push_back(std::abs(r.InsertSize()));

  // track the read length
  readlen_max = std::max(readlen_max, r.Length());

}

 void BamReadGroup::computeStats() {

   if (isize_vec.empty()) {
     isize_mean = 0.0;
     sd_isize   = 0.0;
     return;
   }
   
   // Remove the top 2% of values to filter out extreme outliers
   std::sort(isize_vec.begin(), isize_vec.end());
   size_t n = isize_vec.size();
   size_t keep = static_cast<size_t>(std::floor(n * 0.98));
   if (keep == 0) {
     // If 95% rounds down to 0, keep at least one element
     keep = 1;
   }
   isize_vec.resize(keep);
   
   // Calculate mean
   isize_mean = std::accumulate(isize_vec.begin(), isize_vec.end(), 0.0) / keep;
    
    // Calculate population standard deviation
    double sq_sum = 0.0;
    for (int val : isize_vec) {
      double diff = val - isize_mean;
        sq_sum += diff * diff;
    }
    sd_isize = std::sqrt(sq_sum / keep);
    
    // Clean up memory
    isize_vec.clear();
    
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
