// LearnBamParams.cpp

#include "LearnBamParams.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"

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

LearnBamParams::LearnBamParams(const SvabaSharedConfig& sc,
  const std::string& bamPath) {
  : sc_(sc),
    bam_(bamPath)
{

  // open the BAM for learning
  if (!reader_.Open(bam_)) {
    throw std::runtime_error("Could not open bam file " + bam_ + " in LearnBamParams");
  }
}

BamParamsMap LearnBamParams::learnParams() {
  BamParamsMap result;

  // gather read group ids
  auto hdr = reader_.Header();
  auto groups = extractReadGroups(hdr);

  // init count for each group
  std::map<std::string,size_t> counts;
  for (auto const& g : groups) {
    counts[g] = 0;
  }

  size_t satisfied = 0;
  SeqLib::BamRecord r;
  reader_.Reset();

  // scan through records
  while (reader_.GetNextRecord(r)) {
    // skip flagged or unmapped
    if (r.DuplicateFlag() ||
        r.QCFailFlag()   ||
        r.SecondaryFlag()||
       !r.MappedFlag())
      continue;

    // get read group tag
    std::string rg;
    if (!r.GetZTag("RG", rg)) {
      rg = "NA";
    }

    // ensure we have a count slot
    auto itc = counts.find(rg);
    if (itc == counts.end()) {
      counts[rg] = 0;
      itc = counts.find(rg);
    }

    // only sample up to the per group limit
    if (itc->second < per_rg_limit_) {
      auto& params = result[rg];
      if (params.read_group.empty()) {
        params.read_group = rg;
      }
      // update basic stats
      params.readlen   = std::max(params.readlen,   r.Length());
      params.max_mapq  = std::max(params.max_mapq,  r.MapQuality());
      if (r.NumClip() >= 5) {
        params.num_clip++;
      }
      if (r.InsertSize() > 0 && r.ProperOrientation()) {
        params.isize_vec.push_back(r.FullInsertSize());
      }
      // increment count
      itc->second++;
      if (itc->second == per_rg_limit_) {
        satisfied++;
      }
    }

    // break if all groups have enough reads
    if (satisfied == counts.size()) {
      break;
    }
  }

  // compute derived stats for each group
  for (auto& kv : result) {
    kv.second.collectStats();
  }

  return result;
}

BamLearningResult LearnBamParams::learnAll(
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
