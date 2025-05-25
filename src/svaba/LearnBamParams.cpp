#include "LearnBamParams.h"

#include <numeric>
#include "SeqLib/BamReader.h"
#include "svabaUtils.h"

std::vector<std::string> extractReadGroups(const SeqLib::BamHeader& hdr) {
  std::vector<std::string> rgs;
  std::string text = hdr.AsString();               // full header, newline separated
  std::istringstream lines(text);
  std::string line;

  while (std::getline(lines, line)) {
    // look for lines that start with "@RG"
    if (line.rfind("@RG", 0) == 0) {
      // split on tabs to find the ID: field
      std::istringstream fields(line);
      std::string field;
      while (std::getline(fields, field, '\t')) {
        if (field.rfind("ID:", 0) == 0) {
          rgs.push_back(field.substr(3));  // drop the "ID:" prefix
          break;
        }
      }
    }
  }

  return rgs;
}

std::ostream& operator<<(std::ostream& out, const BamParams& p) {
 
  out << "@@@ READ GROUP " << p.read_group << " Insert Size: " << p.mean_isize 
      << "(" << p.sd_isize << ") [0.025%,97.5%] [" << p.lp << "," 
      << p.hp << "], Mean Coverage: " << p.mean_cov << " Read Length: " 
      << p.readlen << " Max MapQ:" << p.max_mapq;
  // out << "@@@ Estimated fraction of reads that are discordant: " << p.frac_disc << std::endl;
  // out << "@@@ Estimated fraction of reads that are clipped:    " << p.frac_clip << std::endl;
  // out << "@@@ Estimated fraction of low quality alignments:    " << p.frac_bad << std::endl;
  // out << "@@@ Estimated mean coverage: " <<  p.mean_cov << std::endl;
  // out << "@@@ Read length: " << p.readlen << std::endl;
  // out << "@@@ Max mapping quality: " << p.max_mapq;
  return out;
}

void LearnBamParams::learnParams(BamParamsMap& p) {

  // make sure the bam is open
  SeqLib::GRC grv;
  SeqLib::BamReader bwalker;
  assert(bwalker.Open(bam));

  // find all of the read groups in the header
  std::vector<std::string> all_rgs = extractReadGroups(bwalker.Header());

  // initialize a counter for each RG
  for (auto& rg : all_rgs) {
    rg_counts[rg] = 0;
  }
  size_t satisfied = 0; // how many RGs have hit learning limit
  
  // read from the beginning until hit max_count
  bwalker.Reset();
  SeqLib::BamRecord r;  
  while (bwalker.GetNextRecord(r)) {
    ++num_reads_seen;
    process_read(r, p, satisfied);

    // if this is taking a while, let user know
    if (num_reads_seen > 200000000 && num_reads_seen % 200000000 == 0) {
      std::cerr << "...extensive BAM reading for insert size learning, RGs not homogenous - " << std::endl << 
	"    " << SeqLib::AddCommas(num_reads_seen) << " reads viewed on bam: " << bam << std::endl;
    }
    
    // if we hit all the learning limits, break
    if (satisfied == rg_counts.size()) {
      break;
    }
  }

  // get the summary stats
  for (auto& i : p) {
    i.second.collectStats();
  }
  
}

void BamParams::collectStats() {

  if (isize_vec.size() < 100) {
    std::cerr << "not enough paired-end reads to get insert-size distribution. skipping discordant analysis" << std::endl;
    std::cerr << "\read group: " << read_group << " Insert-sizes sampled: " << isize_vec.size() << " visited " << visited << std::endl;
    return;
  }
    
  // sort the isize vec
  std::sort(isize_vec.begin(), isize_vec.end());

  // get the 5% and 95% tiles
  lp = isize_vec[std::floor(isize_vec.size() * 0.025)];
  hp = isize_vec[std::floor(isize_vec.size() * 0.975)];

  // trim to only those in the 5-95% range
  std::vector<int> tmp;
  for (auto& i : isize_vec)
    if (i >= lp && i <= hp)
      tmp.push_back(i);
  isize_vec = tmp;

  // get the median first. Use to threshold
  median_isize = isize_vec.size() ? svabaUtils::CalcMHWScore(isize_vec) : 0;

  double sum = std::accumulate(isize_vec.begin(), isize_vec.end(), 0.0);
  mean_isize = isize_vec.size() > 0 ? sum / isize_vec.size() : 0;

  // get isize stdev
  double mm = mean_isize;
  std::vector<double> diff(isize_vec.size());
  std::transform(isize_vec.begin(), isize_vec.end(), diff.begin(), [mm](double x) { return x - mm; });
  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  sd_isize = std::sqrt(sq_sum / isize_vec.size());

}

void LearnBamParams::process_read(const SeqLib::BamRecord& r, BamParamsMap& p, size_t& satisfied) {

  // exclude some read types
  if (r.DuplicateFlag() || r.QCFailFlag() || r.SecondaryFlag() || !r.MappedFlag())
    return;
  
  // get the RG
  std::string RG;
  if (!r.GetZTag("RG", RG))
    RG = "NA";
  
  // update the read group seen counter
  auto it = rg_counts.find(RG);

  // warn user if read group in read is not in header
  // but can proceed anyway afterwards
  if (it == rg_counts.end()) {
    std::cerr << "[WARNING] saw unexpected read group \"" << RG << "\"; known RGs:";
    for (auto const & kv : rg_counts)
      std::cerr << " " << kv.first;
    std::cerr << "\n";
    
    // add it in since missing
    rg_counts[RG] = 0;
    it = rg_counts.find(RG); // and find the new slot
  }

  // now count the read group and check if past limit
  assert(it != rg_counts.end());
  if (it->second == per_rg_limit)
    return; // we've already learned enough for this RG
  ++(it->second);
  if (it->second == per_rg_limit)
    satisfied++;
  
  // find which RG and make new if needed
  BamParamsMap::iterator ff = p.find(RG);
  if (ff == p.end()) {
    p[RG] = BamParams(RG);
    ff = p.find(RG);
  }
  
  // reacheck the max readlen and mapq for that readgroup
  ff->second.readlen = std::max(r.Length(), ff->second.readlen);
  ff->second.max_mapq = std::max(r.MapQuality(), ff->second.max_mapq);
  
  // count number of clips
  if (r.NumClip() >= 5) ff->second.num_clip++;
  
  // add the insert size to the collection for stats later
  if (r.InsertSize() > 0 && r.ProperOrientation())
    ff->second.isize_vec.push_back(r.FullInsertSize());

  // record visited
  ff->second.visited++;
  
}
