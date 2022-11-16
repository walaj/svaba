#include "LearnBamParams.h"

#include <numeric>
#include "SeqLib/BamReader.h"
#include "svabaUtils.h"

std::ostream& operator<<(std::ostream& out, const BamParams& p) {

  float frac_clip = (float)p.num_clip/p.visited;
  
  out << "\t@@@ READ GROUP: " << p.read_group << " - Insert Size: " << p.mean_isize 
      << "(SD: " << p.sd_isize << ") [0.025%,97.5%] [" << p.lp << "," 
      << p.hp << "]" << std::endl
      << " Mean Coverage: " << p.mean_cov << " Read Length: " 
      << p.readlen << " Max MapQ:" << p.max_mapq << std::endl;
  out << "\t\tReads Sampled: " << SeqLib::AddCommas(p.visited) << " Num Clipped(frac): "
      << SeqLib::AddCommas(p.num_clip) << "(" << frac_clip << ")" << std::endl;
    
  // out << "@@@ Estimated fraction of reads that are discordant: " << p.frac_disc << std::endl;
  // out << "@@@ Estimated fraction of reads that are clipped:    " << p.frac_clip << std::endl;
  // out << "@@@ Estimated fraction of low quality alignments:    " << p.frac_bad << std::endl;
  // out << "@@@ Estimated mean coverage: " <<  p.mean_cov << std::endl;
  // out << "@@@ Read length: " << p.readlen << std::endl;
  // out << "@@@ Max mapping quality: " << p.max_mapq;
  return out;
}

void LearnBamParams::learnParams(BamParamsMap& p, int max_count) {

  SeqLib::BamReader bwalker;
  assert(bwalker.Open(bam));

  SeqLib::BamRecord r;

  int wid = 0; // Total span of genome tested
  double pos1 = 0, pos2 = 0, chr = -1;

  // read from the beginning until hit max_count
  // pos1 is set by process_read as lowest position on chr, pos2 is highest
  // wid is the max span tested and is updated with each chr
  size_t count = 0;
  pos1 = 0; pos2 = 0; chr = -1; wid = 0;
  bwalker.Reset();
  while (bwalker.GetNextRecord(r) && ++count < max_count)
    process_read(r, count, p, pos1, pos2, chr, wid);
  wid += pos2 - pos1; // add in the last chr width
  
  // get the summary stats
  for (auto& i : p) {
    i.second.collectStats();
    i.second.mean_cov = (wid > 0) ? i.second.visited * i.second.readlen / wid : 0;
  }
  
}

// go back through and estimate how many read pairs will be labeled as discordant in this BAM file
void LearnBamParams::estimateDiscordant(BamParamsMap& p, int max_count, int il_cutoff) {

  SeqLib::BamReader bwalker;
  assert(bwalker.Open(bam));

  SeqLib::BamRecord r;
  size_t count = 0;
  
  while (bwalker.GetNextRecord(r) && ++count < max_count) {
    if (r.DuplicateFlag() || r.QCFailFlag() || r.SecondaryFlag() || !r.MappedFlag())
      return;

    std::string RG;
    get_rg(RG, r);
    
    BamParamsMap::iterator ff = p.find(RG);
    
    // new read group, make a new object here
    if (ff == p.end()) {
      p[RG] = BamParams(RG);
      ff = p.find(RG);
    }

    if (r.InsertSize() >= il_cutoff)
      ff->second.num_disc++;
  }
  
}

float BamParams::discFrac() {
  return (float)num_disc / visited;
}

void BamParams::collectStats() {

  if (isize_vec.size() < 100) {
    std::cerr << "not enough paired-end reads to get insert-size distribution. skipping discordant analysis" << std::endl;
    std::cerr << "read group: " << read_group << " Insert-sizes sampled: " << isize_vec.size() << " visited " << visited << std::endl;
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

  //frac_disc = num_disc / (double)visited;
  //frac_bad  = num_bad / (double)visited;
  //frac_clip = num_clip / (double)visited;


}

void LearnBamParams::process_read(const SeqLib::BamRecord& r, size_t count, BamParamsMap& p,
				  double& pos1, double& pos2, double& chr, int& wid) const {
   
   if (r.DuplicateFlag() || r.QCFailFlag() || r.SecondaryFlag() || !r.MappedFlag())
     return;
   
   // get the first read position
   if (count == 1 || r.ChrID() != chr) {
     chr = r.ChrID();
     if (count > 1)
       wid += pos2 - pos1;
     pos1 = r.Position();
   }
   
   // update last read position (stay on same chrom as first)
   if (chr == r.ChrID())
     pos2 = r.Position();
   
   std::string RG;
   get_rg(RG, r);
   
   BamParamsMap::iterator ff = p.find(RG);
   
   // new read group, make a new object here
   if (ff == p.end()) {
     p[RG] = BamParams(RG);
     ff = p.find(RG);
    }
   
   // max readlen and mapq
   ff->second.readlen = std::max(r.Length(), ff->second.readlen);
   ff->second.max_mapq = std::max(r.MapQuality(), ff->second.max_mapq);
   
   // count number of clips
   if (r.NumClip() >= 5) ff->second.num_clip++;
   
   // collect all of the insert sizes
   if (r.InsertSize() > 0 && r.ProperOrientation())
     ff->second.isize_vec.push_back(r.FullInsertSize());
   
   ff->second.visited++;
   
}

void LearnBamParams::get_rg(std::string& rg, const SeqLib::BamRecord& r) const {
  if (!r.GetZTag("RG", rg))
    rg = "NA";
  
  //debug
  //std::cerr << "RG " << rg << std::endl;
  
  // hack for simulated data
  if (rg.find("tumor") != std::string::npos) {
    std::string qn = r.Qname();
    size_t posr = qn.find(":", 0);
    rg = (posr != std::string::npos) ? qn.substr(0, posr) : rg;
  } else {
    // best practice without "tumor" hack
    //rg = r.ParseReadGroup();
  }
  
  //debug
  //std::cerr << "   RG " << rg << std::endl; //debug
  
}
