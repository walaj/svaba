#include "LearnBamParams.h"

#include <numeric>
#include "SnowTools/BamWalker.h"

std::ostream& operator<<(std::ostream& out, const BamParams& p) {
 
  out << "@@@ READ GROUP " << p.read_group << " ISIZE " << p.mean_isize << "(" << p.sd_isize << ") MEAN COV " << p.mean_cov << " READLEN " << p.readlen << " MAX MAPQ " << p.max_mapq;
  // out << "@@@ Estimated fraction of reads that are discordant: " << p.frac_disc << std::endl;
  // out << "@@@ Estimated fraction of reads that are clipped:    " << p.frac_clip << std::endl;
  // out << "@@@ Estimated fraction of low quality alignments:    " << p.frac_bad << std::endl;
  // out << "@@@ Estimated mean coverage: " <<  p.mean_cov << std::endl;
  // out << "@@@ Read length: " << p.readlen << std::endl;
  // out << "@@@ Max mapping quality: " << p.max_mapq;
  return out;
}

//http://stackoverflow.com/questions/2114797/compute-median-of-values-stored-in-vector-c
double CalcMHWScore(std::vector<int>& scores)
{
  double median;
  size_t size = scores.size();

  std::sort(scores.begin(), scores.end());

  if (size  % 2 == 0)
    {
      median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
    }
  else 
    {
      median = scores[size / 2];
    }

  return median;
}

void LearnBamParams::learnParams(BamParamsMap& p, int max_count) {

  SnowTools::GenomicRegionVector grv = {SnowTools::GenomicRegion(0, 1000000,2000000), SnowTools::GenomicRegion(1,1000000,2000000)};
  SnowTools::BamWalker bwalker;
  bwalker.OpenReadBam(bam);
  bwalker.setBamWalkerRegions(grv);

  SnowTools::BamRead r;

  int wid = 0;
  double pos1 = 0, pos2 = 0, chr = -1;
  bool rule;

  // loop through a bunch of reads
  int count = 0;
  while (bwalker.GetNextRead(r, rule) && ++count < max_count) {

    bool qcpass = !r.DuplicateFlag() && !r.QCFailFlag() && !r.SecondaryFlag() && r.MappedFlag();
    if (!qcpass)
      continue;

    if (count % 500000 == 0)
      std::cerr << "...learning from read " << r.Brief(bwalker.header()) << " at count " << SnowTools::AddCommas(count) << " of " << SnowTools::AddCommas(max_count) << std::endl;

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
    
    // try to get the read group tag from qname first
    std::string qn = r.Qname();
    size_t posr = qn.find(":", 0);
    std::string RG = (posr != std::string::npos) ? qn.substr(0, posr) : r.GetZTag("RG");
    
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
      ff->second.isize_vec.push_back(r.InsertSize());

    ff->second.visited++;

  }

  wid += pos2 - pos1;
  // get the summary stats
  for (auto& i : p) {
    i.second.collectStats();
    i.second.mean_cov = (wid > 0) ? i.second.visited * i.second.readlen / wid : 0;
  }

}

void BamParams::collectStats() {

  if (isize_vec.size() < 100) {
    std::cerr << "!!!! NOT ENOUGH READS TO GET STATS !!!!" << std::endl;
    std::cerr << "!!!! " << read_group << " " << isize_vec.size() << " visited " << visited << std::endl;
  }
    
  // sort the isize vec
  std::sort(isize_vec.begin(), isize_vec.end());

  // get the 5% and 95% tiles
  int lp = isize_vec[std::floor(isize_vec.size() * 0.2)];
  int hp = isize_vec[std::floor(isize_vec.size() * 0.8)];

  // trim to only those in the 5-95% range
  std::vector<int> tmp;
  for (auto& i : isize_vec)
    if (i >= lp && i <= hp)
      tmp.push_back(i);
  isize_vec = tmp;

  // get the median first. Use to threshold
  median_isize = isize_vec.size() ? CalcMHWScore(isize_vec) : 0;

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

void LearnBamParams::learnParams(BamParams& p, int max_count) {

  SnowTools::GenomicRegionVector grv = {SnowTools::GenomicRegion(0, 1000000,2000000), SnowTools::GenomicRegion(1,1000000,2000000)};
  SnowTools::BamWalker bwalker;
  bwalker.OpenReadBam(bam);
  bwalker.setBamWalkerRegions(grv);

  SnowTools::BamRead r;

  double frac_clip = 0;
  double frac_disc = 0;
  double frac_bad =0;
  int mapq = 0;
  int readlen = 0;
  double cov_count = 0;
  int count = 0;

  double pos1 = 0, pos2 = 0;
  int chr1 = 0;
  bool rule;

  std::vector<int> isizer;

  while (bwalker.GetNextRead(r, rule) && ++count < max_count) {

    bool qcpass = !r.DuplicateFlag() && !r.QCFailFlag() && !r.SecondaryFlag();
    if (!qcpass)
      continue;

    readlen = std::max(r.Length(), readlen);
    mapq = std::max(r.MapQuality(), mapq);

    // get the first read position
    if (count == 1) {
      pos1 = r.Position();
      chr1 = r.ChrID();
    }
    
    if (abs(r.InsertSize() > 1000))
      ++frac_disc;
    if (r.NumClip() >= 5)
      ++frac_clip;
    if ( r.CigarSize() >= 6 || r.MapQuality() <= 5 || (r.QualitySequence().length() < 0.6 * readlen) ) {
      ++frac_bad;
    }
	
    if (r.InsertSize() > 0 && r.ProperOrientation() && r.InsertSize() < 10000)
      //&& r.ProperPair()) 
      isizer.push_back(r.InsertSize());

    if (count % 500000 == 0)
      std::cerr << "...learning from read " << r.Brief(bwalker.header()) << " at count " << SnowTools::AddCommas(count) << " of " << SnowTools::AddCommas(max_count) << std::endl;
    
    // last read position
    if (chr1 == r.ChrID() && r.MapQuality() > 0 && std::abs(r.InsertSize()) > 10 && std::abs(r.InsertSize()) < 2000) {
      ++cov_count;
      pos2 = r.Position();
    }
  }

  double sum = std::accumulate(isizer.begin(), isizer.end(), 0.0);
  p.mean_isize = isizer.size() > 0 ? sum / isizer.size() : 0;
  double mm = p.mean_isize;

  p.median_isize = isizer.size() ? CalcMHWScore(isizer) : 0;

  // get isize stdev
  std::vector<double> diff(isizer.size());
  std::transform(isizer.begin(), isizer.end(), diff.begin(), [mm](double x) { return x - mm; });
  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  p.sd_isize = std::sqrt(sq_sum / isizer.size());

  p.mean_cov = (pos2 - pos1 > 0) ? cov_count * readlen / (pos2 - pos1) : 0;
  p.frac_disc = frac_disc / (double)count;
  p.frac_bad = frac_bad / (double)count;
  p.frac_clip = frac_clip / (double)count;
  p.readlen = readlen;
  p.max_mapq = mapq;

}
