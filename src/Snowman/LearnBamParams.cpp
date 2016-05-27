#include "LearnBamParams.h"

#include "SnowTools/BamWalker.h"

std::ostream& operator<<(std::ostream& out, const BamParams& p) {
 
  out << "@@@ Estimated fraction of reads that are discordant: " << p.frac_disc << std::endl;
  out << "@@@ Estimated fraction of reads that are clipped:    " << p.frac_clip << std::endl;
  out << "@@@ Estimated fraction of low quality alignments:    " << p.frac_bad << std::endl;
  out << "@@@ Estimated mean coverage: " <<  p.mean_cov << std::endl;
  out << "@@@ Read length: " << p.readlen << std::endl;
  out << "@@@ Max mapping quality: " << p.max_mapq;
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
