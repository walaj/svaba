#include "LearnBamParams.h"

#include "SnowTools/BamWalker.h"

std::ostream& operator<<(std::ostream& out, const BamParams& p) {
  
  out << "@@@ Estimated fraction of reads that are discordant: " << p.frac_disc << std::endl;
  out << "@@@ Estimated fraction of reads that are clipped:    " << p.frac_clip << std::endl;
  out << "@@@ Estimated fraction of reads that are bad:    " << p.frac_bad << std::endl;
  out << "@@@ Estimated mean coverage: " <<  p.mean_cov << std::endl;
  out << "@@@ Read length: " << p.readlen << std::endl;
  out << "@@@ Max mapping quality: " << p.max_mapq;
  return out;
}

void LearnBamParams::learnParams(BamParams& p, int max_count) {

  SnowTools::BamWalker bwalker;
  bwalker.OpenReadBam(bam);

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
    if ( r.CigarSize() >= 6 || r.MapQuality() <= 5 || (r.QualitySequence().length() < 0.6 * readlen) )
      ++frac_bad;
	
    if (count % 500000 == 0)
      std::cerr << "...learning from read " << r.Brief(bwalker.header()) << " at count " << SnowTools::AddCommas(count) << " of " << SnowTools::AddCommas(max_count) << std::endl;
    
    // last read position
    if (chr1 == r.ChrID() && r.MapQuality() > 0 && std::abs(r.InsertSize()) > 10 && std::abs(r.InsertSize()) < 2000) {
      ++cov_count;
      pos2 = r.Position();
    }
  }
  
  p.mean_cov = (pos2-pos1 > 0) ? cov_count * readlen / (pos2 - pos1) : 0;
  p.frac_disc = frac_disc / (double)count;
  p.frac_bad = frac_bad / (double)count;
  p.frac_clip = frac_clip / (double)count;
  p.readlen = readlen;
  p.max_mapq = mapq;
  
}
