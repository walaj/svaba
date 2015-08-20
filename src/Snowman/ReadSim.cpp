#include "ReadSim.h"

#include <iostream>
#include <random>
#include <cassert>

#include "SnowTools/SnowUtils.h"

void ReadSim::addAllele(const std::string& s, double af) {
  
  if (m_seq.size() && m_seq.back().length() != s.length()) {
    std::cerr << "ReadSim::addAllele error. Expecting sequences to be same size" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (af <= 0 || af > 1) {
    std::cerr << "ReadSim::addAllele error. Expecting allelic fraction > 0 and <= 1. value is " << af << std::endl;
    exit(EXIT_FAILURE);
  }
    
  m_seq.push_back(s);
  m_frac.push_back(af);

}

int ReadSim::__random_allele(std::vector<double>& cs) const {
  
  // get a random allele
  size_t al = 0;
  uint32_t rval = rand() % 100000;
  
  double rand_allele = rval % 100000;
  while (al < cs.size()) {
    if (rand_allele <= cs[al] * 100000) 
      break;
    ++al;
  }
  
  return al;

}

double ReadSim::__get_cumsum(std::vector<double>& cs) const {

  double csum = 0;
  for (size_t i = 0; i < m_frac.size(); ++i) {
    csum += m_frac[i];
    cs.push_back(csum);
  }
  
  if (csum < 0.999 || csum > 1.001) {
    std::cerr << "ReadSim::addSequence: Expecting sum of allelic fractions = 1. Sum is " << csum << std::endl;
    return csum;
  }
  return csum;
}

void ReadSim::sampleReadsToCoverage(std::vector<std::string>& reads, int cov, 
				    double error_rate, double ins_error_rate, double del_error_rate, 
				    int readlen) {

  // check validity
  if (m_seq.size() == 0) {
    std::cerr << "ReadSim::sampleReadsToCoverage: No sequences. Add with ReadSim::addSequence" << std::endl;
    return;
  }

  // how many smaples to get desired coverage?
  size_t seqlen = m_seq[0].length();
  int ns = cov * seqlen / readlen;

  uint32_t maxstart = seqlen - readlen;

  // get cumulative sum for sampling
  std::vector<double> m_frac_cumsum;
  __get_cumsum(m_frac_cumsum);

  // loop through num samples and start sampling
  for (int i = 0; i < ns; ++i) {

    // random allele
    int al = __random_allele(m_frac_cumsum);

    // get a random subsequence
    //uint32_t rstart;
    uint32_t rstart = rand() % maxstart;
    std::string s = m_seq[al].substr(rstart, readlen);
    
    // add the errors
    makeSNVErrors(s, error_rate);

    // add the insertion errors
    if (rand() % 100000 > ins_error_rate*100000)
      makeInsErrors(s);

    // only make sometimes
    if (rand() % 100000 > del_error_rate*100000)
      makeDelErrors(s, rstart, m_seq[al]);

    reads.push_back(s);
  }
  
}

void ReadSim::samplePairedEndReadsToCoverage(std::vector<std::string>& reads1, std::vector<std::string>& reads2, 
					     int cov, double error_rate, double ins_error_rate, double del_error_rate, 
					     int readlen, 
					     double mean_isize, double sd_isize) {
  // check validity
  if (m_seq.size() == 0) {
    std::cerr << "ReadSim::samplePairedEndReadsToCoverage: No sequences. Add with ReadSim::addSequence" << std::endl;
    return;
  }

  // how many samples to get desired coverage?
  size_t seqlen = m_seq[0].length();
  int ns = cov * seqlen / readlen / 2;

  // get cumulative sum for sampling
  std::vector<double> m_frac_cumsum;
  __get_cumsum(m_frac_cumsum);

  // setup random number generator for isize
  std::default_random_engine generator(rand());
  std::normal_distribution<double> distribution(mean_isize, sd_isize);
  
  // loop through num samples and start sampling
  for (int i = 0; i < ns; ++i) {

    if (i % 100000 == 0)
      std::cerr << "...sampling read " << i << " of " << ns << std::endl;

    // get a random isize
    int isize = distribution(generator);

    uint32_t maxstart = seqlen - readlen*2 - isize;

    // random allele
    int al = __random_allele(m_frac_cumsum);

    // get random subsequences
    uint32_t rstart = rand() % maxstart;
    std::string s1 = m_seq[al].substr(rstart, readlen);
    std::string s2 = m_seq[al].substr(rstart + readlen + isize, readlen);

    // dont keep reads with N
    if (s1.find("N") != std::string::npos || s2.find("N") != std::string::npos)
      continue;

    SnowTools::rcomplement(s2);

    // add the errors
    //std::cerr << "      " << s1 << " " << rstart << " rlen " << readlen << " ms " << m_seq[al].length() << std::endl;
    
    if (error_rate > 0) {
      makeSNVErrors(s1, error_rate);
      makeSNVErrors(s2, error_rate);
    }

    if (ins_error_rate) {
      if (rand() % 100000 < ins_error_rate*100000)
	makeInsErrors(s1);
      if (rand() % 100000 < ins_error_rate*100000)
	makeInsErrors(s2);
    }

    if (del_error_rate) {
      if (rand() % 100000 < del_error_rate*100000)
	makeDelErrors(s1, rstart, m_seq[al]);
      if (rand() % 100000 < del_error_rate*100000)
	makeDelErrors(s2, rstart + readlen + isize, m_seq[al]);
    }

    reads1.push_back(s1);
    reads2.push_back(s2);
  }  
}

std::ostream& operator<<(std::ostream& out, const Indel& i) {
  out << i.len << "\t" << i.type << "\t" << i.gr.chr << "\t" << i.gr.pos1 << "\t" << i.gr.pos2 << "\t" << i.frag_id << "\t" << (i.lead_base + i.ref_seq) << "\t" << (i.lead_base + i.alt_seq);
  return out;
}

Indel ReadSim::makeDelErrors(std::string& s) {

  const std::vector<int> sizer = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
				  2,2,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,4,5,5,5,6,6,
				  7,8,9,10,11,12,13,14,15,16,17,18,20,22,24,26};

  uint32_t rpos = s.length() + 1;
  // get a random size
  int ds = getRandomIndelSize();

  size_t failsafe = 0;
  while (rpos + ds > s.length() && failsafe < 1000) {
    
    ++failsafe;
    // get a random position in read
    rpos = rand() % (s.length() - 5 - ds);
    rpos += 5;

  }

  // get the replacement sequence
  s = s.substr(0, rpos) + s.substr(rpos + ds, s.length() - rpos - ds); // + refseq.substr(rpos + s.length(), ds);

  return Indel(); //{ds, 'D', rpos};


}
Indel ReadSim::makeDelErrors(std::string& s, int sstart, const std::string& refseq) {
  
  const std::vector<int> sizer = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
				  2,2,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,4,5,5,5,6,6,
				  7,8,9,10,11,12,13,14,15,16,17,18,20,22,24,26};

  int ds = getRandomIndelSize();

  // get a random position in read
  uint32_t rpos = rand() % (s.length() - 5 - ds);
  rpos += 5;

  // make sure we aren't at end
  if (sstart + s.length() + ds >= refseq.length())
    return Indel(); // empty

  // get the replacement sequence
  s = s.substr(0, rpos) + s.substr(rpos + ds, s.length() - rpos - ds) + refseq.substr(sstart + s.length(), ds);

  //return Indel(ds, 'D', );
  return Indel();
  
}

int ReadSim::getRandomIndelSize() const {

  const std::vector<int> s = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
				  2,2,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,4,5,5,5,6,6,
				  7,8,9,10,11,12,13,14,15,16,17,18,20,22,24,26};
  
  size_t rr = rand() % s.size();
  return s[rr];
}

Indel ReadSim::makeInsErrors(std::string& s, bool keep_size) {


  const std::vector<int> sizer = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
				  2,2,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,4,5,5,5,6,6,
				  7,8,9,10,11,12,13,14,15,16,17,18,20,22,24,26};

  const char TCGA[5] = "TCGA";

  // get a random size
  int is = getRandomIndelSize();

  // generate the random insertion piece
  std::string ins(is, 'N');
  for (int i = 0; i < is; ++i) 
    ins[i] = TCGA[rand() % 3];
  
  // get a random position in read
  uint32_t rpos = rand() % (s.length() - 10);
  assert(s.length() - 10  > 0);

  rpos += 10;

  if (keep_size) {
    if (rpos < s.length() / 2) { // trim off back end
      s = s.substr(0, rpos) + ins + s.substr(rpos, s.length() - rpos - is);
    } else { // trim off front end
      s = s.substr(is, rpos - is) + ins + s.substr(rpos, s.length() - rpos + is); 
    }
  } else {
    s = s.substr(0, rpos) + ins + s.substr(rpos, s.length() - rpos);
  }

  return Indel(); //Indel(is, 'I', rpos);
  
}

void ReadSim::makeSNVErrors(std::string& s, double er) {

  const char T[4] = "ACG";
  const char C[4] = "ATG";
  const char G[4] = "ACT";
  const char A[4] = "TCG";

  std::default_random_engine generator(rand());
  std::binomial_distribution<int> distribution(s.length(), er);
  size_t num_errors = distribution(generator);

  for (size_t i = 0; i < num_errors; ++i) {

    uint32_t pos = rand() % s.length();
    size_t w = rand() % 3;

    if (s.at(pos) == 'A')
      s[pos] = A[w];
    else if (s.at(pos) == 'T')
      s[pos] = T[w];
    else if (s.at(pos) == 'C')
      s[pos] = C[w];
    else if (s.at(pos) == 'G')
      s[pos] = G[w];
    else if (s.at(pos) == 'N')
      ;
    else
      std::cerr << "ReadSim::makeSNVErrors: Unexpected character in string. Char: " << s.at(pos) << " main seq " << s <<  std::endl;

    
  }

}

void ReadSim::baseQualityRelevantErrors(std::string& s, const std::string& bq) {


  char TCGA[5] = "TCGA";

  assert(s.length() == bq.length());
  for (size_t i = 0; i < s.length(); ++i) {
    if (bq.at(i) < 35) { // low quality-ish
      s[i] = TCGA[rand() % 4];
    }
  }
  
}

