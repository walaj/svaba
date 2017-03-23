#ifndef SVABA_UTILS_H__
#define SVABA_UTILS_H__

#include <ctime>
#include <sstream>
#include <unordered_map>
#include <map>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/RefGenome.h"

#define SRTAG(r) ((r).GetZTag("SR") + "_" + std::to_string((r).AlignmentFlag()) + "_" + (r).Qname())

namespace svabaUtils {

// make a structure to store timing opt
struct svabaTimer {
  
  svabaTimer();

  std::unordered_map<std::string, double> times;
  std::vector<std::string> s;

  clock_t curr_clock;

  void stop(const std::string& part);

  void start();
  
  // print it
  friend std::ostream& operator<<(std::ostream &out, const svabaTimer st);
};

 double CalcMHWScore(std::vector<int>& scores);
 
 int overlapSize(const SeqLib::BamRecord& query, const SeqLib::BamRecordVector& subject);
 bool hasRepeat(const std::string& seq);
 std::string runTimeString(int num_t_reads, int num_n_reads, int contig_counter, 
			   const SeqLib::GenomicRegion& region, const SeqLib::BamHeader& h, const svabaTimer& st, 
			   const timespec& start);
 int countJobs(const std::string& regionFile, SeqLib::GRC &file_regions, SeqLib::GRC &run_regions, 
	       const SeqLib::BamHeader& h, int chunk, int window_pad);
 
  template <typename T>
  void fopen(const std::string& s, T& o) {
    o.open(s.c_str(), std::ios::out);
  }

  void print(std::stringstream& s, std::ofstream& log, bool cerr);

  std::string fileDateString();

  std::string __bamOptParse(std::map<std::string, std::string>& obam, std::istringstream& arg, int sample_number, const std::string& prefix);

  bool __openWriterBam(const SeqLib::BamHeader& h, const std::string& name, SeqLib::BamWriter& wbam);

  void __open_bed(const std::string& f, SeqLib::GRC& b, const SeqLib::BamHeader& h);

  bool __header_has_chr_prefix(bam_hdr_t * h);

  bool __open_index_and_writer(const std::string& index, SeqLib::BWAWrapper * b, const std::string& wname, SeqLib::BamWriter& writer, SeqLib::RefGenome *& r, SeqLib::BamHeader& bwa_header);

  /** Generate a weighed random integer 
   * @param cs Weighting for each integer (values must sum to one) 
   * @return Random integer bounded on [0,cs.size())
   */
  int weightedRandom(const std::vector<double>& cs);
  
}


#endif
