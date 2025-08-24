#pragma once

#include <ctime>
#include <sstream>
#include <unordered_map>
#include <map>
#include <chrono>
#include <tuple>

#include "svabaLogger.h"
#include "gzstream.h"

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/RefGenome.h"
#include "SeqLib/GenomicRegion.h"

typedef std::pair<size_t, size_t> CountPair; 
typedef std::tuple<size_t, size_t, std::string> Substring;
typedef std::vector<Substring> SubstringList;

#define SRTAG(r) ((r).GetZTag("SR") + "_" + std::to_string((r).AlignmentFlag()) + "_" + (r).Qname())

namespace svabaUtils {

  // make a structure to store timing opt
  struct svabaTimer {

    static const std::string header;
    
    svabaTimer();
    
    std::vector<std::string> s;
    
    //clock_t curr_clock;
    
    std::chrono::time_point<std::chrono::steady_clock> wall_start;
    double wall_elapsed = 0;
    
    void stop(const std::string& part);
    
    void start();
    
    // print it
    friend std::ostream& operator<<(std::ostream &out, const svabaTimer st);

    CountPair weird_read_count = {0,0}; // tumor/normal weird reads
    CountPair mate_read_count = {0,0}; // tumor/normal weird reads
    size_t dc_read_count = 0;
    size_t dc_cluster_count = 0;
    size_t contig_count = 0;
    size_t aligned_contig_count;
    size_t bps_count;

    std::string logRuntime(const SeqLib::BamHeader& h);
    
    // run time
    std::unordered_map<std::string, double> times;
    int pct_r = 0, pct_m = 0, pct_k = 0, pct_as = 0; // pct_pp = 0;

    SeqLib::GenomicRegion gr;
    
  };
  
  std::string myreplace(std::string &s,
			std::string toReplace,
			std::string replaceWith);
  
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

    if (!o) {
      std::cerr << "Failed to open file: " << s << std::endl;
      exit(EXIT_FAILURE);
    }
    
  }

  // Declare specialization for ogzstream
  template <>
  void fopen<ogzstream>(const std::string& s, ogzstream& o);

  std::string fileDateString();
  
  bool __header_has_chr_prefix(bam_hdr_t * h);
  
  /** Generate a weighed random integer 
   * @param cs Weighting for each integer (values must sum to one) 
   * @return Random integer bounded on [0,cs.size())
   */
  int weightedRandom(const std::vector<double>& cs);

  std::vector<std::string> tokenize_delimited(const std::string& str, char delim);

  void checkHeaderCompatibility(const SeqLib::BamHeader& bamHeader,
			      const SeqLib::BamHeader& refHeader,
				SvabaLogger& logger);

  std::vector<std::pair<int, int>> find_repeats(std::string_view seq, size_t single_repeat_count, size_t dinuc_repeat_count);

  std::vector<int> parsePLString(const std::string& pl_str);

  SubstringList find_long_dinuc_repeats(const std::string& s);

  SubstringList find_long_homopolymers(const std::string& s);  
  
}
