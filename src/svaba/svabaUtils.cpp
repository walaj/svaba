#include "svabaUtils.h"

#include <iomanip>
#include <filesystem>  // C++17
namespace fs = std::filesystem;

namespace svabaUtils {

static std::string POLYA = "AAAAAAAAAAAAAAAAAAA";
static std::string POLYT = "TTTTTTTTTTTTTTTTTTT";
static std::string POLYC = "CCCCCCCCCCCCCCCCCCC";
static std::string POLYG = "GGGGGGGGGGGGGGGGGGG";
static std::string POLYAT = "ATATATATATATATATATATATAT";
static std::string POLYCG = "CGCGCGCGCGCGCGCGCGCGCGCG";
static std::string POLYTG = "TGTGTGTGTGTGTGTGTGTGTGTG";
static std::string POLYCA = "CACACACACACACACACACACACA";
static std::string POLYAG = "AGAGAGAGAGAGAGAGAGAGAGAG";
static std::string POLYTC = "TCTCTCTCTCTCTCTCTCTCTCTC";


  const std::string svabaTimer::header =
  "chromosome\tstart\tend\t"
  "tumor_weird_reads\tnormal_weird_reads\t"
  "tumor_mate_reads\tnormal_mate_reads\t"    
  "discordant_reads\tdiscordant_clusters\t"
  "contigs\tcontigs_pass\tbps\truntime_seconds\t"
  "pct_r\tpct_m\tpct_k\tpct_as";

  std::string svabaTimer::logRuntime(const SeqLib::BamHeader& h) {
    double total_time = 0;
    for (const auto& i : times)
      total_time += i.second;
    
    if (total_time == 0)
      total_time = 1.0;  // Avoid division by zero
    
    auto get_pct = [&](const std::string& key) -> int {
      auto it = times.find(key);
      return (it != times.end()) ? SeqLib::percentCalc<double>(it->second, total_time) : 0;
    };
    
    // Store percentages
    pct_r  = get_pct("r");
    pct_m  = get_pct("m");
    pct_k  = get_pct("k");
    pct_as = get_pct("as");
    //pct_pp = get_pct("pp");
    
    std::ostringstream oss;
    oss << gr.ChrName(h) << "\t"
	<< gr.pos1 << "\t" << gr.pos2 << "\t"
	<< weird_read_count.first << "\t" << weird_read_count.second << "\t"
	<< mate_read_count.first << "\t" << mate_read_count.second << "\t"
	<< dc_read_count << "\t" << dc_cluster_count << "\t"
	<< contig_count << "\t" << aligned_contig_count << "\t"
	<< bps_count << "\t" << total_time << "\t"
	<< pct_r  << "\t" << pct_m  << "\t"
	<< pct_k  << "\t" << pct_as;
    
    return oss.str();
  }
  
  
  std::string fileDateString() {
    // set the time string
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    std::stringstream mdate;
    mdate << (now->tm_year + 1900) 
	  << (now->tm_mon + 1 < 10 ? "0": "") << (now->tm_mon + 1) 
	  << (now->tm_mday < 10 ? "0" : "")
	  << now->tm_mday;
    return mdate.str();
  }

  svabaTimer::svabaTimer() {
    s = {"r", "m", "as", "bw", "pp", "k"};
    for (auto& i : s)
      times[i] = 0;
  }

  void svabaTimer::stop(const std::string& part) {
    auto wall_end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(wall_end - wall_start).count();
    times[part] += elapsed;
    wall_elapsed += elapsed;    
    //times[part] += static_cast<double>(clock() - curr_clock) / CLOCKS_PER_SEC;
    //curr_clock = clock();
  }

  void svabaTimer::start() {
    wall_start = std::chrono::steady_clock::now();
    //curr_clock = clock(); 
  }

  std::ostream& operator<<(std::ostream &out, const svabaTimer st) {

    double total_time = 0;
    for (auto& i : st.times)
      total_time += i.second;
    if (total_time == 0)
      return out;

    char buffer[140];
    
    auto itr = st.times.find("r");
    auto itm = st.times.find("m");
    auto ita = st.times.find("as");
    auto itk = st.times.find("k");
    auto itt = st.times.find("t");
    auto itp = st.times.find("pp");

    if (total_time)
      snprintf (buffer, 140, "R: %2d%% M: %2d%% T: %2d%% C: %2d%% A: %2d%% P: %2d%%", 
	       SeqLib::percentCalc<double>(itr->second, total_time),
	       SeqLib::percentCalc<double>(itm->second, total_time),
	       SeqLib::percentCalc<double>(itt->second, total_time),
	       SeqLib::percentCalc<double>(itk->second, total_time),
	       SeqLib::percentCalc<double>(ita->second, total_time),
	       SeqLib::percentCalc<double>(itp->second, total_time));
    else
      snprintf (buffer, 140, "NO TIME");
    out << std::string(buffer);
    return out;
  }


// useful replace function
  std::string myreplace(std::string &s,
			std::string toReplace,
			std::string replaceWith)
  {
    return(s.replace(s.find(toReplace), toReplace.length(), replaceWith));
  }
  
bool hasRepeat(const std::string& seq) {
  
  if ((seq.find(POLYT) == std::string::npos) && 
      (seq.find(POLYA) == std::string::npos) && 
      (seq.find(POLYC) == std::string::npos) && 
      (seq.find(POLYG) == std::string::npos) && 
      (seq.find(POLYCG) == std::string::npos) && 
      (seq.find(POLYAT) == std::string::npos) && 
      (seq.find(POLYTC) == std::string::npos) && 
      (seq.find(POLYAG) == std::string::npos) && 
      (seq.find(POLYCA) == std::string::npos) && 
      (seq.find(POLYTG) == std::string::npos))
    return false;
  return true;
  
}

int overlapSize(const SeqLib::BamRecord& query, const SeqLib::BamRecordVector& subject) {

  // get the amount covered by subjet sequences
  typedef std::pair<int32_t, int32_t> AP;
  typedef std::pair<AP, std::string> AP_wseq;
  std::vector<AP_wseq> align_pos;
  std::string convention_seq = "";
  if (subject.size()) { 
    convention_seq = subject[0].Sequence(); // set the orientation convention
    for(auto& j : subject) {
      if (j.Sequence() == convention_seq)
	align_pos.push_back(AP_wseq(AP(j.AlignmentPosition(),j.AlignmentEndPosition()), j.Sequence()));
      else
	align_pos.push_back(AP_wseq(AP(j.AlignmentPositionReverse(),j.AlignmentEndPositionReverse()), j.Sequence()));
    }
  }
  
  bool same_orientation = query.Sequence() == convention_seq; 
  int al1 = same_orientation ? query.AlignmentPosition() : query.AlignmentPositionReverse();
  int al2 = same_orientation ? query.AlignmentEndPosition() : query.AlignmentEndPositionReverse();
  
  int max_overlap = 0;
  for (auto& k : align_pos) { // check each subject alignment to see if it overlaps query
    AP kk = k.first;
    max_overlap = std::max(max_overlap, std::max(0, std::min(kk.second, al2) - std::max(kk.first,al1)));
    //  out = true;
    //if (al1 >= kk.first && al2 <= kk.second)
    //  query_contained = true;
  }
  return max_overlap;
  
}

  void print(std::stringstream& s, std::ofstream& log, bool cerr) {
    log << s.str();
    if (cerr)
      std::cerr << s.str();
    s.str(std::string());
  }
  

  std::string runTimeString(int num_t_reads, int num_n_reads, int contig_counter, 
			    const SeqLib::GenomicRegion& region, const SeqLib::BamHeader& h, const svabaTimer& st,
			    const timespec& start) {
    
    std::stringstream ss;

    if ( (num_t_reads + num_n_reads) > 0 && !region.IsEmpty()) {
      std::string print1 = SeqLib::AddCommas<int>(region.pos1);
      std::string print2 = SeqLib::AddCommas<int>(region.pos2);
      char buffer[180];
      snprintf (buffer, 180, "***Ran %2s:%11s-%11s | T: %5d N: %5d C: %5d | ", 
	       h.IDtoName(region.chr).c_str(),
	       print1.c_str(),print2.c_str(),
	       (int)num_t_reads, (int)num_n_reads, 
	       (int)contig_counter);
      ss << std::string(buffer) << st << " | ";
#ifndef __APPLE__
      ss << SeqLib::displayRuntime(start);
#endif
    } else if (num_t_reads + num_n_reads > 0) {
      char buffer[180];
      snprintf (buffer, 180, "Ran Whole Genome | T: %5d N: %5d C: %5d | ", 
	       (int)num_t_reads, (int)num_n_reads, 
	       (int)contig_counter);
      
      ss << std::string(buffer) << st << " | ";
#ifndef __APPLE__
      ss << SeqLib::displayRuntime(start);
#endif
    }    
    return ss.str();
  }

  // just get a count of how many jobs to run. Useful for limiting threads. Also set the regions
  int countJobs(const std::string& regionFile, SeqLib::GRC &file_regions, SeqLib::GRC &run_regions, 
		const SeqLib::BamHeader& h, int chunk, int window_pad) {
    
    // open the region file if it exists
    if (!regionFile.empty() && fs::exists(regionFile)) {
      try {
	file_regions = SeqLib::GRC(regionFile, h);
      }
      catch (const std::exception &exc) {
	std::cerr
	  << "Found chromosome in region file " << regionFile
	  << " not in reference genome. Skipping.\n"
	  << "    Caught error: " << exc.what() << "\n";
      }
    }
    
    // parse as a samtools string eg 1:1,000,000-2,000,000
    else if (regionFile.find(":") != std::string::npos && regionFile.find("-") != std::string::npos) {
      file_regions.add(SeqLib::GenomicRegion(regionFile, h));
    }
    
    // it's a single chromosome
    else if (!regionFile.empty()) {
      SeqLib::GenomicRegion gr(regionFile, "1", "1", h);
      if (gr.chr == -1 || gr.chr >= h.NumSequences()) {
	std::cerr << "ERROR: Trying to match chromosome " << regionFile << " to one in header, but no match found" << std::endl;
	exit(EXIT_FAILURE);
      } else {
	gr.pos2 = h.GetSequenceLength(gr.chr); //get_()->target_len[gr.chr];
	file_regions.add(gr);
      }
    }

    // its empty, so cover whole genome
    else { 
      // add all chromosomes
      for (int i = 0; i < h.NumSequences(); i++) {
	file_regions.add(SeqLib::GenomicRegion(i, 1, h.GetSequenceLength(i))); //h.get()_->target_len[i]));
      }
    }
    
    // check if the mask was successful
    if (file_regions.size() == 0) {
      std::cerr << "ERROR: Cannot read region file: " << regionFile << 
	" or something wrong with bam header ('chr' prefix mismatch?)" << std::endl;
      exit(EXIT_FAILURE);
    }
    
    // divide it up
    if (chunk > 0) { // if <= 0, whole genome at once
      for (auto& r : file_regions) {
	SeqLib::GRC test(chunk, window_pad, r);
	run_regions.Concat(test);
      }
    }

    // now clear file regions, to signal that it was empty
    if (regionFile.empty())
      file_regions.clear(); 
    
    return run_regions.size();
      
  }
  

  std::string __bamOptParse(std::map<std::string, std::string>& obam, std::istringstream& arg, int sample_number, const std::string& prefix) {
    std::stringstream ss;
    std::string bam;
    arg >> bam;
    ss.str(std::string());
    ss << prefix << std::setw(3) << std::setfill('0') << sample_number;
    obam[ss.str()] = bam;
    return bam;
  }

  bool __openWriterBam(const SeqLib::BamHeader& h, const std::string& name, SeqLib::BamWriter& wbam) {
    
    if (!wbam.Open(name))
      return false;

    wbam.SetHeader(h);;
    if (!wbam.WriteHeader())
      return false;
    
    return true;
  }

  /*  bool __header_has_chr_prefix(bam_hdr_t * h) {
    for (int i = 0; i < h.NumSequences()) //->n_targets; ++i) 
      if (h->target_name[i] && std::string(h->target_name[i]).find("chr") != std::string::npos) 
	return true;
    return false;
  }
  */

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

  int weightedRandom(const std::vector<double>& cs) {
    
    // get a weighted random number
    size_t al = 0;
    double rand_val = rand() % 1000;
    while (al < cs.size()) {
      if (rand_val <= cs[al] * 1000) 
	return al;
      ++al;
    }
    return al;
  }

  std::vector<std::string> tokenize_delimited(const std::string& str, char delim) {
    std::vector<std::string> tokens;
    size_t start = 0;
    size_t end = str.find(delim);
    
    while (end != std::string::npos) {
      tokens.push_back(str.substr(start, end - start));
      start = end + 1;
      end = str.find(delim, start);
    }
    
    // Add the last token
    tokens.push_back(str.substr(start));
    
    return tokens;
  }
  
/// Check that the BAM header and BWA reference header line up.
void checkHeaderCompatibility(const SeqLib::BamHeader& bamHeader,
			      const SeqLib::BamHeader& refHeader,
			      SvabaLogger& logger)  {
  
  bool triggerExplain = false;

  if (bamHeader.NumSequences() != refHeader.NumSequences()) {
    triggerExplain = true;
    logger.log(
	       true, true,  //toerr, tolog
	       "!!!!!!!!!!! WARNING !!!!!!!!!!!\n",
	       "!!!!!! Number of sequences in BAM header mismatches reference\n",
	       "!!!!!! BAM: ", bamHeader.NumSequences(),
	       " -- Ref: ", refHeader.NumSequences()
	       );
  }

  
  if (bamHeader.NumSequences() != refHeader.NumSequences()) {
    triggerExplain = true;
    logger.log(true, true, //toerr, tolog
	       "!!!!!!!!!!! WARNING !!!!!!!!!!!\n",
	       "!!!!!! Number of sequences in BAM header mismatches reference\n"
	       "!!!!!! BAM: ", bamHeader.NumSequences(),
	       " -- Ref: ", refHeader.NumSequences());
  }
  
  const int n = std::min(bamHeader.NumSequences(), refHeader.NumSequences());
  for (int i = 0; i < n; ++i) {
    auto bamName = bamHeader.IDtoName(i);
    auto refName = refHeader.IDtoName(i);
    if (bamName != refName) {
      triggerExplain = true;
      logger.log(true, true,
		 "!!!!!! BAM sequence id ", i, ": \"", bamName, "\"",
		 " -- Ref sequence id ", i, ": \"", refName, "\"");
    }
  }
  
  if (triggerExplain) { 
    logger.log(true, true, std::string(100, '!'), "\n",
	       "!!! SvABA is being run with different reference genome than the reads were mapped to.\n",
	       "!!! This can cause a massive failure in variant detection!\n",
	       "!!! If you are *sure* that the two references are functionally equivalent (e.g. chr1 vs 1)\n",
	       "!!! and that the order of the chromosomes is equivalent between the two,\n",
	       "!!! you can override this error with option \"--override-reference-check\"\n",
	       std::string(100, '!'));
    std::exit(EXIT_FAILURE);
  }
}


  std::vector<std::pair<int, int>> find_repeats(std::string_view seq, size_t single_repeat_count = 0, size_t dinuc_repeat_count = 0) {
    std::vector<std::pair<int, int>> result;
    size_t n = seq.size();

    // Check for single-base repeats (e.g., "AAAA")
    if (single_repeat_count > 0 && n >= single_repeat_count) {
        for (size_t i = 0; i <= n - single_repeat_count; ++i) {
            char base = seq[i];
            size_t j = i + 1;
            while (j < n && seq[j] == base) ++j;

            size_t len = j - i;
            if (len >= single_repeat_count) {
                result.emplace_back(i, j - 1);
                i = j - 1;  // advance to end of this repeat block
            }
        }
    }

    // Check for dinucleotide repeats (e.g., "AGAGAG")
    if (dinuc_repeat_count > 0 && n >= 2 * dinuc_repeat_count) {
        for (size_t i = 0; i <= n - 2 * dinuc_repeat_count; ++i) {
            std::string_view unit = seq.substr(i, 2);
            size_t j = i + 2;
            size_t count = 1;

            while (j + 1 < n && seq.substr(j, 2) == unit) {
                ++count;
                j += 2;
            }

            if (count >= dinuc_repeat_count) {
                result.emplace_back(i, j - 1);
                i = j - 2;  // advance to end of this repeat block
            }
        }
    }

    return result;
  }
  
  std::vector<int> parsePLString(const std::string& pl_str) {
    std::vector<int> out;
    out.reserve(3);
    std::stringstream ss(pl_str);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
      try {
	out.push_back(std::stoi(tok));
      }
      catch (const std::exception& e) {
	throw std::runtime_error("Invalid integer in PL string: " + tok);
      }
    }
    if (out.size() != 3) {
      throw std::runtime_error("PL string must have exactly 3 comma-separated values");
    }
    return out;
  }
  
  // thresholds
  constexpr size_t HOMOPOLYMER_MIN  = 3;  // e.g. 16 the same base
  constexpr size_t DINUC_REPEAT_MIN = 2;   // e.g. 8 a nt motif
  
  SubstringList find_long_homopolymers(const std::string& s) {
    
    SubstringList results;
    
    if (s.size() < HOMOPOLYMER_MIN) return results;
    
    size_t run = 1;
    for (size_t i = 1; i < s.size(); ++i) {
        if (s[i] == s[i - 1]) {
            ++run;
        } else {
            if (run >= HOMOPOLYMER_MIN) {
                size_t end = i - 1;
                size_t start = end - run + 1;
                results.emplace_back(start, end, std::string(run, s[i - 1]));
            }
            run = 1;
        }
    }

    // Catch repeat at end of string
    if (run >= HOMOPOLYMER_MIN) {
      size_t end = s.size() - 1;
      size_t start = end - run + 1;
      results.emplace_back(start, end, std::string(run, s.back()));
    }
    
    return results;
  }
  
  SubstringList find_long_dinuc_repeats(const std::string& s) {
    
    SubstringList results;
    
    if (s.size() < 2 * DINUC_REPEAT_MIN) return results;
    
    for (size_t i = 0; i + 2 * DINUC_REPEAT_MIN <= s.size(); ++i) {
        std::string unit = s.substr(i, 2);
        size_t reps = 1;

        while (i + 2 * (reps + 1) <= s.size() &&
               s.substr(i + 2 * reps, 2) == unit) {
            ++reps;
        }

        if (reps >= DINUC_REPEAT_MIN) {
            size_t start = i;
            size_t end = i + 2 * reps - 1;

            std::string repeated_seq;
            for (size_t r = 0; r < reps; ++r)
                repeated_seq += unit;

            results.emplace_back(start, end, repeated_seq);

            i = end - 1; // Skip ahead past this repeat
        }
    }

    return results;
}

  
} // end namespace svabaUtils
