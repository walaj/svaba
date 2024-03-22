#include "svabaUtils.h"

#include <iomanip>

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
    s = {"r", "m", "as", "bw", "pp", "t", "k"};
    for (auto& i : s)
      times[i] = 0;
    curr_clock = clock();
  }

  void svabaTimer::stop(const std::string& part) { 
    times[part] += (clock() - curr_clock); 
    curr_clock = clock();
  }

  void svabaTimer::start() { 
    curr_clock = clock(); 
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
      snprintf (buffer, 180, "Ran %2s:%11s-%11s | T: %5d N: %5d C: %5d | ", 
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
    bool rgfile = SeqLib::read_access_test(regionFile);
    if (rgfile) {
      try {
	file_regions = SeqLib::GRC(regionFile, h);
      } catch (const std::exception &exc) {
	std::cerr << "Found chromosome in region file " << regionFile << " not in reference genome. Skipping." << std::endl;
	std::cerr << "     Caught error: " << exc.what() << std::endl;
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
    else { 
      // add all chromosomes
      for (int i = 0; i < h.NumSequences(); i++) {
	//if (i < 23) // don't add outsdie of 1-X
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

  void __open_bed(const std::string& f, SeqLib::GRC& b, const SeqLib::BamHeader& h) {
    if (f.empty())
      return;
    b = SeqLib::GRC(f, h);
    b.CreateTreeMap();
  }
  
  bool __open_index_and_writer(const std::string& index, SeqLib::BWAWrapper * b, const std::string& wname, SeqLib::BamWriter& writer, SeqLib::RefGenome *& r, SeqLib::BamHeader& bwa_header) {
    
    // load the BWA index
    if (!b->LoadIndex(index))
      return false;

    // load the same index, but for querying seq from ref
    if (!r->LoadIndex(index))
      return false;

    // get the dictionary from reference
    bwa_header = b->HeaderFromIndex();
    
    // open the bam for writing  
    writer.SetHeader(bwa_header);
    if (!writer.Open(wname)) // open and write header
      return false;
    writer.WriteHeader();

    return true;
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
  
  
}
