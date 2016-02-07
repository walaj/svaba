#include "SnowmanUtils.h"
#include <iomanip>

namespace SnowmanUtils {

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
    std::stringstream month;
    std::stringstream mdate;
    if ( (now->tm_mon+1) < 10)
      month << "0" << now->tm_mon+1;
    else 
      month << now->tm_mon+1;
    mdate << (now->tm_year + 1900) << month.str() <<  now->tm_mday;
    return mdate.str();
  }

  SnowTimer::SnowTimer() {
    s = {"r", "m", "as", "bw", "pp", "k"};
    for (auto& i : s)
      times[i] = 0;
    curr_clock = clock();
  }

  void SnowTimer::stop(const std::string& part) { 
    times[part] += (clock() - curr_clock); 
    curr_clock = clock();
  }

  void SnowTimer::start() { 
    curr_clock = clock(); 
  }

  std::ostream& operator<<(std::ostream &out, const SnowTimer st) {

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
    auto itp = st.times.find("pp");

    if (total_time)
      sprintf (buffer, "R: %2d%% M: %2d%% K: %2d%% A: %2d%% P: %2d%%", 
	       SnowTools::percentCalc<double>(itr->second, total_time),
	       SnowTools::percentCalc<double>(itm->second, total_time),
	       SnowTools::percentCalc<double>(itk->second, total_time),
	       SnowTools::percentCalc<double>(ita->second, total_time),
	       //SnowTools::percentCalc<double>(itb->second, total_time),
	       SnowTools::percentCalc<double>(itp->second, total_time));
    else
      sprintf (buffer, "NO TIME");
    out << std::string(buffer);
    return out;
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

int overlapSize(const SnowTools::BamRead& query, const SnowTools::BamReadVector& subject) {

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
			    const SnowTools::GenomicRegion& region, const bam_hdr_t * h, const SnowTimer& st,
			    const timespec& start) {
    
    std::stringstream ss;

    if ( (num_t_reads + num_n_reads) > 0 && !region.isEmpty()) {
      std::string print1 = SnowTools::AddCommas<int>(region.pos1);
      std::string print2 = SnowTools::AddCommas<int>(region.pos2);
      char buffer[180];
      sprintf (buffer, "Ran %2s:%11s-%11s | T: %5d N: %5d C: %5d | ", 
	       h->target_name[region.chr],
	       print1.c_str(),print2.c_str(),
	       (int)num_t_reads, (int)num_n_reads, 
	       (int)contig_counter);
      ss << std::string(buffer) << st << " | ";
      ss << SnowTools::displayRuntime(start);
      ss << std::endl;
    } else if (num_t_reads + num_n_reads > 0) {
      char buffer[180];
      sprintf (buffer, "Ran Whole Genome | T: %5d N: %5d C: %5d | ", 
	       (int)num_t_reads, (int)num_n_reads, 
	       (int)contig_counter);
      
      ss << std::string(buffer) << st << " | ";
      ss << SnowTools::displayRuntime(start);
      ss << std::endl;
      
    }
    
    return ss.str();
  }

  // just get a count of how many jobs to run. Useful for limiting threads. Also set the regions
  int countJobs(const std::string& regionFile, SnowTools::GRC &file_regions, SnowTools::GRC &run_regions, 
		bam_hdr_t * h, int chunk, int window_pad) {
    
    // open the region file if it exists
    bool rgfile = SnowTools::read_access_test(regionFile);
    if (rgfile)
      file_regions.regionFileToGRV(regionFile, 0);
    
    // parse as a samtools string eg 1:1,000,000-2,000,000
    else if (regionFile.find(":") != std::string::npos && regionFile.find("-") != std::string::npos)
      file_regions.add(SnowTools::GenomicRegion(regionFile, h));
    
    // it's a single chromosome
    else if (!regionFile.empty()) {
      SnowTools::GenomicRegion gr(regionFile, "1", "1", h);
      if (gr.chr == -1 || gr.chr >= h->n_targets) {
	std::cerr << "ERROR: Trying to match chromosome " << regionFile << " to one in header, but not match found" << std::endl;
	exit(EXIT_FAILURE);
      } else {
	gr.pos2 = h->target_len[gr.chr];
	file_regions.add(gr);
      }
    }
    // add all chromosomes
    else {
      for (int i = 0; i < h->n_targets; i++) {
	int region_id = bam_name2id(h, h->target_name[i]);
	
	//if (opt::verbose > 1)
	//  std::cerr << "chr id from header " << region_id << " name " << h->target_name[i] << " len " << h->target_len[i] << std::endl;
	
	if (region_id < 23) // don't add outsdie of 1-X
	  file_regions.add(SnowTools::GenomicRegion(region_id, 1, h->target_len[i]));
      }
    }
    
    //if (opt::verbose > 1)
    //for (auto& i : file_regions)
    //  std::cerr << "file regions " << i << std::endl;
    
    // check if the mask was successful
    if (file_regions.size() == 0) {
      std::cerr << "ERROR: Cannot read region file: " << regionFile << " or something wrong with tumor bam header" << std::endl;
      exit(EXIT_FAILURE);
    }
    
    // divide it up
    if (chunk > 0) // if <= 0, whole genome at once
      for (auto& r : file_regions) {
	SnowTools::GRC test(chunk, window_pad, r);
	run_regions.concat(test);
      }
    return run_regions.size();
    
  }
    

  std::string __bamOptParse(std::unordered_map<std::string, std::string>& obam, std::istringstream& arg, int sample_number, const std::string& prefix) {
    std::stringstream ss;
    std::string bam;
    arg >> bam;
    ss.str(std::string());
    ss << prefix << std::setw(3) << std::setfill('0') << sample_number;
    obam[ss.str()] = bam;
    return bam;
  }

  void __openWriterBam(const SnowTools::BamWalker& bwalker, const std::string& name, SnowTools::BamWalker& wbam) {
    bam_hdr_t * r2c_hdr = bam_hdr_dup(bwalker.header());
    wbam.SetWriteHeader(r2c_hdr);
    wbam.OpenWriteBam(name);
  }

  bool __header_has_chr_prefix(bam_hdr_t * h) {
    for (int i = 0; i < h->n_targets; ++i) 
      if (h->target_name[i] && std::string(h->target_name[i]).find("chr") != std::string::npos) 
	return true;
    return false;
  }

  void __open_bed(const std::string& f, SnowTools::GRC& b, bam_hdr_t* h) {
    //blacklist.add(SnowTools::GenomicRegion(1,33139671,33143258)); // really nasty region
    if (f.empty())
      return;
    b.regionFileToGRV(f, 0, h, __header_has_chr_prefix(h));
    b.createTreeMap();
  }
  
  faidx_t * __open_index_and_writer(const std::string& index, SnowTools::BWAWrapper * b, const std::string& wname, SnowTools::BamWalker& writer, faidx_t * findex) {

    if (!SnowTools::read_access_test(index))
      return findex;

    b->retrieveIndex(index);
    
    // open the bam for writing  
    writer.SetWriteHeader(b->HeaderFromIndex());
    writer.OpenWriteBam(wname); // open and write header

    findex  = fai_load(index.c_str());  // load the sequence
    return findex;
  }

  
}


