#ifndef ALIGNED_CONTIG
#define ALIGNED_CONTIG

#include "Util.h"
#include "Read2Contig.h"
#include "VCFRecord.h"
#include "GenomicRegion.h"
#include "seqan_tools.h"
#include <algorithm>
#include "SVBamReader.h"

using namespace std;

typedef vector<BamTools::BamAlignment> BAVec;
typedef vector<BamTools::CigarOp> CigarOpVec;

class AlignedContig;
typedef unordered_map<string, AlignedContig> ContigMap;

struct SBlat {

  // constructor that takes a line from a BLAT hit
  SBlat(string input) {
    
    string s_match, s_mismatch, s_repmatch, s_Ns, s_Qgap_count, s_Qgap_bases, 
      s_Tgap_count, s_Tgap_bases, s_strand, query_name, s_query_size, s_query_start,
      s_query_end, hit_name, s_hit_size, s_hit_start, s_hit_end, s_blockCount, s_blockSizes;

    istringstream iss(input);

    if (!(iss >> s_match >> s_mismatch >> s_repmatch >> s_Ns >> s_Qgap_count >>
	  s_Qgap_bases >> s_Tgap_count >> s_Tgap_bases >> s_strand >> query_name >>
	  s_query_size >> s_query_start >> s_query_end >> hit_name >> s_hit_size >> 
	  s_hit_start >> s_hit_end >> s_blockCount >> s_blockSizes))
      cerr << "Error in making SBlat object" << endl; 

    match = stoi(s_match);
    mismatch = stoi(s_mismatch);
    repmatch = stoi(s_repmatch);
    Ns = stoi(s_Ns);
    Qgap_count = stoi(s_Qgap_count);	  
    Qgap_bases = stoi(s_Qgap_bases);
    Tgap_count = stoi(s_Tgap_count);
    Tgap_bases = stoi(s_Tgap_bases);
    query_size = stoi(s_query_size);
    query_start = stoi(s_query_start);
    query_end   = stoi(s_query_end);
    hit_size = stoi(s_hit_size);
    hit_start = stoi(s_hit_start);
    hit_end = stoi(s_hit_end);
    blockCount = stoi(s_blockCount);
    blockSizes = stoi(s_blockSizes);
  }

  int match, mismatch, repmatch, Ns, Qgap_count, Qgap_bases, Tgap_count, Tgap_bases;
  char strand;
  string query_name;
  int query_size, query_start, query_end;
  string hit_name;
  int hit_size, hit_start, hit_end;
  int blockCount, blockSizes; 

  // define how these are to be sorted. Sort by biggest match first
  bool operator < (const SBlat& b) const { return (match > b.match); }

};

typedef vector<SBlat> SBlatVec;

struct RepeatMasker {

  RepeatMasker(string input) {

    string s_sw, s_perc_div, s_perc_del, s_perc_ins, s_contig, 
      s_query_begin, s_query_end, s_query_left, dum, s_repeat_class, s_repeat_name, 
      s_repeat_begin, s_repeat_end, s_repeat_left, s_id, dum2;

    istringstream iss(input);
    if (!(iss >> s_sw >> s_perc_div >> s_perc_del >> s_perc_ins >> s_contig >> s_query_begin 
	  >> s_query_end >> s_query_left >> dum >> s_repeat_name >> s_repeat_class >> 
	  s_repeat_begin >> s_repeat_end >> s_repeat_left >> s_id)) 
      cerr << "Error in making RepeatMasker object" << endl; 

    sw = stoi(s_sw);
    perc_div = stod(s_perc_div);
    perc_del = stod(s_perc_del);
    perc_ins = stod(s_perc_ins);
    query_begin = stoi(s_query_begin);
    query_end = stoi(s_query_end);
    contig = s_contig;
    repeat_class = s_repeat_class;
    repeat_name = s_repeat_name;
    
  }

  //RepeatMasker(int tsw, string tcontig, string trep_name, string trep_class, 
  //	       double tdiv, double tdel, double tins, int beg, int end) :
  // sw(tsw), contig(tcontig), repeat_name(trep_name), repeat_class(trep_class), perc_div(tdiv), perc_del(tdel), perc_ins(tins), 
  // query_begin(beg), query_end(end) {}

  int sw;
  string contig;
  string repeat_name;
  string repeat_class;
  
  double perc_div;
  double perc_del;
  double perc_ins;
  
  int query_begin;
  int query_end;

  // define how these are to be sorted. Sort by biggest sw first
  bool operator < (const RepeatMasker& rep) const { return (sw > rep.sw); }

  // print it out
  string toString() {
    stringstream ss;
    ss << "SW: " << sw << " Contig: " << contig << " RepeatName: " << repeat_name << 
      " RepeatClass: " << repeat_class << " %Div: " << perc_div << " %Del: " << perc_del << 
      " %Ins: " << perc_ins << " QueryBegin: " << query_begin << " QueryEnd: " << query_end;
    return ss.str();
  }

};

typedef vector<RepeatMasker> RepeatMaskerVec;

/*struct Window {

  Window(int trefID, int tpos1, int tpos2) : refID(trefID), pos1(tpos1), pos2(tpos2) {}
  Window() {};
  ~Window() {}
  int refID;
  int pos1;
  int pos2;

  string toString() const;

  };*/

struct CAlignment {

  BamTools::BamAlignment align;
  CigarOpVec cigar;
  string cigstring; 
  int break1;
  int break2;
  int gbreak1;
  int gbreak2;
  unsigned start;

  unsigned nsplit1 = 0;
  unsigned tsplit1 = 0;
  unsigned nsplit2 = 0;
  unsigned tsplit2 = 0;

  size_t ndisc = 0;
  
  bool ca_local = false; // is an alignmend to anchor window

  CAlignment(BamTools::BamAlignment talign, CigarOpVec tcigar, string tcigstring) : align(talign), cigar(tcigar), cigstring(tcigstring) {}

  // define how to sort these 
  bool operator < (const CAlignment& str) const { return (start < str.start); }

};


// define a way to order the contigs by start
struct AlignmentOrdering {
  inline bool operator() (const CAlignment& struct1, const CAlignment& struct2) {
    return (struct1.start < struct2.start);
  }
};

struct BreakPoint {

  bool isBest = false; // marked for best breakpoint when mulitple are redundant
 
  size_t disco_tum = 0;
  size_t disco_norm = 0;

  string idcommon = "";
  string pairid = "";
  
  bool discovar = false;

  DiscordantCluster dc;

  unsigned pos1;
  unsigned pos2;

  unsigned cpos1;  
  unsigned cpos2;

  unsigned refID1;
  unsigned refID2;

  string seq;

  string cname;

  string insertion;
  string homology;

  string id1;
  string id2;
  int matchlen1 = 0;
  int matchlen2 = 0;
  
  char strand1;
  char strand2;

  bool isSomatic = false;
  bool isGermline = false;

  unsigned mapq1; 
  unsigned mapq2; 

  unsigned tsplit1 = 0;
  unsigned tsplit2 = 0;

  unsigned nsplit1 = 0;
  unsigned nsplit2 = 0;

  size_t nsplit = 0;
  size_t tsplit = 0;

  unsigned tall = 0;
  unsigned nall = 0; 

  int nm1 = 0;
  int nm2 = 0;

  unsigned num_dups = 0;
   
  //Window window;
  GenomicRegion window;

  int span;

  unsigned num_align = 0;

  bool part_of_local = false;

  bool local1 = false;
  bool local2 = false;

  string evidence = "";
  string confidence = "";

  BreakPoint(DiscordantCluster tdc);
  BreakPoint() {}

  static string BreakPointHeader();

  string toString() const; 
 
  bool sameBreak(BreakPoint &bp) const;

  void order();

  // return whether a bp is good to move on
  bool isGoodSomatic(int mapq, size_t tsplit_cutoff, size_t nsplit_cutoff) const;

  bool hasDiscordant() const;

  // return whether a bp is good to move on
  bool isGoodGermline(int mapq, size_t allsplit) const;

  // define how to sort these 
  bool operator < (const BreakPoint& bp) const { 
    return bp.refID1 > refID1 || 
       (bp.refID1 == refID1 && bp.pos1 > pos1) || // low pos is first
       (bp.refID1 == refID1 && bp.pos1 == pos1 && bp.pos2 == pos2 && nsplit1 > bp.nsplit1) || // if same, check nsplit
       (bp.refID1 == refID1 && bp.pos1 == pos1 && bp.pos2 == pos2 && nsplit1 == bp.nsplit1 && tsplit1 > bp.tsplit1); // if also same, check tsplit
  }
  friend ostream& operator<<(std::ostream& out, const BreakPoint& bp) { out << bp.toString(); return out; }

  // print to file
  void printToFile(ofstream &of, ContigMap * contigs);

  // print to VCF
  string printToVCF(int split_cut, int mapq_cut, int uniq) const;

  // return VCFRecord
  VCFRecordVector getVCFRecord(int uniq, const BamAlignmentVector &bamreads) const;

};

typedef vector<CAlignment> AlignVec;
typedef vector<BreakPoint> BPVec;
typedef unordered_map<string, BreakPoint> BPMap;

class AlignedContig {

 public:  

  // constructor taking in a BAM record from the contig BAM
  AlignedContig(const BamTools::BamAlignment align); 
  AlignedContig() {}
  ~AlignedContig() {}

  BPVec m_breaks;  
  BreakPoint m_farbreak;
  BreakPoint m_farbreak_filt; // global breakpoint, but remove bad fragments (e.g. 60, 60, 0 mapqs, keep 60-60 connection)
  vector<BamAlignment> m_bamreads;
  AlignVec m_align;

  // add a new contig alignment
  void addAlignment(const BamTools::BamAlignment align);

  void fillExtraReads();

  void addDiscordantCluster(DiscordantCluster dc) { m_dc.push_back(dc); } 

  int numDiscordantClusters() const { return m_dc.size(); }

  string printDiscordantClusters() const;

  // return the name of the contigs
  string getContigName() const { return m_align[0].align.Name; }
 
  // return the number of alignments
  int getNumPrimaryAlign() const { return m_align.size(); }
  int getNumSecondaryAlign() const { return m_align_second.size(); }

  int getMaxMapq() const { return m_maxmapq; }
  int getMinMapq() const { return m_minmapq; }

  bool hasLocal() const { return m_local; }
  R2CVec getReads() const { return  m_reads; }

  //
  int getNumReads() const { return m_reads.size(); };

  int getNumReads(const char type) const { 
    int countr = 0;
    for (R2CVec::const_iterator it = m_reads.begin(); it != m_reads.end(); it++) 
      if (it->rname.at(0) == type)
	countr++;
    return countr; 
  };

  // sort the read2contigs
  void sortReads();

  // make work function for getting PER-ALIGNMENT breaks from BWA-MEM
  void setBreaks(CAlignment &align);

  void addRead2Contig(Read2Contig rc);

  void splitCoverage();

  string printForR() const;

  // set whether this alignment has a somatic breakpoint
  //void setSomatic();

  void printAlignments(ofstream &ostream) const;
  void printContigFasta(ofstream &ostream) const;

  bool isGermline() const { return m_germline; };
  bool isSomatic() const { return m_somatic; };

  bool intersect(const AlignedContig *al) const; 

  // flips the cigar if the contig is aligned to the opposite strand
  CigarOpVec orientCigar(const BamTools::BamAlignment align);

  // converts the BamTools cigar format to string
  string cigarToString(CigarOpVec cig);

  // parses the contig file name to determine where the anchor window was
  //void setWindow(const string s);

  // 
  //Window getWindow() const { return m_window; }

  // find the breakpoint pairs by looping through ALL the alignments
  void getBreakPairs();
   
  // get the break pairs, but update incase things have changed
  BPVec getBreaks() const { 
    return m_breaks; 
  }

  // add masked seq
  void addMaskedSeq(string seq) { m_masked_seq  = seq; }

  // add repeat masker
  void addRepeatMasker(RepeatMasker rp) { 
    m_rep_vec.push_back(rp); 
    sort(m_rep_vec.begin(), m_rep_vec.end());
  }

  // add SBlat alignment
  void addSBlat(SBlat b) { 
    m_blat_vec.push_back(b); 
    sort(m_blat_vec.begin(), m_blat_vec.end());
  }
  
  // run all the algorithms for updating the breakpoints, in case other parts changed
  void updateBreakpointData(bool skip_realign, bool no_r2c_matched) {

    // realign the reads to the contigs
    if (!skip_realign) 
      realignReads();
    else
      readR2Creads();

    // sort the reads for better visualization
    if (!no_r2c_matched)
      sortReads(); 

    // get the split coverage
    splitCoverage();

    // get the break pairs
    getBreakPairs();

    //debug
    if (m_breaks.size() == 0)
      cerr << "m_align size: "<< m_align.size() << endl;

    // if we don't need to write the output, clear most of the read data
    if (no_r2c_matched) {
      for (BamAlignmentVector::iterator it = m_bamreads.begin(); it != m_bamreads.end(); it++) {
	it->QueryBases = "";
	it->Qualities = "";
	it->RemoveTag("JW");
	it->RemoveTag("TS");
      }
    }

  }

  void readR2Creads();

  BreakPoint getGlobalBreak() const { return m_farbreak; }

  // doing a more stringent alignment of reads to the contig
  // this is to remove normals that ruin somatic calls
  void realignReads();

  string getName() const { return m_align[0].align.Name; }

  size_t getNumTmpAlign() const { return m_tmpalign.size(); }
  
  AlignVec getAlignments() const { return m_align; }

  void setSomatic(const bool somatic) { m_somatic = somatic; }

  void setGermline(const bool germline) { m_germline = germline; }

  void addTmpAlignment(const BamTools::BamAlignment align) { m_tmpalign.push_back(align); }

  void addBams(string tum, string norm, string pan, string r2c) { 
    tbam = tum;
    nbam = norm;
    pbam = pan;
    rbam = r2c;
  }

  vector<BamAlignment> getBamReads() const { return m_bamreads; }

  void settleContigs();

  double SWalign(Read2Contig &r, TSequence &contig, bool revcomp);
  double SWalign(TSequence &contig, bool revcomp, int32_t &pos, string &rseq, int32_t &score);

  string getSequence() const { return m_align[0].align.QueryBases; }
  
 private:

  AlignVec m_align_second;
  R2CVec m_reads;

  bool m_local = false;
  bool m_somatic = false;
  bool m_germline = false;
  int m_maxmapq = 0;
  int m_minmapq = 61;
  GenomicRegion m_window;
  string m_seq = "";
  string m_masked_seq = "";
  RepeatMaskerVec m_rep_vec;
  SBlatVec m_blat_vec;
  int mapq_threshold = 60;
  // store the raw alignments. Move into AlignVec if keeping the contig
  BAVec m_tmpalign; 

  // the bam files that
  string tbam, nbam, pbam, rbam;
  //SVBamReader treader, nreader, preader;

  vector<DiscordantCluster> m_dc;

};



#endif


