#pragma once

#include <cstdlib>
#include <string>
#include <map>
#include <tuple>
#include <unordered_map>

#include "DiscordantCluster.h"
#include "SeqLib/BamRecord.h"
#include "svabaRead.h"

enum class LocalAlignment {
  NOTSET,          
  NONVAR_LOCAL_REALIGNMENT,   // contig aligns to local region without variant
  FROM_LOCAL_REGION,          // at least one end of contig bp maps to assembled region, this is good
  FROM_DISTANT_REGION          // neither end maps to local assembly region, this is bad
};

enum class SVType {
  NOTSET,
  TSI_LOCAL,          // TSI_LOCAL
  TSI_GLOBAL,         // TSI_GLOBAL
  ASSMB,              // SIMPLE ASSEMBLY
  ASDIS,
  DSCRD,         // DISCORDANT ONLY
  INDEL               // INDEL
};

// forward declares
class STCoverage;
namespace SeqLib {
  class RefGenome;
  class BamHeader;
}
class AlignmentFragment;
class AlignedContig;

using SeqLib::BamHeader;
using SeqLib::BamRecord;
using SeqLib::RefGenome;
using SeqLib::BamRecordPtr;
using SeqLib::GenomicRegion;
using SeqLib::GRC;

struct BreakPoint;
struct BreakEnd;

typedef std::vector<std::string> HashVector;
typedef std::vector<BreakPoint> BPVec;
typedef std::shared_ptr<BreakPoint> BreakPointPtr;
typedef std::vector<BreakPointPtr> BreakPointPtrVector;

typedef std::unordered_set<std::string> ReadNameSet;

/*struct ReducedBreakEnd {
  
  
  ReducedBreakEnd(): mapq(0), sub_n(0), nm(0) {}
  
  ReducedBreakEnd(const SeqLib::GenomicRegion& g, int mq);

  friend std::ostream& operator<<(std::ostream& os, const ReducedBreakEnd& rbe);
  
  SeqLib::GenomicRegion gr;
  uint32_t mapq:8, sub_n:8, nm:16;  
};
*/

 struct BreakEnd {
   
   BreakEnd() = default;

   void transferContigAlignmentData(const AlignmentFragment* f,
				    bool isleft);
   
   void setLocal(const GenomicRegion& window);

   std::string printSimple(const BamHeader& h) const;

   std::string hash(int offset) const;

   std::string id;
   GenomicRegion gr;
   
   // contig level informaiton set by AlignmentFragment
   int mapq = -1;
   int cpos = -1;
   int nm = -1;
   int matchlen = -1;
   int sub_n = -1;
   double as_frac= 0;
   int as = -1;
   LocalAlignment local = LocalAlignment::NOTSET;

   // for high-confidence reads
   std::unordered_map<std::string, int> split;  
   //std::unordered_map<std::string, double> af;

   //   friend std::ostream& operator<<(std::ostream& out, const BreakEnd& b);
 };

class BreakPoint {

public:
  
  // for discordant clusters
  BreakPoint(DiscordantCluster& tdc,
	     DiscordantClusterMap& dmap, 
	     const GenomicRegion& region,
	     const SvabaSharedConfig* _sc);

  // for rearrangements (including complex)
  BreakPoint(const SvabaSharedConfig* _sc,
	     const AlignmentFragment* left,
	     const AlignmentFragment* right,
	     const AlignedContig* alc);

  // for indels
  BreakPoint(const AlignmentFragment* f,
	     const int idx, // index for which indel on this f alignment we are making
	     const SvabaSharedConfig* _sc);  

  // for readback from file
  BreakPoint(const std::string &line,
	     const BamHeader& h,
	     const SvabaSharedConfig* _sc);

  // enforce one of above modalities
  BreakPoint() = delete;
  
  // Disable copy
  BreakPoint(const BreakPoint&) = delete;
  BreakPoint& operator=(const BreakPoint&) = delete;
  
  // Disable move
  BreakPoint(BreakPoint&&) = delete;
  BreakPoint& operator=(BreakPoint&&) = delete;
  
  static std::string header() {
    return std::string(
		       "#chr1\tpos1\tstrand1\tchr2\tpos2\tstrand2\tref\talt\t"
		       "span\tmapq1\tmapq2\tnm1\tnm2\tdmq1\tdmq2\tsplit\t"
		       "cigar\talt\tcov\tdcn\tdct\tsub_n1\tsub_n2\thomol\tinsert\t"
		       "contig_and_region\tnumalign\tconf\ttype\tqual\t2ndary\t"
		       "somscore\tsomlod\tmaxlod\tdbsnp"
		       );
  }
  
  double somatic_score = -1;
  
  // LogOdds that variant not in normal   
  double somatic_lod = -1; 
  
  std::string seq, cname, rs,
    insertion, homology, repeat_seq,
    confidence, ref, alt;
  
  // the evidence per break-end
  BreakEnd b1, b2;
  
  int imprecise = -1; //false;
  
  // number of matched bases to left and right 
  // dont want indel matches where there is not confidence alignment
  // because indel is too close to end
  int left_match = -1, right_match = -1;
  
  // reads spanning this breakpoint (these are read-to-genome, but have an r2c)
  svabaReadPtrVector reads;
  
  // discordant reads supporting this aseembly bp
  DiscordantCluster dc;
  
  int quality = -1;
  int secondary = -1;
  int pass = -1; //false;
  int pon = 0;
  int num_align = 0;
  SVType svtype = SVType::NOTSET;
  LocalAlignment local = LocalAlignment::NOTSET;

  const SvabaSharedConfig* sc;
  
  double error_rate = 1e-4;
  
  // keep track of how much of contig is covered by split
  std::pair<int,int> split_cov_bounds = std::pair<int, int>(1e5, -1); // dummy to extreme opposite vals

  GenomicRegion BreakEndAsGenomicRegionLeft() const;
  
  GenomicRegion BreakEndAsGenomicRegionRight() const;

  std::string printDeletionMarksForAlignmentsFile() const;
  
  bool isIndel() const;
  
  HashVector getBreakEndHashes();
  
  void __combine_alleles();
  
  void __rep(int rep_num, std::string& rseq, bool fwd = true);
  
   void setLocal(const GenomicRegion& window);
  
  void score_somatic(double NODBCUTOFF, double DBCUTOFF);
  
  void addCovs(const std::unordered_map<std::string, STCoverage*>& covs);

   /** Retrieve the reference sequence at a breakpoint and determine if 
    * it lands on a repeat */
   void repeatFilter();

   std::string printSimple(const BamHeader& h) const;

   void CombineWithDiscordantClusterMap(DiscordantClusterMap& dmap);
   
   /*! @function determine if the breakpoint has split read support
    * @param reference to a vector of read smart pointers that have been aligned to a contig
    * @discussion Note: will cause an error if the AL tag not filled in for the reads. 
    * The AL tag is filled in by AlignedContig::alignReadsToContigs.
    */
   void splitCoverage(svabaReadPtrVector& bav);
   
   /*! Determines if the BreakPoint overlays a blacklisted region. If 
    * and overlap is found, sets the blacklist bool to true.
    *
    * Note that currently this only is set for the pos1 of indels.
    * If the BreakPoint object is not an indel, no action is taken. 
    * @param grm An interval tree map created from a BED file containing blacklist regions
    */
   void checkBlacklist(GRC &grv);
   
   /*! Score a breakpoint with a QUAL score, and as somatic or germline
    */
   void scoreBreakpoint(double LOD_CUTOFF, double LOD_CUTOFF_DBSNP, double LOD_CUTOFF_SOMATIC, double LOD_CUTOFF_SOMATIC_DBSNP, double scale_errors, int min_dscrd_size);
   
   /*! Compute the allelic fraction (tumor and normal) for this BreakPoint.
    *
   * The allelic fraction is computed by taking the base-pair level coverage
   * as the denominator, and the max of number of split reads and number of 
   * cigar supporting reads as the numerator. Note that because, theoretically
   * but rarely, the number of split reads could be > 0 while the bp-level coverage
   * at a variant could be exactly zero. This is because unmapped reads could be called split
   * reads but are not counted in the coverage calculation. In such a case, the allelic fraction is
   * set to -1. By the same argument, the allelic fraction could rarely be > 1.
   * @param t_cov Base-pair level Coverage object, with coverage for all reads from Tumor bam(s).
   * @param n_cov Base-pair level Coverage object, with coverage for all reads from Normal bam(s).
   */
   //void addAllelicFraction(STCoverage * t_cov, STCoverage * n_cov);
  
  /*! @function get the span of the breakpoints (in bp). -1 for interchrom
   * @return int distance between breakpoints
   */
   int getSpan() const;

   /*! @function get a unique string representation of this breakpoint.
    * Format for indel is chr_breakpos_type (eg. 0_134134_I)
    * @return string with breakpoint info
    */
  std::string getHashString() const;
  
  bool hasMinimal() const;
  
  bool sameBreak(BreakPoint &bp) const;
  
  void order();
  
  bool isEmpty() const { return (b1.gr.pos1 == 0 && b2.gr.pos1 == 0); }
  
  std::string toFileString(const BamHeader& header);
  
  bool hasDiscordant() const;
  
  //bool operator==(const BreakPoint& bp) const;
  
  // define how to sort these
  bool operator<(const BreakPoint& o) const {
    return std::tie(
		    b1.gr,      // first key
		    b2.gr,      // then
		    n.split,
		    t.split,
		    dc.ncount,
		    dc.tcount,
		    cname       // last key
		    ) < std::tie(
				 o.b1.gr,
				 o.b2.gr,
				 o.n.split,
				 o.t.split,
				 o.dc.ncount,
				 o.dc.tcount,
				 o.cname
				 );
  }
  
  void score_dscrd(int min_dscrd_size);
  void score_assembly_only();
  void score_assembly_dscrd();
  void score_indel(double LOD_CUTOFF, double LOD_CUTOFF_DBSNP);
  void set_homologies_insertions();
  bool valid() const;
  
  void setRefAlt(const RefGenome* main_rg,
		 const BamHeader& header);
  
  
  
  private:
  
  struct SampleInfo {

    
    // no default construction:
    SampleInfo() = default;

    int split = 0;
    int cigar = 0;
    int alt = 0;
    int cov = 0;
    int disc = 0;
    
    // genotype info
    double GQ = 0;
    double NH_GQ = 0; // GQ of 0/0. Higher is more likely to be not hom ref
    std::string PL;
    std::string genotype;
    std::vector<double> genotype_likelihoods = {0,0,0};
    
    double af = 0;
    double error_rate = 1e-4;
    
    double LO = 0; // log odds of variant vs error
    double SLO = 0; // MAPQ scaled log odds of variant vs error
    double LO_n = 0; // log odds of variant at af=0.5 vs ref (af=0) with errors
    
    // UniqueName set of supporting reads
    // not for read tracking but actually needed for breakpoint scoring
    ReadNameSet supporting_reads; 
    
    static double LogLikelihood(double ref, double alt,
				double f, double e_fwd,
				double e_back);
    
    void modelSelection(double err, int readlen);
    
    double __genotype_likelihoods(int g, double er, int alt, int cov);
     
    std::string toFileString(SVType svtype) const;
     
    void FillFromString(const std::string& s, SVType svtype);
     
     void UpdateAltCounts();
     
   };

    // Friend declarations so these free functions can see the private nested type
  friend std::ostream& operator<<(
    std::ostream& out,
    const BreakPoint::SampleInfo& a
  );
  friend BreakPoint::SampleInfo operator+(
    const BreakPoint::SampleInfo& a1,
    const BreakPoint::SampleInfo& a2
  );
  
public:

  

  //tumor allele, normal allele, all allele
  SampleInfo t, n, a;
  
  // ordered to keep in alphabetical order by prefix (e.g. n001)
  std::map<std::string, SampleInfo> allele; 
  
};

 // struct ReducedDiscordantCluster {
 //   uint32_t mapq1:8, mapq2:8, tcount:8, ncount:8;
 // };
 
 // struct ReducedBreakPoint {

 //   // some helper functions
 //   char* __string_alloc2char(const std::string& str, char * p) {
 //     if (!str.empty() && str != "x") {
 //       p = (char*)malloc(str.length() + 1);
 //       strcpy(p, str.c_str());
 //       return p;
 //     } else {
 //       return nullptr;
 //     }
 //   }
   
 //   inline void smart_check_free(char * p) {
 //     if (p)
 //       free(p);
 //   }

 //   int getSpan() const {
 //     if (indel && !insertion) // deletion
 //       return (abs((int)b1.gr.pos1 - (int)b2.gr.pos1) - 1);
 //     if (indel) // insertion
 //       return (strlen(insertion)); // insertion
 //     if (b1.gr.chr == b2.gr.chr)
 //       return abs((int)b1.gr.pos1-(int)b2.gr.pos1);
 //     else
 //       return -1;

 //   }

 //   // define how to sort these  
 //   bool operator<(const ReducedBreakPoint& bp) const;

 //   // print it with the correct chromsome string
 //   std::string print(const BreakPoint& b, const SeqLib::BamHeader& h) const;

 //   ReducedBreakPoint() {}
 //   ~ReducedBreakPoint() {
 //     smart_check_free(ref);
 //     smart_check_free(alt);
 //     smart_check_free(cname);
 //     smart_check_free(homology);
 //     smart_check_free(insertion);
 //     smart_check_free(evidence);
 //     smart_check_free(confidence);
 //     smart_check_free(repeat);
 //   }
 //   ReducedBreakPoint(const std::string &line, const SeqLib::BamHeader& h);

 //   char * ref;
 //   char * alt;
 //   char * cname;
 //   char * evidence;
 //   char * confidence;
 //   char * insertion;
 //   char * homology;
 //   char * repeat;

 //   std::string read_names, bxtable;

 //   std::vector<std::string> format_s;

 //   ReducedBreakEnd b1, b2;
 //   double somatic_score = 0;
 //   double somatic_lod = 0; // LogOdds that variant not in normal
 //   double true_lod = 0;

 //   //uint32_t nsplit:8, tsplit:8, 
 //   //uint32_t tcov_support:8, ncov_support:8, tcov:8, ncov:8;
 //   uint32_t cov:16, af_n:7, num_align:5, secondary:1, dbsnp:1, pass:1, blacklist:1, indel:1, imprecise:1;
 //   uint32_t tcigar:8, ncigar:8, dummy:8, af_t:8; 
 //   float quality;
 //   uint8_t pon;

 //   ReducedDiscordantCluster dc;

 // };

