#pragma once

#include <cstdlib>
#include <string>
#include <map>
#include <unordered_map>

#include "DiscordantCluster.h"
#include "svabaRead.h"

// forward declares
class STCoverage;
namespace SeqLib {
  class RefGenome;
  class BamHeader;
}
struct BreakPoint;

 struct SampleInfo {

   bool indel = false;

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

   int readlen = 0;

   double af = 0;
   double error_rate = 1e-4;

   double LO = 0; // log odds of variant vs error
   double SLO = 0; // MAPQ scaled log odds of variant vs error
   double LO_n = 0; // log odds of variant at af=0.5 vs ref (af=0) with errors

   std::set<std::string> supporting_reads; // holds SR tags (not qnames)

   friend std::ostream& operator<<(std::ostream& out, const SampleInfo& a);

   friend SampleInfo operator+(const SampleInfo& a1, const SampleInfo& a2);

   double __log_likelihood(double ref, double alt, double f, double e);

   void modelSelection(double err);

   double __genotype_likelihoods(int g, double er, int alt, int cov);

   std::string toFileString() const;

   void fromString(const std::string& s);

   void adjust_alt_counts();
   
 };

   
struct ReducedBreakEnd {


  ReducedBreakEnd(): mapq(0), sub_n(0), nm(0) {}
  
  ReducedBreakEnd(const SeqLib::GenomicRegion& g, int mq, const std::string & chr_n);

  friend std::ostream& operator<<(std::ostream& os, const ReducedBreakEnd& rbe);
  
  std::string chr_name;
  SeqLib::GenomicRegion gr;
  uint32_t mapq:8, sub_n:8, nm:16;  
};


 struct BreakEnd {
   
   BreakEnd() { mapq = 0; sub_n = 0; nm = 0; }

   BreakEnd(const SeqLib::GenomicRegion& g, int mq, const std::string & chr_n);
   
   BreakEnd(const SeqLib::BamRecord& b);
   
   void checkLocal(const SeqLib::GenomicRegion& window);

   std::string print(const SeqLib::BamHeader& h) const;

   std::string hash(int offset = 0) const;

   std::string id;
   std::string chr_name;
   SeqLib::GenomicRegion gr;

   int mapq = -1;
   int cpos = -1;
   int nm = -1;
   int matchlen = -1;
   int simple = 0;

   std::unordered_map<std::string, int> split;  // for high-confidence reads
   std::unordered_map<std::string, int> splitI; // for informative reads
   std::unordered_map<std::string, double> af;

   int sub_n = -1;
   double as_frac= 0;
   bool local = false;

   //   friend std::ostream& operator<<(std::ostream& out, const BreakEnd& b);
 };

 struct BreakPoint {
   
   static std::string header() { 
     return "#chr1\tpos1\tstrand1\tchr2\tpos2\tstrand2\tref\talt\tspan\tmapq1\tmapq2\tnm1\tnm2\tdisc_mapq1\tdisc_mapq2\tsplit\tcigar\talt\tcov\tsub_n1\tsub_n2\thomology\tinsertion\tcontig\tnumalign\tconfidence\tevidence\tquality\tsecondary_alignment\tsomatic_score\tsomatic_lod\ttlod\tpon_samples\trepeat_seq\tgraylist\tDBSNP\treads\tbxtags"; 
   }

   double somatic_score = 0;
   double somatic_lod = 0; // LogOdds that variant not in normal

   int readlen = 0;

   int aligned_covered = 0;
   
   std::string seq, cname, rs, insertion, homology, repeat_seq, evidence, confidence, ref, alt, read_names, bxtable;   

   // count of unique bx tags
   size_t bx_count = 0;

   // the evidence per break-end
   BreakEnd b1, b2;

   // number of matched bases to left and right 
   // dont want indel matches where there is not confidence alignment
   // because indel is too close to end
   int left_match = 0, right_match = 0;
   
   SampleInfo t, n, a;

   // reads spanning this breakpoint
   svabaReadVector reads;

   // store if it has a non-clipped local alignment
   bool has_local_alignment = false;

   //int t_reads = 0, n_reads = 0;

   // discordant reads supporting this aseembly bp
   DiscordantCluster dc;
   
   int quality = 0;

   // total coverage at that position
   std::map<std::string, SampleInfo> allele; // ordered to keep in alphabetical order by prefix (e.g. n001)

   bool secondary = false;

   //std::unordered_set<std::string> split_reads, qnames;
   
   int pon = 0;
   int num_align = 0;

   bool complex = false;
   bool complex_local = false; // a local piece (e.g. AB) of a complex break (ABCD)
   
   bool isindel = false;
   bool blacklist = false;

   double error_rate = 1e-4;

   // keep track of how much of contig is covered by split
   std::pair<int,int> split_cov_bounds = std::pair<int, int>(1e5, -1); // dummy to extreme opposite vals

   void __combine_alleles();

   void __rep(int rep_num, std::string& rseq, bool fwd = true);
   
   /** Construct a breakpoint from a cluster of discordant reads
    */
   BreakPoint(DiscordantCluster& tdc,
	      DiscordantClusterMap& dmap, 
	      const SeqLib::GenomicRegion& region,
	      const SeqLib::BamHeader& h,
	      SvabaSharedConfig& sc);
     
   BreakPoint() {}
   
   BreakPoint(const std::string &line, const SeqLib::BamHeader& h);

   void checkLocal(const SeqLib::GenomicRegion& window);

   void score_somatic(double NODBCUTOFF, double DBCUTOFF);

   void addCovs(const std::unordered_map<std::string, STCoverage*>& covs);

   /** Retrieve the reference sequence at a breakpoint and determine if 
    * it lands on a repeat */
   void repeatFilter();

   std::string print(const SeqLib::BamHeader& h) const;

   void __combine_with_discordant_cluster(DiscordantClusterMap& dmap);
   
   /*! @function determine if the breakpoint has split read support
    * @param reference to a vector of read smart pointers that have been aligned to a contig
    * @discussion Note: will cause an error if the AL tag not filled in for the reads. 
    * The AL tag is filled in by AlignedContig::alignReadsToContigs.
    */
   //void splitCoverage(SeqLib::BamRecordVector &bav);
   void splitCoverage(svabaReadVector &bav);
   
   /*! Determines if the BreakPoint overlays a blacklisted region. If 
    * and overlap is found, sets the blacklist bool to true.
    *
    * Note that currently this only is set for the pos1 of indels.
    * If the BreakPoint object is not an indel, no action is taken. 
    * @param grm An interval tree map created from a BED file containing blacklist regions
    */
   void checkBlacklist(SeqLib::GRC &grv);
   
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
   
   std::string toFileString(bool noreads = false);
   
   bool hasDiscordant() const;
   
   bool operator==(const BreakPoint& bp) const;
   
   // define how to sort these 
   bool operator < (const BreakPoint& bp) const { 

     if (b1.gr < bp.b1.gr)
       return true;
     else if (bp.b1.gr < b1.gr)
       return false;
     
     if (b2.gr < bp.b2.gr)
       return true;
     else if (bp.b2.gr < b2.gr)
       return false;
     
     if (n.split > bp.n.split) 
       return true;
     else if (n.split < bp.n.split)
       return false;
     
     if (t.split > bp.t.split)
       return true;
     else if (t.split < bp.t.split)
       return false;
     
     if (t.split > bp.t.split)
       return true;
     else if (t.split < bp.t.split)
       return false;
     
     if (dc.ncount > bp.dc.ncount)
       return true;
     else if (dc.ncount < bp.dc.ncount)
       return false;
     
     if (dc.tcount > bp.dc.tcount)
       return true;
     else if (dc.tcount < bp.dc.tcount)
       return false;
     
     if (cname > bp.cname)
       return true;
     else if (cname < bp.cname)
       return false;
     
     return false;
  }

   //   friend std::ostream& operator<<(std::ostream& out, const BreakPoint& bp);
   
   void score_dscrd(int min_dscrd_size);
   void score_assembly_only();
   void score_assembly_dscrd();
   void score_indel(double LOD_CUTOFF, double LOD_CUTOFF_DBSNP);
   void format_readname_string();
   void set_homologies_insertions();
   void set_evidence();
   bool valid() const;
   void format_bx_string();

   void setRefAlt(const SeqLib::RefGenome* main_rg);


};

 struct ReducedDiscordantCluster {
   uint32_t mapq1:8, mapq2:8, tcount:8, ncount:8;
 };
 
 struct ReducedBreakPoint {

   // some helper functions
   char* __string_alloc2char(const std::string& str, char * p) {
     if (!str.empty() && str != "x") {
       p = (char*)malloc(str.length() + 1);
       strcpy(p, str.c_str());
       return p;
     } else {
       return nullptr;
     }
   }
   
   inline void smart_check_free(char * p) {
     if (p)
       free(p);
   }

   int getSpan() const {
     if (indel && !insertion) // deletion
       return (abs((int)b1.gr.pos1 - (int)b2.gr.pos1) - 1);
     if (indel) // insertion
       return (strlen(insertion)); // insertion
     if (b1.gr.chr == b2.gr.chr)
       return abs((int)b1.gr.pos1-(int)b2.gr.pos1);
     else
       return -1;

   }

   // define how to sort these  
   bool operator<(const ReducedBreakPoint& bp) const;

   // print it with the correct chromsome string
   std::string print(const BreakPoint& b, const SeqLib::BamHeader& h) const;

   ReducedBreakPoint() {}
   ~ReducedBreakPoint() {
     smart_check_free(ref);
     smart_check_free(alt);
     smart_check_free(cname);
     smart_check_free(homology);
     smart_check_free(insertion);
     smart_check_free(evidence);
     smart_check_free(confidence);
     smart_check_free(repeat);
   }
   ReducedBreakPoint(const std::string &line, const SeqLib::BamHeader& h);

   char * ref;
   char * alt;
   char * cname;
   char * evidence;
   char * confidence;
   char * insertion;
   char * homology;
   char * repeat;

   std::string read_names, bxtable;

   std::vector<std::string> format_s;

   ReducedBreakEnd b1, b2;
   double somatic_score = 0;
   double somatic_lod = 0; // LogOdds that variant not in normal
   double true_lod = 0;

   //uint32_t nsplit:8, tsplit:8, 
   //uint32_t tcov_support:8, ncov_support:8, tcov:8, ncov:8;
   uint32_t cov:16, af_n:7, num_align:5, secondary:1, dbsnp:1, pass:1, blacklist:1, indel:1, imprecise:1;
   uint32_t tcigar:8, ncigar:8, dummy:8, af_t:8; 
   float quality;
   uint8_t pon;

   ReducedDiscordantCluster dc;

 };

 
typedef std::vector<BreakPoint> BPVec;
