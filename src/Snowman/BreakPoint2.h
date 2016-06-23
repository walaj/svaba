#ifndef SNOWTOOLS_BREAKPOINT_H__
#define SNOWTOOLS_BREAKPOINT_H__

#include <cstdlib>
#include <string>
#include <map>
#include <unordered_map>
#include <unordered_set>

#include "htslib/faidx.h"

#include "SnowTools/BWAWrapper.h"
#include "SnowTools/STCoverage.h"

#include "PONFilter.h"
#include "DiscordantCluster.h"

namespace SnowTools {

  // forward declares
  struct BreakPoint;
  
  typedef std::vector<BreakPoint> BPVec;
 
 struct ReducedBreakEnd {
   
   ReducedBreakEnd() {}

   ReducedBreakEnd(const GenomicRegion& g, int mq, const std::string & chr_n);
   
   std::string chr_name;
   //char * chr_name;
   GenomicRegion gr;
   int32_t mapq:8, sub_n:8, nm:16;

 };

 struct BreakEnd {
   
   BreakEnd() { mapq = 0; sub_n = 0; nm = 0; }

   BreakEnd(const GenomicRegion& g, int mq, const std::string & chr_n);
   
   BreakEnd(const BamRead& b);
   
   void checkLocal(const GenomicRegion& window);

   std::string hash(int offset = 0) const;

   std::string id;
   std::string chr_name;
   GenomicRegion gr;

   int mapq = -1;
   int cpos = -1;
   int nm = -1;
   int matchlen = -1;

   std::unordered_map<std::string, int> split;
   std::unordered_map<std::string, double> af;

   int sub_n = -1;
   double as_frac= 0;
   bool local;

   friend std::ostream& operator<<(std::ostream& out, const BreakEnd& b);
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
   
   void __smart_check_free(char * p) {
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

   ReducedBreakPoint() {}
   ~ReducedBreakPoint() {
     __smart_check_free(ref);
     __smart_check_free(alt);
     __smart_check_free(cname);
     __smart_check_free(homology);
     __smart_check_free(insertion);
     __smart_check_free(evidence);
     __smart_check_free(confidence);
     __smart_check_free(repeat);
     //__smart_check_free(read_names);
   }
   ReducedBreakPoint(const std::string &line, bam_hdr_t* h);

   char * ref;
   char * alt;
   char * cname;
   char * evidence;
   char * confidence;
   char * insertion;
   char * homology;
   char * repeat;
   //char * read_names;
   std::string read_names;

   std::vector<std::string> format_s;

   //std::string ref;
   //std::string alt;
   //std::string cname;
   //std::string evidence;
   //std::string confidence;
   //std::string insertion;
   //std::string homology;

   ReducedBreakEnd b1, b2;
   double somatic_score = 0;
   double somatic_lod = 0; // LogOdds that variant not in normal
   double true_lod = 0;

   uint32_t nsplit:8, tsplit:8, af_n:7, num_align:5, secondary:1, dbsnp:1, pass:1, blacklist:1, indel:1, imprecise:1;
   uint32_t tcov_support:8, ncov_support:8, tcov:8, ncov:8;
   uint32_t tcigar:8, ncigar:8, quality:8, af_t:8; 
   uint8_t pon;

   ReducedDiscordantCluster dc;

 };

 struct SampleInfo {

   bool indel;

   int split = 0;
   int cigar = 0;
   int alt =0;
   int clip_cov = 0;
   int cov = 0;
   int disc = 0;
   
   // genotype info
   double GQ = 0;
   double PL = 0;
   std::string genotype;

   double af = 0;
   double error_rate = 1e-4;

   double LO = 0; // log odds of variant vs error
   double SLO = 0; // MAPQ scaled log odds of variant vs error
   double LO_n = 0; // log odds of variant at af=0.5 vs ref (af=0) with errors

   std::set<std::string> supporting_reads; // holds SR tags (not qnames)

   friend std::ostream& operator<<(std::ostream& out, const SampleInfo& a);

   friend SampleInfo operator+(const SampleInfo& a1, const SampleInfo& a2);

   double __log_likelihood(int ref, int alt, double f, double e);

   void modelSelection(double err);

   std::string toFileString() const;

   void fromString(const std::string& s);

   void __adjust_alt_counts();
   
 };
 
 struct BreakPoint {
   
   static std::string header() { 
     return "chr1\tpos1\tstrand1\tchr2\tpos2\tstrand2\tref\talt\tspan\tmapq1\tmapq2\tnm1\tnm2\tdisc_mapq1\tdisc_mapq2\tsub_n1\tsub_n2\thomology\tinsertion\tcontig\tnumalign\tconfidence\tevidence\tquality\tsecondary_alignment\tsomatic_score\tsomatic_lod\ttrue_lod\tpon_samples\trepeat_seq\tgraylist\tDBSNP\treads"; 
   }

   double somatic_score = 0;
   double somatic_lod = 0; // LogOdds that variant not in normal

   std::string seq, cname, rs, insertion, homology, repeat_seq, evidence, confidence, ref, alt, read_names;   

   // the evidence per break-end
   BreakEnd b1, b2;
   
   SampleInfo t, n, a;

   // reads spanning this breakpoint
   BamReadVector reads;

   //int t_reads = 0, n_reads = 0;

   // discordant reads supporting this aseembly bp
   DiscordantCluster dc;
   
   int quality = 0;

   // total coverage at that position
   std::map<std::string, SampleInfo> allele; // ordered to keep in alphabetical order by prefix (e.g. n001)

   bool secondary = false;

   std::unordered_set<std::string> split_reads, qnames;
   
   int pon = 0;
   int num_align = 0;

   bool complex = false;
   
   bool isindel = false;
   bool blacklist = false;

   // keep track of how much of contig is covered by split
   std::pair<int,int> split_cov_bounds = std::pair<int, int>(1e5, -1); // dummy to extreme opposite vals

   void __combine_alleles();

   void __rep(int rep_num, std::string& rseq, bool fwd = true);
   
   /** Construct a breakpoint from a cluster of discordant reads
    */
   BreakPoint(DiscordantCluster& tdc, const BWAWrapper * bwa, SnowTools::DiscordantClusterMap& dmap);
     
   BreakPoint() {}
   
   BreakPoint(const std::string &line, bam_hdr_t* h);

   void checkLocal(const GenomicRegion& window);

   void __set_total_reads();
   
   void __score_somatic(double NODBCUTOFF, double DBCUTOFF);

   void addCovs(const std::unordered_map<std::string, STCoverage*>& covs, const std::unordered_map<std::string, STCoverage*>& clip_covs);

   /** Retrieve the reference sequence at a breakpoint and determine if 
    * it lands on a repeat */
   void repeatFilter();

   void __combine_with_discordant_cluster(DiscordantClusterMap& dmap);
   
   /*! @function determine if the breakpoint has split read support
    * @param reference to a vector of read smart pointers that have been aligned to a contig
    * @discussion Note: will cause an error if the AL tag not filled in for the reads. 
    * The AL tag is filled in by AlignedContig::alignReadsToContigs.
    */
   void splitCoverage(BamReadVector &bav);
   
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
   void scoreBreakpoint(double LOD_CUTOFF, double DBCUTOFF, double NODBCUTOFF, double LRCUTOFF, int min_dscrd_size);
   
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
   void addAllelicFraction(STCoverage * t_cov, STCoverage * n_cov, STCoverage * n_clip_cov);
  
  /*! @function get the span of the breakpoints (in bp). -1 for interchrom
   * @return int distance between breakpoints
   */
   int getSpan() const;

   /*! @function check the breakpoint against a panel of normals
    * @param Panel of normals hash
    * @return number of normal samples with this variant
    */
   int checkPon(const SnowTools::PONFilter * p);
   
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

   friend std::ostream& operator<<(std::ostream& out, const BreakPoint& bp);
   
   void __score_dscrd(int min_dscrd_size);
   void __score_assembly_only();
   void __score_assembly_dscrd();
   void __score_indel(double LOD_CUTOFF);
   std::string __format_readname_string();
   void __set_homologies_insertions();
   void __set_evidence();
   bool valid() const;
   
   double __sv_is_somatic() const;
   double __indel_is_somatic() const;

   void setRefAlt(faidx_t * main_findex, faidx_t * viral_findex);

};

}

#endif
