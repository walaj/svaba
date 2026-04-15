#pragma once

#include <cstdlib>
#include <string>
#include <map>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include "DiscordantCluster.h"
#include "SeqLib/BamRecord.h"
#include "SvabaRead.h"

using CigarMapMap = std::unordered_map<std::string, SeqLib::CigarMap>;

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

enum class SomaticState {
  NOTSET,
  SOMATIC_LOD,    // somatic by LOD cutoff
  NORMAL_LOD,     // normal by LOD cutoff
  FAILED,         // normal - failed hard filter
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
  int mapq = -1; ///< mapping quality of alignment
  int cpos = -1; ///< breakpoint position on the contig (0-based)
  int nm = -1;   ///< number of mismatching bases/indels bases
  int matchlen = -1; /// number of matching bases in the alignment
  int as = -1;  ///< primary alignment score  
  int sub = -1; ///< best secondary alignment score

  // Svaba2.0: composite reliability of the contig-alignment fragment backing
  // this breakend, in [0, 1]. For an SV BreakPoint the two ends can come
  // from different BWA fragments of the same split-aligned contig, so
  // each side stores its own score. Default 1.0 so legacy code paths
  // that never set it don't accidentally demote calls.
  double contig_conf = 1.0;
  
  LocalAlignment local = LocalAlignment::NOTSET;
  
  // for high-confidence reads
  std::unordered_map<std::string, int> split;  
  
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
		       "span\tsplit\talt\tcov\tcigar\tcigar_near\t"
		       "dmq1\tdmq2\tdcn\tdct\t"
		       "mapq1\tmapq2\tnm1\tnm2\tas1\tas2\tsub1\tsub2\t"
		       "homol\tinsert\trepeat\t"
		       "contig_and_region\tnaln\tconf\ttype\tqual\t2ndary\t"
		       "somatic\tsomlod\tmaxlod\tdbsnp\tcontig_conf1\tcontig_conf2"
		       );
  }
  
  SomaticState somatic = SomaticState::NOTSET;

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
  
  int quality = -1;   // QUAL score
  int secondary = -1; // is this a secondary
  int pass = -1;      //false;
  int num_align = 0;  // number of alignments for contigs that generated this

  // SvABA2.0: for indels only, track whether the underlying AlignmentFragment
  // had its contig sequence flipped (i.e. the primary contig alignment is on
  // the reverse strand of the genome). b1.cpos / b2.cpos are computed from
  // the un-flipped BAM CIGAR (so they are BAM/genome-forward coordinates),
  // but m_seq in the AlignedContig and the read-to-contig positions from the
  // r2c BWA index are both in assembly-native orientation. When flipped is
  // true, these two coordinate systems are mirror images and we must convert
  // cpos to m_seq coordinates before comparing with r2c positions or
  // rendering onto m_seq. Use cpos_on_m_seq() below to get the converted
  // (b1, b2) pair.
  int  contig_len        = 0;
  bool flipped_on_contig = false;

  // Return (b1.cpos, b2.cpos) converted to m_seq / r2c (assembly-native)
  // orientation. When flipped_on_contig is false this is the identity.
  // When flipped, the coordinate system reverses and the "before" / "after"
  // sense of an indel swaps, so we mirror AND swap.
  std::pair<int,int> cpos_on_m_seq() const {
    if (!flipped_on_contig || contig_len <= 0)
      return {b1.cpos, b2.cpos};
    const int L = contig_len;
    return {L - 1 - b2.cpos, L - 1 - b1.cpos};
  }
  
  double LO_s = 0; // log odds of variant being somatic (see svabaModels.cpp - SomaticLOD)
  double max_lod = 0; // the highest LOD across all samples
  
  SVType svtype = SVType::NOTSET;
  LocalAlignment local = LocalAlignment::NOTSET;

  const SvabaSharedConfig* sc;
  
  // keep track of how much of contig is covered by split
  // first is left-most position on contig that has a read aligned to it
  // second is right-most
  std::pair<int,int> split_cov_bounds;

  GenomicRegion BreakEndAsGenomicRegionLeft() const;
  
  GenomicRegion BreakEndAsGenomicRegionRight() const;

  std::string printDeletionMarksForAlignmentsFile() const;

  // SvABA2.0: union of UniqueNames of split-supporting reads across all
  // samples. Used by AlignedContig::printToAlignmentsFile to tag each
  // read line with its variant-support kind. SampleInfo is a private
  // nested type, so we expose this method instead of forcing callers
  // to iterate `allele` and reach into SampleInfo directly.
  std::unordered_set<std::string> getAllSupportingReads() const {
    std::unordered_set<std::string> out;
    for (const auto& kv : allele)
      for (const auto& un : kv.second.supporting_reads)
        out.insert(un);
    return out;
  }

  bool isIndel() const;
  
  HashVector getBreakEndHashes();
  
  void __combine_alleles();
  
  void __rep(int rep_num, std::string& rseq, bool fwd = true);
  
  void setLocal(const GenomicRegion& window);
  
  void score_somatic(double error_fwd);

  void indelCigarCheck(const CigarMapMap& cmap);
  
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
  void scoreBreakpoint(); 
   
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
  
  bool isEmpty() const { return (b1.gr.pos1 == 0 && b2.gr.pos1 == 0); }
  
  std::string toFileString(const BamHeader& header) const;
  
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
  
  void score_dscrd(); 
  void score_assembly_only();
  void score_assembly_dscrd();
  void score_indel(); 
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
    int cigar_near = 0; // cigar matches close (but not same) 
    int alt = 0;
    int cov = 0;
    int disc = 0;
    
    // genotype info
    //NB: PL (the Phred-scaled -10*log_10(GL - max(GL))
    //        rounded to int by convention
    double GQ = 0; // GQ take, max 99, quality of the genotype tag
    double NH_GQ = 0; // GQ of 0/0. Higher is more likely to be not hom ref
    std::string genotype; // GT tag
    std::vector<double> genotype_likelihoods = {0,0,0}; // GL tag
    std::vector<int>    phred_likelihoods = {0,0,0}; // PL tag
    
    double af = 0;

    // log odds
    double LO = 0; // log odds of variant vs error
    double SLO = 0; // MAPQ scaled log odds of variant vs error
    double LO_n = 0; // log odds of variant at af=0.5 vs ref (af=0) with errors
    
    // UniqueName set of supporting reads
    // not for read tracking but actually needed for breakpoint scoring
    ReadNameSet supporting_reads; 
    
    void modelSelection(double err, int readlen);
     
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

