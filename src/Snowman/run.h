#ifndef GRAB_READS_H
#define GRAB_READS_H

#include <string>
#include <iostream>
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "ReadTable.h"
#include <pthread.h>
#include "workqueue.h"
#include <vector>
#include "contigs.h"
#include "api/algorithms/Sort.h"
#include "GenomicRegion.h"
//#include "SVBamReader.h"
#include "AlignedContig.h"
#include "VariantBamReader.h"
#include "DiscordantCluster.h"
#include "BamAndReads.h"
#include "EncodedString.h"

using namespace std;
using namespace BamTools;

struct EncodedBA {
  
  BamAlignment a;
  DNAEncodedString s;

};
typedef vector<EncodedBA> EncodedBAVector;


//typedef unordered_map<string, EncodedBA> ReadMap;
typedef unordered_map<string, BamAlignmentUP> ReadMap;
typedef unordered_map<string, string> RMap;
typedef pair<BamAlignmentUP, string> ReadAndCluster;
typedef unordered_map<BamAlignmentUP, string> ReadAndClusterMap;
typedef unordered_map<string, DiscordantCluster> ClusterMap; // tumor / normal

typedef vector<BamAlignment*> BamAlignmentVectorP;

typedef unordered_map<string, BamAndReads> BARMap;

struct IterPair {
  IterPair(BamAlignmentVector::const_iterator s, BamAlignmentVector::const_iterator e) : start(s), end(e) {}
  IterPair() {}
  ~IterPair() {}
  BamAlignmentVector::const_iterator start;
  BamAlignmentVector::const_iterator end;  
  size_t size() const { return end - start; }
};

typedef map<std::string, IterPair> IterPairMap;
typedef vector<std::string> StringVec;

//int runAll(vector<IterPair>& ivec, string name, AlignedContigVec * cont_out, DMap * dmap, BamAlignmentVectorAncParPairVector &vec);
void addDiscordantPairsBreakpoints(BPVec &bp, DMap * dmap);
//void clusterReads(BamAlignmentVector &bav, GenomicRegionVector &grv, RMap &rmap, char anchor_strand, char partner_strand);
void clusterReads(BamAlignmentUPVector &bav, ReadMap &discordant_reads, GenomicRegionVector &grv, char anchor_strand, char partner_strand);
void finalizeCluster(GenomicRegionVector &grv, BamAlignmentUPVector &reads_buffer, 
       ReadMap &disc_reads, int pos, GenomicRegion anc, GenomicRegion par, 
		     int dtcount, int dncount);
GenomicRegionVector calculateClusters(BamAlignmentUPVector &bav);
DMap clusterDiscordantReads(BamAlignmentUPVector &bav);
bool grabReads(int refID, int pos1, int pos2, AlignedContigVec * cont_out);
//void SGAassemble(stringstream &asqg_stream, int minOverlap, int cutoff, string prefix, ContigVector &contigs);
bool runSnowman(int argc, char** argv);
void parseRunOptions(int argc, char** argv);
void writeR2C(bool makeIndex = false);
bool _cluster(vector<BamAlignmentUPVector> &main, BamAlignmentUPVector &buff, pair<int,int> &last_info, pair<int, int> this_info, BamAlignmentUP &a);


//bool read_lt_mate(BamAlignmentUP &a) {
//  return ( (a->RefID < a->MateRefID) || ((a->RefID == a->MateRefID) && (a->Position < a->MatePosition)) );
//}

//void chunkReadsForAssembly(const int refID, const int pos1, const int chunk, const int pad, 
//			   AlignedContigVec * cont_out, BamAlignmentVector * tbav, BamAlignmentVector * nbav,
//			   DMap * dmap);
//void addDiscCluster(BamTools::BamAlignment a1, BamTools::BamAlignment a2, size_t cluster);
//bool parseRegionFile(GenomicRegionVector &gr);
//void getChunkReads(const BamAlignmentVector * srv, const unsigned refID, const unsigned pos1, const unsigned chunk, const unsigned pad, IterPairMap &mmap);
//void deduplicateReadsPos(const BamAlignmentVector &inbav, BamAlignmentVector &outbav);
//void matchReads2Contigs(ContigVector * contigs, BamAlignmentVector &bav, ContigVector * cont_out);
void doAssembly(ReadTable *pRT, std::string name, ContigVector &contigs, int pass);
//void combineR2C(EncodedBAVector &read_in, ReadMap &read_out);
//void grabPairmateReads(vector<IterPair>& ivec, const GenomicRegion window, DMap * dmap, BamAlignmentVectorAncParPairVector &vec);
int countJobs(GenomicRegionVector &file_regions, GenomicRegionVector &run_regions);
//void learnParameters();
//void _learn_params(BamTools::BamReader &reader, vector<double> &mapq_result, vector<double> &isize_result,
//		   double &inter_cov, double &clip_cov, GenomicRegion &gr, int &readlen);
//void writeDiscBam(BamAlignmentVector * disc);
//void runBWA();
bool isDiscordant(const BamAlignmentUP &a, bool ancrev, bool parrev);
bool isTumorRead(const BamAlignmentUP &a);
void cleanR2C();
void cleanR2CBam();
//void cleanDiscBam();
//void writeContigFasta(AlignedContigVec *ct);
//void writeReadsBam(AlignedContigVec *ct); // deprecated
//void writeReadsBam(EncodedBAVector *reads);
//void handleDiscordant(BamAlignmentVector &bavd, string name, GenomicRegion gr, DMap * dmap);
//void clearMemWriteData();
//void ContigsToReadTable(const ContigVector &contigs, ReadTable &pRT);
//void BamAlignmentVectorToReadTable(const BamAlignmentVector &bav, ReadTable &pRT);
void combineContigsWithDiscordantClusters(DMap this_dmap, AlignedContigVec * cont_out);

class SnowmanWorkItem {

private:
  int m_refid; 
  int m_pos1;
  int m_pos2;
  int m_number;
  AlignedContigVec * m_cont;
  DMap * m_disc;
   
public:
  SnowmanWorkItem(int refid, int start, int end, int number, AlignedContigVec * cont_out)  
    : m_refid(refid), m_pos1(start), m_pos2(end), m_number(number), m_cont(cont_out){}
  ~SnowmanWorkItem() {}
 
  int getNumber() { return m_number; }
  int getRefID() { return m_refid; }
  int getPos1() { return m_pos1; }
  int getPos2() { return m_pos2; }

  bool run() { return grabReads(m_refid, m_pos1, m_pos2, m_cont); }
  AlignedContigVec* output() { return m_cont; }

};

template<typename T> inline double calc_sd(vector<T> vec) {
  double mean  = accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
  double sqsum = inner_product(vec.begin(), vec.end(), vec.begin(), 0.0);
  double stdev = sqrt(sqsum / vec.size() - mean * mean);
  return stdev;
}


// define a sorter to sort by the Mate Position
typedef std::binary_function<BamAlignmentUP, BamAlignmentUP, bool> AlignmentSortBase;
struct ByMatePosition : public AlignmentSortBase {
  
  // ctor
  ByMatePosition(const Algorithms::Sort::Order& order = Algorithms::Sort::AscendingOrder)
    : m_order(order) { }
  
  // comparison function
  bool operator()(const BamAlignmentUP& lhs, const BamAlignmentUP& rhs) {
    
    // force unmapped aligmnents to end
    if ( lhs->MateRefID == -1 ) return false;
    if ( rhs->MateRefID == -1 ) return true;
    
    // if on same reference, sort on position
    if ( lhs->MateRefID == rhs->MateRefID )
      return Algorithms::Sort::sort_helper(m_order, lhs->MatePosition, rhs->MatePosition);
    
    // otherwise sort on reference ID
    return Algorithms::Sort::sort_helper(m_order, lhs->MateRefID, rhs->MateRefID);
  }
  // used by BamMultiReader internals
  static inline bool UsesCharData(void) { return false; }
  
  // data members
  private:
  const Algorithms::Sort::Order m_order;
};

// define a sorter to sort by the Mate Position
typedef std::binary_function<BamAlignmentUP, BamAlignmentUP, bool> AlignmentSortBase;
struct ByReadPosition : public AlignmentSortBase {
  
  // ctor
  ByReadPosition(const Algorithms::Sort::Order& order = Algorithms::Sort::AscendingOrder)
    : m_order(order) { }
  
  // comparison function
  bool operator()(const BamAlignmentUP& lhs, const BamAlignmentUP& rhs) {
    
    // force unmapped aligmnents to end
    if ( lhs->RefID == -1 ) return false;
    if ( rhs->RefID == -1 ) return true;
    
    // if on same reference, sort on position
    if ( lhs->RefID == rhs->RefID )
      return Algorithms::Sort::sort_helper(m_order, lhs->Position, rhs->Position);
    
    // otherwise sort on reference ID
    return Algorithms::Sort::sort_helper(m_order, lhs->RefID, rhs->RefID);
  }
  // used by BamMultiReader internals
  static inline bool UsesCharData(void) { return false; }
  
  // data members
  private:
  const Algorithms::Sort::Order m_order;
};

// define a sorter to sort by the Mate Position
typedef std::binary_function<BamAlignmentUP, BamAlignmentUP, bool> AlignmentSortBase;
struct ByReadAndMatePosition : public AlignmentSortBase {
  
  // ctor
  ByReadAndMatePosition(const Algorithms::Sort::Order& order = Algorithms::Sort::AscendingOrder)
    : m_order(order) { }

  // comparison function
  bool operator()(const BamAlignmentUP& lhs, const BamAlignmentUP& rhs) {

    // set the lower nums
    int Lref, Lpos, Rref, Rpos;

    // read is lower
    if ( (lhs->RefID < lhs->MateRefID) || (lhs->RefID == lhs->MateRefID && lhs->Position < lhs->MatePosition)) {
      Lref = lhs->RefID;
      Lpos = lhs->Position;
    // mate is lower
    } else {
      Lref = lhs->MateRefID;
      Lpos = lhs->MatePosition;
    }
    
    // read is lower
    if ( (rhs->RefID < rhs->MateRefID) || (rhs->RefID == rhs->MateRefID && rhs->Position < rhs->MatePosition)) {
      Rref = rhs->RefID;
      Rpos = rhs->Position;
    // mate is lower
    } else {
      Rref = rhs->MateRefID;
      Rpos = rhs->MatePosition;
    }

    // force unmapped aligmnents to end
    if ( Lref == -1) return false;
    if ( Rref == -1) return true;
    
    // if on same reference, sort on position
    if ( Lref == Rref )
      return Algorithms::Sort::sort_helper(m_order, Lpos, Rpos);
    
    // otherwise sort on reference ID
    return Algorithms::Sort::sort_helper(m_order, Lref, Rref);
  }
  // used by BamMultiReader internals
  static inline bool UsesCharData(void) { return false; }
  
  // data members
  private:
  const Algorithms::Sort::Order m_order;
};


#endif


