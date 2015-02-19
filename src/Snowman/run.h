#ifndef SNOWMAN_RUN_H
#define SNOWMAN_RUN_H

#include <string>
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "ReadTable.h"
#include <pthread.h>
#include "workqueue.h"
#include <vector>
#include "contigs.h"
#include "api/algorithms/Sort.h"
#include "GenomicRegion.h"
#include "AlignedContig.h"
#include "DiscordantCluster.h"
#include "BamAndReads.h"
#include <time.h>

using namespace std;
using namespace BamTools;

typedef unordered_map<string, BamAlignmentUP> ReadMap;
typedef unordered_map<string, unique_ptr<BamAndReads> > BARMap;

void initializeFiles();
void addDiscordantPairsBreakpoints(BPVec &bp, DMap& dmap);
GenomicRegionVector calculateClusters(BamAlignmentUPVector &bav);
DMap clusterDiscordantReads(BamAlignmentUPVector &bav);
bool grabReads(int refID, int pos1, int pos2);
bool runSnowman(int argc, char** argv);
void parseRunOptions(int argc, char** argv);
void writeR2C(ReadMap &r2c);
bool _cluster(vector<BamAlignmentUPVector> &cvec, BamAlignmentUPVector &clust, BamAlignmentUP &a, bool mate);
void _convertToDiscordantCluster(DMap &dd, vector<BamAlignmentUPVector> cvec, BamAlignmentUPVector &bav);
void doAssembly(ReadTable *pRT, std::string name, ContigVector &contigs, int pass);
int countJobs(GenomicRegionVector &file_regions, GenomicRegionVector &run_regions);
bool isDiscordant(const BamAlignmentUP &a, bool ancrev, bool parrev);
bool isTumorRead(const BamAlignmentUP &a);
void cleanR2CBam();
void combineContigsWithDiscordantClusters(DMap &dm, AlignedContigVec &contigs);

/** @brief p-thread work item that calls Snowman on a small region

    Detailed description follows here.
    @author X. XYZ, DESY
    @date March 2008
*/
class SnowmanWorkItem {

private:
  int m_refid; 
  int m_pos1;
  int m_pos2;
  int m_number;
   
public:
  SnowmanWorkItem(int refid, int start, int end, int number)  
    : m_refid(refid), m_pos1(start), m_pos2(end), m_number(number) {}
  ~SnowmanWorkItem() {}
 
  int getNumber() { return m_number; }
  int getRefID() { return m_refid; }
  int getPos1() { return m_pos1; }
  int getPos2() { return m_pos2; }

  bool run() { return grabReads(m_refid, m_pos1, m_pos2); }

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

// make a structure to store timing opt
struct SnowTimer {

  SnowTimer() {
    s = {"r", "m", "as", "bw", "cl", "wr", "sw"};
    for (auto& i : s)
      times[i] = 0;
    curr_clock = clock();
  }

  unordered_map<string, double> times;
  vector<string> s;

  clock_t curr_clock;

  void stop(string part) { 
    times[part] += (clock() - curr_clock); 
    curr_clock = clock();
  }
  void start() { curr_clock = clock(); }

  // print it
  friend ostream& operator<<(ostream &out, const SnowTimer st) {

    double total_time = 0;
    for (auto& i : st.times)
      total_time += i.second;
    if (total_time == 0)
      return out;

    char buffer[140];
    
    auto itr = st.times.find("r");
    auto itm = st.times.find("m");
    auto itc = st.times.find("cl");
    auto ita = st.times.find("as");
    auto itb = st.times.find("bw");
    auto its = st.times.find("sw");

    sprintf (buffer, "R: %2d%% M: %2d%% D: %2d%% A: %2d%% B: %2d%% S: %2d%%", 
	     SnowUtils::percentCalc<double>(itr->second, total_time),
	     SnowUtils::percentCalc<double>(itm->second, total_time),
	     SnowUtils::percentCalc<double>(itc->second, total_time),
	     SnowUtils::percentCalc<double>(ita->second, total_time),
	     SnowUtils::percentCalc<double>(itb->second, total_time),
	     SnowUtils::percentCalc<double>(its->second, total_time));
    out << string(buffer);
    return out;
  }

};


#endif


