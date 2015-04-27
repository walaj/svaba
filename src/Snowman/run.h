#ifndef SNOWMAN_RUN_H
#define SNOWMAN_RUN_H

#include <string>
#include <vector>
#include <pthread.h>
#include <ctime>

#include "ReadTable.h"
#include "workqueue.h"

#include "BamToolsUtils.h"
#include "SnowTools/GenomicRegion.h"
#include "SnowTools/GenomicRegionCollection.h"
#include "SnowTools/SnowUtils.h"

#include "contigs.h"
#include "AlignedContig.h"
#include "DiscordantCluster.h"
#include "SnowmanBamWalker.h"

#include "BWAWrapper.h"

// needed for seq record vector
#include "Util.h" 

using namespace std;
using SnowTools::GRC;
using SnowTools::GR;

typedef unordered_map<string, Read> ReadMap;
//typedef unordered_map<string, unique_ptr<BamAndReads> > BARMap;

void initializeFiles();
void addDiscordantPairsBreakpoints(BPVec &bp, DMap& dmap);
GRC calculateClusters(ReadVec &bav);
DMap clusterDiscordantReads(ReadVec &bav, const GR& interval);
bool _cluster(vector<ReadVec> &cvec, ReadVec &clust, Read &a, bool mate);
bool grabReads(const GR &egion, bwaidx_t* idx);
bool runSnowman(int argc, char** argv);
void parseRunOptions(int argc, char** argv);
void writeR2C(ReadMap &r2c);
void _convertToDiscordantCluster(DMap &dd, vector<ReadVec> cvec, ReadVec &bav);
void doAssembly(ReadTable *pRT, std::string name, ContigVector &contigs, int pass);
int countJobs(GRC& file_regions, GRC& run_regions);
void combineContigsWithDiscordantClusters(DMap &dm, AlignedContigVec &contigs);
void learnParameters();
//SnowTools::GenomicRegionCollection<SnowTools::GenomicRegion> checkReadsMateRegions(SnowTools::GenomicRegionCollection<SnowTools::GenomicRegion> mate_final, unique_ptr<BARMap>& bar);
SeqRecordVector toSeqRecordVector(ReadVec &bav);

//bool findOverlapBlocksExactSnow(const string &w, const BWT* pBWT,
//				const BWT* pRevBWT);


/** @brief p-thread work item that calls Snowman on a small region

    Detailed description follows here.
    @author X. XYZ, DESY
    @date March 2008
*/
class SnowmanWorkItem {

private:
  SnowTools::GenomicRegion m_gr;
  int m_number;  
  bwaidx_t * m_idx;
   
public:
  SnowmanWorkItem(const SnowTools::GenomicRegion& gr, int number, bwaidx_t* idx)  
    : m_gr(gr), m_number(number), m_idx(idx) {}
  ~SnowmanWorkItem() {}
 
  int getNumber() { return m_number; }

  bool run() { return grabReads(m_gr, m_idx); }

};

template<typename T> inline double calc_sd(vector<T> vec) {
  double mean  = accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
  double sqsum = inner_product(vec.begin(), vec.end(), vec.begin(), 0.0);
  double stdev = sqrt(sqsum / vec.size() - mean * mean);
  return stdev;
}


// define a sorter to sort by the Mate Position
typedef std::binary_function<Read, Read, bool> AlignmentSortBase;
struct ByMatePosition : public AlignmentSortBase {
  
  // ctor
  ByMatePosition(const BamTools::Algorithms::Sort::Order& order = BamTools::Algorithms::Sort::AscendingOrder)
    : m_order(order) { }
  
  // comparison function
  bool operator()(const Read& lhs, const Read& rhs) {
    
    // force unmapped aligmnents to end
    if ( r_mid(lhs) == -1 ) return false;
    if ( r_mid(rhs) == -1 ) return true;
    
    // if on same reference, sort on position
    if ( r_mid(lhs) == r_mid(rhs) )
      return BamTools::Algorithms::Sort::sort_helper(m_order, r_mpos(lhs), r_mpos(rhs));
    
    // otherwise sort on reference ID
    return BamTools::Algorithms::Sort::sort_helper(m_order, r_mid(lhs), r_mid(rhs));
  }
  // used by BamMultiReader internals
  static inline bool UsesCharData(void) { return false; }
  
  // data members
  private:
  const BamTools::Algorithms::Sort::Order m_order;
};

// define a sorter to sort by the Mate Position
typedef std::binary_function<Read, Read, bool> AlignmentSortBase;
struct ByReadPosition : public AlignmentSortBase {
  
  // ctor
  ByReadPosition(const BamTools::Algorithms::Sort::Order& order = BamTools::Algorithms::Sort::AscendingOrder)
    : m_order(order) { }
  
  // comparison function
  bool operator()(const Read& lhs, const Read& rhs) {
    
    // force unmapped aligmnents to end
    if ( r_id(lhs) == -1 ) return false;
    if ( r_id(rhs) == -1 ) return true;
    
    // if on same reference, sort on position
    if ( r_id(lhs) == r_id(rhs) )
      return BamTools::Algorithms::Sort::sort_helper(m_order, r_pos(lhs), r_pos(rhs));
    
    // otherwise sort on reference ID
    return BamTools::Algorithms::Sort::sort_helper(m_order, r_id(lhs), r_id(rhs));
  }
  // used by BamMultiReader internals
  static inline bool UsesCharData(void) { return false; }
  
  // data members
  private:
  const BamTools::Algorithms::Sort::Order m_order;
};

// define a sorter to sort by the Mate Position
typedef std::binary_function<Read, Read, bool> AlignmentSortBase;
struct ByReadAndMatePosition : public AlignmentSortBase {
  
  // ctor
  ByReadAndMatePosition(const BamTools::Algorithms::Sort::Order& order = BamTools::Algorithms::Sort::AscendingOrder)
    : m_order(order) { }

  // comparison function
  bool operator()(const Read& lhs, const Read& rhs) {

    // set the lower nums
    int Lref, Lpos, Rref, Rpos;

    // read is lower
    if ( (r_id(lhs) < r_mid(lhs)) || (r_id(lhs) == r_mid(lhs) && r_pos(lhs) < r_mpos(lhs))) {
      Lref = r_id(lhs);
      Lpos = r_pos(lhs);
    // mate is lower
    } else {
      Lref = r_mid(lhs);
      Lpos = r_mpos(lhs);
    }
    
    // read is lower
    if ( (r_id(rhs) < r_mid(rhs)) || (r_id(rhs) == r_mid(rhs) && r_pos(rhs) < r_mpos(rhs))) {
      Rref = r_id(rhs);
      Rpos = r_pos(rhs);
    // mate is lower
    } else {
      Rref = r_mid(rhs);
      Rpos = r_mpos(rhs);
    }

    // force unmapped aligmnents to end
    if ( Lref == -1) return false;
    if ( Rref == -1) return true;
    
    // if on same reference, sort on position
    if ( Lref == Rref )
      return BamTools::Algorithms::Sort::sort_helper(m_order, Lpos, Rpos);
    
    // otherwise sort on reference ID
    return BamTools::Algorithms::Sort::sort_helper(m_order, Lref, Rref);
  }
  // used by BamMultiReader internals
  static inline bool UsesCharData(void) { return false; }
  
  // data members
  private:
  const BamTools::Algorithms::Sort::Order m_order;
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
	     SnowTools::percentCalc<double>(itr->second, total_time),
	     SnowTools::percentCalc<double>(itm->second, total_time),
	     SnowTools::percentCalc<double>(itc->second, total_time),
	     SnowTools::percentCalc<double>(ita->second, total_time),
	     SnowTools::percentCalc<double>(itb->second, total_time),
	     SnowTools::percentCalc<double>(its->second, total_time));
    out << string(buffer);
    return out;
  }

};


#endif


