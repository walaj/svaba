#ifndef SNOWMAN_RUN_H__
#define SNOWMAN_RUN_H__

#include <pthread.h>

#include <string>
#include <vector>
#include <ostream>
#include <unordered_map>
#include <map>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/RefGenome.h"

#include "SnowmanUtils.h"
#include "AlignedContig2.h"
#include "SnowmanBamWalker.h"
#include "DiscordantCluster.h"
#include "SnowmanAssemblerEngine.h"

#include "workqueue.h"

void learnBamParams(SeqLib::BamReader& walk, std::string id);
void makeVCFs();
int overlapSize(const SeqLib::BamRecord& query, const SeqLib::BamRecordVector& subject);
bool hasRepeat(const std::string& seq);
void parseRunOptions(int argc, char** argv);
void runSnowman(int argc, char** argv);
void learnParameters(const SeqLib::GRC& regions);
int countJobs(SeqLib::GRC &file_regions, SeqLib::GRC &run_regions);
void sendThreads(SeqLib::GRC& regions_torun);
bool runBigChunk(const SeqLib::GenomicRegion& region); 
SeqLib::GRC makeAssemblyRegions(const SeqLib::GenomicRegion& region);
void alignReadsToContigs(SeqLib::BWAWrapper& bw, const SeqLib::UnalignedSequenceVector& usv, SeqLib::BamRecordVector& bav_this, std::vector<AlignedContig>& this_alc, SeqLib::RefGenome * rg);
SnowmanBamWalker __make_walkers(const std::string& p, const std::string& b, const SeqLib::GenomicRegion& region, int& tcount, int& ncount);
MateRegionVector __collect_normal_mate_regions(std::map<std::string, SnowmanBamWalker>& walkers);
MateRegionVector __collect_somatic_mate_regions(std::map<std::string, SnowmanBamWalker>& walkers, MateRegionVector& bl);
SeqLib::GRC __get_exclude_on_badness(std::map<std::string, SnowmanBamWalker>& walkers, const SeqLib::GenomicRegion& region);
void run_test_assembly();
/** @brief p-thread work item that calls Snowman on a small region

    Detailed description follows here.
    @author X. XYZ, DESY
    @date March 2008
*/

class SnowmanWorkItem {

 private:
  SeqLib::GenomicRegion m_gr;
  int m_number;  

 public:
  SnowmanWorkItem(const SeqLib::GenomicRegion& gr, int number)  
    : m_gr(gr), m_number(number) {}
    ~SnowmanWorkItem() {}
    
    int getNumber() { return m_number; }
    
    bool run() { return runBigChunk(m_gr); }
    
};


#endif
