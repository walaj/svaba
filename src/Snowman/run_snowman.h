#ifndef SNOWMAN_RUN_H__
#define SNOWMAN_RUN_H__

#include <string>
#include <vector>
#include <unordered_map>
#include <pthread.h>
#include <ostream>

#include "SnowTools/GenomicRegionCollection.h"
#include "SnowTools/BWAWrapper.h"
#include "SnowTools/AlignedContig2.h"
#include "SnowmanAssemblerEngine.h"
#include "SnowTools/DiscordantCluster.h"
#include "SnowTools/BamWalker.h"

#include "SnowmanUtils.h"

#include "SnowmanBamWalker.h"
#include "workqueue.h"

void learnBamParams(SnowTools::BamWalker& walk, std::string id);
void makeVCFs();
int overlapSize(const BamRead& query, const BamReadVector& subject);
bool hasRepeat(const std::string& seq);
void parseRunOptions(int argc, char** argv);
void runSnowman(int argc, char** argv);
void learnParameters(const SnowTools::GRC& regions);
int countJobs(SnowTools::GRC &file_regions, SnowTools::GRC &run_regions);
void sendThreads(SnowTools::GRC& regions_torun);
bool runBigChunk(const SnowTools::GenomicRegion& region); 
SnowTools::GRC makeAssemblyRegions(const SnowTools::GenomicRegion& region);
void alignReadsToContigs(SnowTools::BWAWrapper& bw, const SnowTools::USeqVector& usv, BamReadVector& bav_this, std::vector<SnowTools::AlignedContig>& this_alc);
SnowmanBamWalker __make_walkers(const std::string& p, const std::string& b, const SnowTools::GenomicRegion& region, int& tcount, int& ncount);
MateRegionVector __collect_normal_mate_regions(std::unordered_map<std::string, SnowmanBamWalker>& walkers);
MateRegionVector __collect_somatic_mate_regions(std::unordered_map<std::string, SnowmanBamWalker>& walkers, MateRegionVector& bl);
SnowTools::GRC __get_exclude_on_badness(std::unordered_map<std::string, SnowmanBamWalker>& walkers, const SnowTools::GenomicRegion& region);
/** @brief p-thread work item that calls Snowman on a small region

    Detailed description follows here.
    @author X. XYZ, DESY
    @date March 2008
*/

class SnowmanWorkItem {

 private:
  SnowTools::GenomicRegion m_gr;
  int m_number;  

 public:
  SnowmanWorkItem(const SnowTools::GenomicRegion& gr, int number)  
    : m_gr(gr), m_number(number) {}
    ~SnowmanWorkItem() {}
    
    int getNumber() { return m_number; }
    
    bool run() { return runBigChunk(m_gr); }
    
};


#endif
