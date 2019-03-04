#ifndef SVABA_RUN_H__
#define SVABA_RUN_H__

#include <pthread.h>

#include <string>
#include <vector>
#include <ostream>
#include <unordered_map>
#include <unordered_set>
#include <map>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/RefGenome.h"
#include "SeqLib/BFC.h"

#include "svabaUtils.h"
#include "AlignedContig.h"
#include "svabaBamWalker.h"
#include "DiscordantCluster.h"
#include "svabaAssemblerEngine.h"

#include "workqueue.h"

#define SVABA_VERSION "1.1.0"

// typedefs
typedef std::map<std::string, std::string> BamMap;
typedef std::pair<int, int> CountPair;
typedef std::map<std::string, SeqLib::SharedHTSFile> HTSFileMap;

void learnBamParams(SeqLib::BamReader& walk, std::string id);
void makeVCFs();
int overlapSize(const SeqLib::BamRecord& query, const SeqLib::BamRecordVector& subject);
bool hasRepeat(const std::string& seq);
void parseRunOptions(int argc, char** argv);
void runsvaba(int argc, char** argv);
void learnParameters(const SeqLib::GRC& regions);
int countJobs(SeqLib::GRC &file_regions, SeqLib::GRC &run_regions);
void sendThreads(SeqLib::GRC& regions_torun);
bool runWorkItem(const SeqLib::GenomicRegion& region, svabaThreadUnit& wu, long unsigned int thread_id);
SeqLib::GRC makeAssemblyRegions(const SeqLib::GenomicRegion& region);
void alignReadsToContigs(SeqLib::BWAWrapper& bw, const SeqLib::UnalignedSequenceVector& usv, SeqLib::BamRecordVector& bav_this, std::vector<AlignedContig>& this_alc, const SeqLib::RefGenome * rg);
void set_walker_params(svabaBamWalker& walk);
MateRegionVector __collect_normal_mate_regions(WalkerMap& walkers);
MateRegionVector __collect_somatic_mate_regions(WalkerMap& walkers, MateRegionVector& bl);
SeqLib::GRC __get_exclude_on_badness(std::map<std::string, svabaBamWalker>& walkers, const SeqLib::GenomicRegion& region);
void correct_reads(std::vector<char*>& learn_seqs, svabaReadVector& brv);
void run_assembly(const SeqLib::GenomicRegion& region, svabaReadVector& bav_this, std::vector<AlignedContig>& master_alc, 
		  SeqLib::BamRecordVector& master_contigs, SeqLib::BamRecordVector& master_microbial_contigs, DiscordantClusterMap& dmap,
		  std::unordered_map<std::string, SeqLib::CigarMap>& cigmap, SeqLib::RefGenome* refg);
void remove_hardclips(svabaReadVector& brv);
CountPair collect_mate_reads(WalkerMap& walkers, const MateRegionVector& mrv, int round, SeqLib::GRC& this_bad_mate_regions);
CountPair run_mate_collection_loop(const SeqLib::GenomicRegion& region, WalkerMap& wmap, SeqLib::GRC& badd);
void collect_and_clear_reads(WalkerMap& walkers, svabaReadVector& brv, std::vector<char*>& learn_seqs, std::unordered_set<std::string>& dedupe);
void WriteFilesOut(svabaThreadUnit& wu); 
void run_test_assembly();

class svabaWorkItem {

 private:
  SeqLib::GenomicRegion m_gr;
  int m_number;  

 public:
  svabaWorkItem(const SeqLib::GenomicRegion& gr, int number)  
    : m_gr(gr), m_number(number) {}
    ~svabaWorkItem() {}
    
    int getNumber() { return m_number; }
    
    bool run(svabaThreadUnit& wu, long unsigned int thread_id) { 
      return runWorkItem(m_gr, wu, thread_id);
    }
};


#endif
