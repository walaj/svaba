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
#include "SVBamReader.h"
#include "AlignedContig.h"

using namespace std;
using namespace BamTools;

typedef unordered_map<string, EncodedBA> ReadMap;
typedef unordered_map<string, string> RMap;
typedef unordered_map<string, size_t> GMap;

struct IterPair {
  IterPair(BamAlignmentVector::const_iterator s, BamAlignmentVector::const_iterator e) : start(s), end(e) {}
  IterPair() {}
  ~IterPair() {}
  BamAlignmentVector::const_iterator start;
  BamAlignmentVector::const_iterator end;  
};

typedef map<std::string, IterPair> IterPairMap;
typedef vector<std::string> StringVec;

void addDiscordantPairsBreakpoints(BPVec &bp, DMap * dmap);
void clusterReads(BamAlignmentVector &bav, GenomicRegionVector &grv, RMap &rmap, char anchor_strand, char partner_strand);
void finalizeCluster(GenomicRegionVector &grv, RMap &rmap, int pos, GenomicRegion anc, GenomicRegion par, int dtcount, int dncount);
bool grabReads(int refID, int pos1, int pos2, AlignedContigVec * cont_out, DMap * bav_disc);
//void SGAassemble(stringstream &asqg_stream, int minOverlap, int cutoff, string prefix, ContigVector &contigs);
int runAll(SeqRecordVector &srv, AlignedContigVec * cont_out);
bool runTaiga(int argc, char** argv);
void parseTaigaOptions(int argc, char** argv);
void chunkReadsForAssembly(const int refID, const int pos1, const int chunk, const int pad, 
			   AlignedContigVec * cont_out, BamAlignmentVector * tbav, BamAlignmentVector * nbav,
			   DMap * disc_out);
//void addDiscCluster(BamTools::BamAlignment a1, BamTools::BamAlignment a2, size_t cluster);
//bool parseRegionFile(GenomicRegionVector &gr);
void getChunkReads(const BamAlignmentVector * srv, const unsigned refID, const unsigned pos1, const unsigned chunk, const unsigned pad, IterPairMap &mmap);
//void deduplicateReadsPos(const BamAlignmentVector &inbav, BamAlignmentVector &outbav);
//void matchReads2Contigs(ContigVector * contigs, BamAlignmentVector &bav, ContigVector * cont_out);
void doAssembly(ReadTable *pRT, std::string name, ContigVector &contigs, int pass);
void combineR2C(EncodedBAVector &read_in, ReadMap &read_out);
void grabPairmateReads(BamAlignmentVector &bav, const GenomicRegion window, DMap * bav_disc);
int countJobs(GenomicRegionVector &file_regions, GenomicRegionVector &run_regions);
void learnParameters();
void _learn_params(BamTools::BamReader &reader, vector<double> &mapq_result, vector<double> &isize_result,
		   double &inter_cov, double &clip_cov, GenomicRegion &gr, int &readlen);
//void writeDiscBam(BamAlignmentVector * disc);
//void runBWA();
void cleanR2C();
void cleanDiscBam();
//void writeContigFasta(AlignedContigVec *ct);
//void writeReadsBam(AlignedContigVec *ct); // deprecated
//void writeReadsBam(EncodedBAVector *reads);
void handleDiscordant(BamAlignmentVector &bavd, string name, GenomicRegion gr, DMap * bav_disc);
//void clearMemWriteData();
//void ContigsToReadTable(const ContigVector &contigs, ReadTable &pRT);
//void BamAlignmentVectorToReadTable(const BamAlignmentVector &bav, ReadTable &pRT);
void combineContigsWithDiscordantClusters(DMap * dmap, AlignedContigVec * cont_out);

class TaigaWorkItem {

private:
  int m_refid; 
  int m_pos1;
  int m_pos2;
  int m_number;
  AlignedContigVec * m_cont;
  DMap * m_disc;
   
public:
  TaigaWorkItem(int refid, int start, int end, int number, AlignedContigVec * cont_out, DMap * disc_out)  
    : m_refid(refid), m_pos1(start), m_pos2(end), m_number(number), m_cont(cont_out), m_disc(disc_out){}
  ~TaigaWorkItem() {}
 
  int getNumber() { return m_number; }
  int getRefID() { return m_refid; }
  int getPos1() { return m_pos1; }
  int getPos2() { return m_pos2; }

  bool run() { return grabReads(m_refid, m_pos1, m_pos2, m_cont, m_disc); }
  AlignedContigVec* output() { return m_cont; }

};

template<typename T> inline double calc_sd(vector<T> vec) {
  double mean  = accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
  double sqsum = inner_product(vec.begin(), vec.end(), vec.begin(), 0.0);
  double stdev = sqrt(sqsum / vec.size() - mean * mean);
  return stdev;
}

#endif


