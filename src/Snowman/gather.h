#ifndef SNOW_GATHER_H
#define SNOW_GATHER_H

//#include "Util.h"
#include "AlignedContig.h"
#include <unordered_map>
#include <vector>
#include <string>
#include "api/BamReader.h"
#include "Read2Contig.h"
#include "GenomicRegion.h"

using namespace std; 

typedef vector<AlignedContig> AlignedContigVector;
//typedef unordered_map<string, AlignedContig> ContigMap;
typedef unordered_map<string, BamTools::BamAlignment> BAMap;
typedef vector<string> StringVec;
typedef unordered_map<string, string> StrStrMap;
typedef unordered_map<string, ATuple> AcountMap;

void addDiscordantPairsBreakpoints(BPVec &bp);
void writeDiscordantBreakpoints();
void combineContigsWithDiscordantClusters();
void writeContigFasta(ContigMap * contigs);
void readDiscordantBam();
void parseDiscovarAcount();
bool runConcat(int argc, char** argv);
void parseConcatOptions(int argc, char** argv);
StringVec &csplit(const string &s, char delim, StringVec &elems);
StringVec csplit(const string &s, char delim);
void convertToVCF(BPVec &bp, ContigMap * contigs);
void markDuplicates(BPVec &bp);
void removeDuplicates(BPVec &in, BPVec &out);
//void runBWA();
void readContigBAM(ContigMap * contigs);
void readR2C(ContigMap * contigs);
void addReadsToContigs(ContigMap * contigs, BamAlignmentVector &r2c_vec);
void realignReadsToContigs(ContigMap * contigs);
void sendRepeatMasker(ContigMap * contigs);
void sendBlat(ContigMap * contigs);
StrStrMap readFasta(string file);
//void generateClassMaps(const BPVec &bpvec, ClassMap &som_map, ClassMap &ger_map);
void writeContigs(const ContigMap * contigs);
void writeAsciiPlots(const ContigMap * contigs);
void writeBreakFiles( BPVec &fine,  BPVec &glob);
void writeForRPlots(ContigMap * contigs);
void writeFinalBams(const ContigMap * contigs);
bool runUpdater(ContigMap *contigs);
//bool runContigReader(GenomicRegion reg);
//void readThreadContigs();
void realignThreadReads(ContigMap *contigs);
//void trimContigsUsingReadSupport(ContigMap * contigs_in, ContigMap *contigs_out);
//void trimContigs(ContigMap *contigs, ContigMap *contigs_after);
void writeR2CmatchedBAM(ContigMap * contigs);
void setBreakSomaticGermline(BPVec &bpvec);

/*class ReadContigWorkItem {
  GenomicRegion m_reg;
  public:
  ReadContigWorkItem(GenomicRegion reg) : m_reg(reg) {}
    ~ReadContigWorkItem() {}
    bool run() { return runContigReader(m_reg); }
};*/

class RealignWorkItem {
  ContigMap * m_cont;
  public:
    RealignWorkItem(ContigMap * contig) : m_cont(contig){}
    ~RealignWorkItem() {}
    bool run() { return runUpdater(m_cont); }
};

#endif



