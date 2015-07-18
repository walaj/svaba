#ifndef SNOWMAN_VCF_GEN_H
#define SNOWMAN_VCF_GEN_H

#include <string>
#include <unordered_map>
#include <vector>
#include "SnowTools/GenomicRegion.h"
#include "SnowTools/AlignedContig.h"
//#include "faidx.h"

using namespace std;

typedef unordered_map<string, string> InfoMap;
typedef unordered_map<string, string> FormatMap;
typedef unordered_map<string, string> FilterMap;
typedef unordered_map<string, string> SampleMap;
typedef unordered_map<string, string> ContigFieldMap;

typedef unordered_map<string, bool> SupportingReadsMap;
typedef pair<string,string> FormatPair;
typedef unordered_map<string, pair<string,string>> FormatRecordMap;

void runVCF(int argc, char** argv);
void parseVCFOptions(int argc, char** argv);
string getRefSequence(const SnowTools::GenomicRegion &gr, faidx_t *fi);

// structure to store a VCF header
struct VCFHeader {

  VCFHeader() {}
  ~VCFHeader() {}

  VCFHeader(string file);

  string fileformat = "VCFv4.2";
  string filedate;
  string source;
  string reference = "hg19";
  string colnames = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

  InfoMap infomap; 
  FilterMap filtermap;
  FormatMap formatmap;
  SampleMap samplemap;

  ContigFieldMap contigfieldmap;
  
  // output it to a string
  friend ostream& operator<<(ostream& out, const VCFHeader& v);

  //set the filedate string to the current date
  void setCurrentDate();

  //add an info field
  void addInfoField(string field, string number, string type, string description);

  void addFilterField(string field, string description);
  void addFormatField(string field, string number, string type, string description);
  void addSampleField(string field);

  void addContigField(string id, string assembly, string length, string species);

};

struct VCFEntry {

  VCFEntry() {}
  ~VCFEntry() {}
  VCFEntry(string line, string method);

  int chr = 0;
  int pos = 0;
  string id;
  string ref;
  string alt;
  string qual = ".";
  string filter = "NA";
  string format;
  string samp1;
  string samp2;  
  
  string idcommon;
  //string method;

  unordered_map<string, string> info_fields;
  FormatRecordMap format_fields;

  // output it to a string
  friend ostream& operator<<(ostream& out, const VCFEntry& v);

  // define how to sort
  bool operator<(const VCFEntry &v) const;

  bool operator==(const VCFEntry &v) const;

};

typedef vector<VCFEntry> VCFEntryVec;

struct VCFEntryPair {

  VCFEntryPair(VCFEntry t1, VCFEntry t2) : e1(t1), e2(t2) {}
  VCFEntryPair(VCFEntry t1) { e1 = t1; e2 = t1; e2.idcommon = "NULL"; }
  VCFEntryPair(VCFEntryPair &v1, VCFEntryPair &v2);
  VCFEntryPair() {};
  ~VCFEntryPair() {};

  VCFEntry e1;
  VCFEntry e2;

  string method;
  string idcommon;

  vector<string> samples;

  string overlap_partner = "";

  SupportingReadsMap supp_reads;

  bool getOverlaps(int pad, VCFEntryPair &v);

  bool hasOverlap() const { return overlap_partner != ""; }
  
  void addCommonInfoTag(string tag, string value);

  string toCSVString() const;

  int tsplit = 0;
  int nsplit = 0;
  int tdisc = 0;
  int ndisc = 0;

  bool indel = false;

  // output it to a string
  friend ostream& operator<<(ostream& out, const VCFEntryPair& v);

};

typedef vector<VCFEntryPair> VCFEntryPairVec;
typedef unordered_map<string, VCFEntryPair> VCFEntryPairMap;
typedef unordered_map<string, VCFEntry> VCFEntryMap;

// declare a structure to hold the entire VCF
struct VCFFile {

  VCFFile() {}
  ~VCFFile() {}

  VCFFile(string file, string tmethod);

  // create a VCFFile from a csv
  VCFFile(string file, const char* index, char sep, string analysis_id);

  string filename;
  string method;

  unordered_map<string, bool> dups;

  //  VCFHeader header;
  VCFHeader indel_header;
  VCFHeader sv_header;
  VCFEntryPairMap entry_pairs;
  VCFEntryMap indels;
  

  // output it to a string
  friend ostream& operator<<(ostream& out, const VCFFile& v);
  
  // write to file
  bool write(string basename) const;

  // write to csv file 
  bool writeCSV() const;

  //
  void deduplicate();
  
  //
  void writeIndels(string basename, bool zip) const;
  void writeSVs(string basename, bool zip) const;
  

};

// 
VCFFile mergeVCFFiles(VCFFile const &v1, VCFFile const &v2);
VCFHeader mergeVCFHeaders(VCFHeader const &h1, VCFHeader const &h2);
template<typename T> T mergeHeaderMaps(T const &m1, T const &m2);
SupportingReadsMap ReadIDToReads(string readid);
InfoMap mergeInfoFields(InfoMap const &m1, InfoMap const &m2);
FormatRecordMap FormatStringToFormatRecordMap(string format, string samp1, string samp2);

#endif
