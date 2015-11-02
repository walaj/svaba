#ifndef SNOWMAN_VCF_GEN_H
#define SNOWMAN_VCF_GEN_H

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "SnowTools/GenomicRegion.h"
#include "SnowTools/BreakPoint2.h"
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

size_t ChrStringToNumber(const std::string& str);
void runVCF(int argc, char** argv);
void parseVCFOptions(int argc, char** argv);
string getRefSequence(const std::string& chr_string, const SnowTools::GenomicRegion& gr, faidx_t *fi);
string formatReadString(const std::string& readid, char type);

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

  void addContigField(string id, int len);

};

struct VCFEntry {

  VCFEntry() {}
  ~VCFEntry() {}
  VCFEntry(string line, string method);
  VCFEntry(const SnowTools::BreakEnd& b);

  // data
  SnowTools::ReducedBreakPoint* bp;
  uint32_t id:30, id_num:2;

  std::string getRefString() const;
  std::string getAltString() const;
  std::string getIdString() const;
  std::pair<std::string, std::string> getSampStrings() const;

  // output it to a string
  friend ostream& operator<<(ostream& out, const VCFEntry& v);

  // define how to sort
  bool operator<(const VCFEntry &v) const;

  bool operator==(const VCFEntry &v) const;

  std::unordered_map<std::string, std::string> fillInfoFields() const;

};

typedef vector<VCFEntry> VCFEntryVec;

struct VCFEntryPair {

  VCFEntryPair(SnowTools::ReducedBreakPoint * b);
  VCFEntryPair() {};
  ~VCFEntryPair() {};

  // data
  VCFEntry e1, e2;
  SnowTools::ReducedBreakPoint * bp;

  //SupportingReadsMap supp_reads;

  bool getOverlaps(int pad, VCFEntryPair &v);

  void addCommonInfoTag(string tag, string value);

  string toCSVString() const;

  // output it to a string
  friend ostream& operator<<(ostream& out, const VCFEntryPair& v);

};

typedef vector<VCFEntryPair> VCFEntryPairVec;
typedef unordered_map<int, VCFEntryPair*> VCFEntryPairMap;
typedef unordered_map<int, VCFEntry> VCFEntryMap;

// declare a structure to hold the entire VCF
struct VCFFile {

  VCFFile() {}
  ~VCFFile() {}

  VCFFile(string file, string tmethod);

  // create a VCFFile from a csv
  VCFFile(string file, string id, bam_hdr_t * h, const VCFHeader& vheader);

  string filename;
  string method;

  string analysis_id; 

  unordered_set<int> dups;

  //  VCFHeader header;
  VCFHeader indel_header;
  VCFHeader sv_header;
  VCFEntryPairMap entry_pairs;
  VCFEntryPairMap indels;
  
  bool include_nonpass = false;

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
