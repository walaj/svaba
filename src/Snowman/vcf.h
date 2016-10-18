#ifndef SNOWMAN_VCF_GEN_H
#define SNOWMAN_VCF_GEN_H

#include <memory>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "SeqLib/GenomicRegion.h"

#include "BreakPoint.h"

typedef std::unordered_map<std::string, std::string> InfoMap;
typedef std::unordered_map<std::string, std::string> FormatMap;
typedef std::unordered_map<std::string, std::string> FilterMap;
typedef std::unordered_map<std::string, std::string> SampleMap;
typedef std::unordered_map<std::string, std::string> ContigFieldMap;

typedef std::unordered_map<std::string, bool> SupportingReadsMap;
typedef std::pair<std::string,std::string> FormatPair;
typedef std::unordered_map<std::string, std::pair<std::string,std::string>> FormatRecordMap;

size_t ChrStringToNumber(const std::string& str);
void runVCF(int argc, char** argv);
void parseVCFOptions(int argc, char** argv);
std::string formatReadString(const std::string& readid, char type);

// structure to store a VCF header
struct VCFHeader {

  VCFHeader() {}
  ~VCFHeader() {}

  VCFHeader(std::string file);

  std::string fileformat = "VCFv4.2";
  std::string filedate;
  std::string source;
  std::string reference = "hg19";
  std::string colnames = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

  InfoMap infomap; 
  FilterMap filtermap;
  FormatMap formatmap;
  SampleMap samplemap;

  ContigFieldMap contigfieldmap;
  
  // output it to a string
  friend std::ostream& operator<<(std::ostream& out, const VCFHeader& v);

  //set the filedate string to the current date
  void setCurrentDate();

  //add an info field
  void addInfoField(std::string field, std::string number, std::string type, std::string description);

  void addFilterField(std::string field, std::string description);
  void addFormatField(std::string field, std::string number, std::string type, std::string description);
  void addSampleField(std::string field);

  void addContigField(std::string id, int len);

};

struct VCFEntry {

  VCFEntry() {}
  ~VCFEntry() {}
  VCFEntry(std::string line, std::string method);
  VCFEntry(const BreakEnd& b);

  // data
  std::shared_ptr<ReducedBreakPoint> bp;
  uint32_t id:30, id_num:2;

  std::string getRefString() const;
  std::string getAltString() const;
  std::string getIdString() const;
  std::pair<std::string, std::string> getSampStrings() const;

  // output it to a string
  friend std::ostream& operator<<(std::ostream& out, const VCFEntry& v);

  // define how to sort
  bool operator<(const VCFEntry &v) const;

  bool operator==(const VCFEntry &v) const;

  std::unordered_map<std::string, std::string> fillInfoFields() const;

};

typedef std::vector<VCFEntry> VCFEntryVec;

struct VCFEntryPair {

  VCFEntryPair(std::shared_ptr<ReducedBreakPoint>& b);
  VCFEntryPair() {};
  ~VCFEntryPair() {};

  // data
  VCFEntry e1, e2;
  std::shared_ptr<ReducedBreakPoint> bp;

  //SupportingReadsMap supp_reads;

  bool getOverlaps(int pad, VCFEntryPair &v);

  void addCommonInfoTag(std::string tag, std::string value);

  std::string toCSVString() const;

  // output it to a string
  friend std::ostream& operator<<(std::ostream& out, const VCFEntryPair& v);

};

typedef std::vector<VCFEntryPair> VCFEntryPairVec;
typedef std::unordered_map<int, std::shared_ptr<VCFEntryPair>> VCFEntryPairMap;
typedef std::unordered_map<int, VCFEntry> VCFEntryMap;

// declare a structure to hold the entire VCF
struct VCFFile {

  VCFFile() {}
  ~VCFFile() {}

  VCFFile(std::string file, std::string tmethod);

  // create a VCFFile from a csv
  VCFFile(std::string file, std::string id, const SeqLib::BamHeader& h, const VCFHeader& vheader, bool nopass);

  std::string filename;
  std::string method;

  std::string analysis_id; 

  std::unordered_set<int> dups;

  //  VCFHeader header;
  VCFHeader indel_header;
  VCFHeader sv_header;
  VCFEntryPairMap entry_pairs;
  VCFEntryPairMap indels;
  
  bool include_nonpass = false;

  // output it to a string
  friend std::ostream& operator<<(std::ostream& out, const VCFFile& v);
  
  // write to file
  bool write(std::string basename) const;

  // write to csv file 
  bool writeCSV() const;

  //
  void deduplicate();
  
  //
  void writeIndels(std::string basename, bool zip, bool onefile) const;
  void writeSVs(std::string basename, bool zip, bool onefile) const;
  

};

// 
VCFFile mergeVCFFiles(VCFFile const &v1, VCFFile const &v2);
VCFHeader mergeVCFHeaders(VCFHeader const &h1, VCFHeader const &h2);
template<typename T> T mergeHeaderMaps(T const &m1, T const &m2);
SupportingReadsMap ReadIDToReads(std::string readid);
InfoMap mergeInfoFields(InfoMap const &m1, InfoMap const &m2);
FormatRecordMap FormatStringToFormatRecordMap(std::string format, std::string samp1, std::string samp2);

#endif
