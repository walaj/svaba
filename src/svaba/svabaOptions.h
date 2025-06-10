#pragma once

#include "svaba_params.h"
#include <string>
#include <vector>
#include <map>

class SvabaLogger;

class SvabaOptions {

 public:
  // high-level flags
  bool   help        = false;
  int    verbose     = 0;
  int    numThreads  = 1;
  std::string analysisId;
  bool hp            = false;
  int perRgLearnLimit = 50'000;
  size_t weird_read_limit = 50'000;

  // dumping
  bool dump_weird_reads      = false; //true;
  bool dump_discordant_reads = false; //true;
  bool dump_corrected_reads  = false; //true;
  
    // inputs
  std::vector<std::string> caseBams;
  std::vector<std::string> controlBams;
  std::string refGenome;
  std::string regionFile;

  // make the log verbose, but will invoke a lot of mutex locks
  bool verbose_log = false;

  int windowpad          = 500;

  std::string main_bam;
  
  // mode flags
  bool singleEnd         = false;
  bool allContigs        = false;
  bool discClusterOnly   = false;
  bool overrideRefCheck  = false;

  // numeric thresholds
  double sdDiscCutoff       = 3.92;
  int    chunkSize          = 25000;
  int32_t maxReadsPerAssem  = -1;
  double lod                = 8.0;
  double lodDb              = 6.0;
  double lodSomatic         = 6.0;
  double lodSomaticDb       = 10.0;
  double scaleError         = 1.0;
  int    maxCov             = 100;

  // SGA / assembly
  int    sgaMinOverlap      = 0;
  float  sgaErrorRate       = 0.0f;
  int    sgaNumRounds       = 3;

  // error correction
  std::string ecCorrectType = "f";   //s, f, or 0
  double      ecSubsample   = 0.50;

  // BWA-MEM tuning
  int   bwaGapOpen       = 32;
  int   bwaGapExt        = 1;
  int   bwaMismatch      = 18;
  int   bwaMatchScore    = 2;
  int   bwaZdrop         = 100;
  int   bwaBandwidth     = 1000;
  float bwaReseedTrigger = 1.5f;
  int   bwaClip3         = 5;
  int   bwaClip5         = 5;

  // filtering
  std::string rulesJson;
  size_t      mateLookupMin       = 3;
  size_t      mateRegionLookupLim = 400;
  bool        noBadAvoid          = true;

  // external DBs
  std::string blacklistFile;
  std::string germlineSvFile;
  std::string dbsnpVcf;

  // BAM map (idpath)
  // e.g. t001, /path/to/
  std::map<std::string,std::string> bams;

  // ----------------------------------------------------
  // Print help/usage
  static void printUsage();

  void printLogger(SvabaLogger& logger) const;  

  void addFRRule(const std::string &rgName, int N);
  
  // Parse argc/argv into an SvabaOptions; throws on error
  static SvabaOptions parse(int argc, char** argv);

};
