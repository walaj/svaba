#include "svabaOptions.h"

#include <sstream>           // for std::ostringstream
#include <iomanip>           // for std::setw and std::setfill
#include <getopt.h>
#include <iostream>
#include <stdexcept>

#include "svabaLogger.h"

void SvabaOptions::printUsage() {
  std::cout << R"(
Usage: svaba run [OPTIONS]

General:
  -h, --help              Show this message
  -v, --verbose <N>       Verbosity level (0-4), default 0
  -p, --threads <N>       Number of threads, default 1
  -a, --analysis-id <ID>  Identifier for this run

Input:
  -t, --case-bam <FILE>     Tumor BAM; may repeat
  -n, --control-bam <FILE>  Normal BAM; may repeat
  -G, --reference-genome <FILE>
                           Indexed reference FASTA
  -k, --region-file <FILE>  BED or samtools-style regions

Mode:
      --single-end        Single-end mode (no mate lookup)
      --all-contigs       Output all assembled contigs
      --discordant-only   Skip assembly, only discordants
      --override-ref-check
                         Skip BAM vs REF compatibility check

Assembly:
      --min-overlap <bp>  Minimum read overlap for SGA
      --error-rate <f>    SGA fractional error rate
      --rounds <N>        Number of assembly rounds

Error-correction:
      --ec-type <s|f|0>   s=SGA k-mer; f=Fermi BFC; 0=off
      --ec-subsample <f>  Fraction to sample for EC learning

Discordant clustering:
      --disc-sd <f>       SD cutoff, default 3.92

Filtering:
      --max-cov <N>       Max coverage to assemble, default 100
      --mate-min <N>      Min reads to trigger mate lookup, default 3
      --mate-lim <N>      Max reads in mate lookup, default 400

BWA-MEM tuning:
      --bw-op <N>         Gap-open pen, default 32
      --bw-ep <N>         Gap-ext pen, default 1
      --bw-mm <N>         Mismatch pen, default 18
      --bw-ms <N>         Match score, default 2
      --bw-zd <N>         Z-drop, default 100
      --bw-bw <N>         Bandwidth, default 1000
      --bw-rt <f>         Reseed trigger, default 1.5
      --bw-c3 <N>         3' clip pen, default 5
      --bw-c5 <N>         5' clip pen, default 5

Output & DBs:
      --blacklist <FILE>  BED of blacklisted regions
      --germline-sv <FILE>
                         BED of known germline SVs
      --dbsnp <VCF>       DBSNP VCF of known variants
)" << "\n";
}

SvabaOptions SvabaOptions::parse(int argc, char** argv) {
  static const char* shortOpts = "hv:p:a:t:n:G:k:";
  static struct option longOpts[] = {
    {"help",                 no_argument,       nullptr, 'h'},
    {"verbose",     required_argument, nullptr, 'v'},
    {"threads",     required_argument, nullptr, 'p'},
    {"analysis-id", required_argument, nullptr, 'a'},
    {"case-bam",    required_argument, nullptr, 't'},
    {"control-bam", required_argument, nullptr, 'n'},
    {"reference-genome", required_argument, nullptr, 'G'},
    {"region-file", required_argument, nullptr, 'k'},
    {"single-end",       no_argument,       nullptr,  1001},
    {"all-contigs",      no_argument,       nullptr,  1002},
    {"discordant-only",  no_argument,       nullptr,  1003},
    {"override-ref-check",no_argument,      nullptr,  1004},
    {"min-overlap", required_argument,     nullptr,  1100},
    {"error-rate", required_argument,      nullptr,  1101},
    {"rounds",     required_argument,      nullptr,  1102},
    {"ec-type",    required_argument,      nullptr,  1200},
    {"ec-subsample",required_argument,     nullptr,  1201},
    {"disc-sd",    required_argument,      nullptr,  1300},
    {"max-cov",    required_argument,      nullptr,  1400},
    {"mate-min",   required_argument,      nullptr,  1401},
    {"mate-lim",   required_argument,      nullptr,  1402},
    {"chunk-size", required_argument,      nullptr,  1500},
    {"bw-op",      required_argument,      nullptr,  1600},
    {"bw-ep",      required_argument,      nullptr,  1601},
    {"bw-mm",      required_argument,      nullptr,  1602},
    {"bw-ms",      required_argument,      nullptr,  1603},
    {"bw-zd",      required_argument,      nullptr,  1604},
    {"bw-bw",      required_argument,      nullptr,  1605},
    {"bw-rt",      required_argument,      nullptr,  1606},
    {"bw-c3",      required_argument,      nullptr,  1607},
    {"bw-c5",      required_argument,      nullptr,  1608},
    {"blacklist",  required_argument,      nullptr,  1700},
    {"germline-sv",required_argument,      nullptr,  1701},
    {"dbsnp",      required_argument,      nullptr,  1702},
    {nullptr,0,nullptr,0}
  };

  SvabaOptions o;
  int idx, c;
  while ((c = getopt_long(argc, argv, shortOpts, longOpts, &idx)) != -1) {
    switch (c) {
      case 'h':   o.help = true;  return o;
      case 'v':   o.verbose      = std::stoi(optarg); break;
      case 'p':   o.numThreads   = std::stoi(optarg); break;
      case 'a':   o.analysisId   = optarg;            break;
      case 't':   o.caseBams   .push_back(optarg);    break;
      case 'n':   o.controlBams.push_back(optarg);    break;
      case 'G':   o.refGenome    = optarg;            break;
      case 'k':   o.regionFile   = optarg;            break;

      case 1001: o.singleEnd        = true; break;
      case 1002: o.allContigs       = true; break;
      case 1003: o.discClusterOnly  = true; break;
      case 1004: o.overrideRefCheck = true; break;

      case 1100: o.sgaMinOverlap    = std::stoi(optarg); break;
      case 1101: o.sgaErrorRate     = std::stof(optarg); break;
      case 1102: o.sgaNumRounds     = std::stoi(optarg); break;

      case 1200: o.ecCorrectType    = optarg;           break;
      case 1201: o.ecSubsample      = std::stod(optarg); break;

      case 1300: o.sdDiscCutoff     = std::stod(optarg); break;

      case 1400: o.maxCov           = std::stoi(optarg); break;
      case 1401: o.mateLookupMin    = std::stoul(optarg);break;
      case 1402: o.mateRegionLookupLim = std::stoul(optarg); break;

      case 1500: o.chunkSize        = std::stoi(optarg); break;

      case 1600: o.bwaGapOpen       = std::stoi(optarg); break;
      case 1601: o.bwaGapExt        = std::stoi(optarg); break;
      case 1602: o.bwaMismatch      = std::stoi(optarg); break;
      case 1603: o.bwaMatchScore    = std::stoi(optarg); break;
      case 1604: o.bwaZdrop         = std::stoi(optarg); break;
      case 1605: o.bwaBandwidth     = std::stoi(optarg); break;
      case 1606: o.bwaReseedTrigger = std::stof(optarg); break;
      case 1607: o.bwaClip3         = std::stoi(optarg); break;
      case 1608: o.bwaClip5         = std::stoi(optarg); break;

      case 1700: o.blacklistFile  = optarg; break;
      case 1701: o.germlineSvFile = optarg; break;
      case 1702: o.dbsnpVcf       = optarg; break;

      case '?':
      default:
        throw std::runtime_error("Unknown or malformed option; see --help");
    }
  }

  // post validation
  if (!o.help) {
    if (o.caseBams.empty())   throw std::runtime_error("Need at least one --case-bam");
    if (o.refGenome.empty())  throw std::runtime_error("Must supply --reference-genome");
    if (o.analysisId.empty()) o.analysisId = "no_id";
    if (o.chunkSize <= 0)      o.chunkSize = 25000;
  }

  // make the bam map from the file names
  int tumorIdx = 1;
  for (const auto& path : o.caseBams) {
    std::ostringstream ss;
    ss << 't'
       << std::setw(3) << std::setfill('0')
       << tumorIdx++;
    o.bams[ ss.str() ] = path;
  }
  
  int normalIdx = 1;
  for (const auto& path : o.controlBams) {
    std::ostringstream ss;
    ss << 'n'
       << std::setw(3) << std::setfill('0')
       << normalIdx++;
    o.bams[ ss.str() ] = path;
  }
  
  // set the rules to skip read learning if doing stdin
  if (o.bams["t001"] == "-" && !o.rulesJson.empty())
    o.rulesJson = R"json(
{
  "global": {
    "duplicate": false,
    "qcfail":   false
  },
  "": {
    "rules": [
      { "isize": 2000 },
      { "rr":    true },
      { "ff":    true },
      { "rf":    true },
      { "ic":    true },
      { "clip":  5, "length": 30 },
      { "ins":   true },
      { "del":   true },
      { "mapped":      true, "mate_mapped": false },
      { "mate_mapped": true, "mapped":      false },
      { "nm": [3, 0] }
    ]
  }
}
)json";
  
  // if read stream, treat as single-end
  if (o.bams["t001"] == "-")
    o.singleEnd = true;

  // set the "main bam"
  o.main_bam = o.bams["t001"];
  
  return o;
}

void SvabaOptions::printLogger(SvabaLogger& logger) const {
  logger.log(verbose > 1, true, "***************************** PARAMS ****************************");
  logger.log(verbose > 1, true, "    DBSNP Database file: ",            dbsnpVcf);
  logger.log(verbose > 1, true, "    Max cov to assemble: ",            maxCov);
  logger.log(verbose > 1, true, "    Error correction mode: ",          ecCorrectType);
  logger.log(verbose > 1, true, "    Subsample-rate for correction learning: ", ecSubsample);
  logger.log(verbose > 1, true, "    ErrorRate: ", 
             (sgaErrorRate < 0.001f ? std::string("EXACT (0)") 
                                        : std::to_string(sgaErrorRate)));
  logger.log(verbose > 1, true, "    Num assembly rounds: ",           sgaNumRounds);
  logger.log(verbose > 1, true, "    Discordant read extract SD cutoff: ",  sdDiscCutoff);
  logger.log(verbose > 1, true, "    Discordant cluster std-dev cutoff: ",  sdDiscCutoff);
  logger.log(verbose > 1, true, "    Minimum number of reads for mate lookup: ", mateLookupMin);
  logger.log(verbose > 1, true, "    LOD cutoff (non-REF): ",           lod);
  logger.log(verbose > 1, true, "    LOD cutoff (non-REF, at DBSNP): ",  lodDb);
  logger.log(verbose > 1, true, "    LOD somatic cutoff: ",            lodSomatic);
  logger.log(verbose > 1, true, "    LOD somatic cutoff (at DBSNP): ", lodSomaticDb);
  logger.log(verbose > 1, true, "    BWA-MEM params:");
  logger.log(verbose > 1, true, "      Gap open penalty: ",           bwaGapOpen);
  logger.log(verbose > 1, true, "      Gap extension penalty: ",      bwaGapExt);
	     logger.log(verbose > 1, true, "      Mismatch penalty: ",           bwaMismatch);
  logger.log(verbose > 1, true, "      Sequence match score: ",       bwaMatchScore);
  logger.log(verbose > 1, true, "      Z-dropoff: ",                  bwaZdrop);
  logger.log(verbose > 1, true, "      Alignment bandwidth: ",        bwaBandwidth);
  logger.log(verbose > 1, true, "      Clip 3 penalty: ",             bwaClip3);
  logger.log(verbose > 1, true, "      Clip 5 penalty: ",             bwaClip5);
  logger.log(verbose > 1, true, "      Reseed trigger: ",             bwaReseedTrigger);

  if (discClusterOnly) {
    logger.log(true, true, "    ######## ONLY DISCORDANT READ CLUSTERING. NO ASSEMBLY ##############");
  }
}
