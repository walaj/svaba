#include "SvabaOptions.h"

#include <sstream>           // for std::ostringstream
#include <iomanip>           // for std::setw and std::setfill
#include <getopt.h>
#include <iostream>
#include <stdexcept>

#include "SvabaLogger.h"
#include "SvabaAssemblerConfig.h"  // SVABA_ASSEMBLER_FERMI for SGA guards

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
)"
#if !SVABA_ASSEMBLER_FERMI
  // SGA-specific knobs. Only surface them to --help when svaba was
  // compiled with SGA as the active assembler (see
  // SvabaAssemblerConfig.h). The option codes themselves still parse
  // under fermi — they just get ignored.
            << R"(
Assembly (SGA):
      --min-overlap <bp>  Minimum read overlap for SGA
      --error-rate <f>    SGA fractional error rate
      --rounds <N>        Number of assembly rounds
)"
#endif
            << R"(
Error-correction:
      --ec-type <s|f|0>   s=SGA k-mer; f=Fermi BFC; 0=off
      --ec-subsample <f>  Fraction to sample for EC learning

Discordant clustering:
      --disc-sd <f>       SD cutoff, default 3.92

Filtering:
      --max-cov <N>       Max coverage to assemble, default 100
      --mate-min <N>      Min reads to trigger somatic mate lookup, default 3
      --mate-min-count <N>
                          Min reads to form a candidate mate region, default 2
      --mate-lim <N>      Max reads in mate lookup, default 400
      --min-mate-mapq <N> Min MAPQ on primary read for mate-region candidate,
                          default -1 (no gate). Set to e.g. 1 to exclude
                          MAPQ=0 multi-mappers from mate lookup.
      --max-mate-chr <N>  Exclude mate regions on chromosomes with ChrID > N,
                          default 23 (through chrY; skips chrM/alt/decoy).
                          Set to -1 for no limit (see --non-human).
      --non-human         Remove hardcoded human genome assumptions.
                          Currently sets --max-mate-chr -1 (allow all chroms)
                          and samples all contigs during insert-size learning.
      --no-nm             Skip high-NM read salvage (NM/len > 0.02).
                          Faster: fewer reads enter r2c and correction.
                          Trades away rare NM-only SV sensitivity.

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
      --dump-reads        Emit per-read debug/visualization outputs
                          (off by default; large output on deep samples).
                          Four files are produced by this flag together:
                            ${ID}.corrected.bam
                            ${ID}.discordant.bam
                            ${ID}.alignments.txt.gz
                            ${ID}.r2c.txt.gz
                          Without this flag, svaba writes the bps.txt.gz
                          / VCF / contigs.bam / runtime.txt summaries but
                          none of the per-read detail. Weird-reads BAM is
                          compile-time only; see
                          SvabaOptions.h::dump_weird_reads.
      --always-realign-corrected
                          Force re-alignment of every corrected read to
                          the reference, even if BFC didn't modify it.
                          By default, reads whose sequence is unchanged
                          reuse the input BAM's CIGAR/NM (valid when the
                          BAM was aligned with BWA). Use this flag when
                          the input was aligned with a non-BWA aligner.
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
    {"no-nm",      no_argument,            nullptr,  1403},
    {"min-mate-mapq",   required_argument, nullptr,  1404},
    {"max-mate-chr",    required_argument, nullptr,  1405},
    {"non-human",       no_argument,       nullptr,  1406},
    {"mate-min-count",  required_argument, nullptr,  1407},
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
    // --dump-reads: runtime opt-in for corrected + discordant-reads BAMs.
    // Single flag, both side effects — see comment in SvabaOptions.h.
    // Weird-reads BAM is deliberately not on this flag (compile-time only).
    {"dump-reads", no_argument,            nullptr,  1800},
    {"always-realign-corrected", no_argument, nullptr, 1801},
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
      case 1403: o.noNmSalvage          = true; break;
      case 1404: o.minMateMAPQ          = std::stoi(optarg); break;
      case 1405: o.maxMateChrID         = std::stoi(optarg); break;
      case 1406: o.nonHuman             = true; break;
      case 1407: o.mateRegionMinCount   = std::stoi(optarg); break;

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

      case 1700: o.blacklistFile.push_back(optarg); break;
      case 1701: o.germlineSvFile = optarg; break;
      case 1702: o.dbsnpVcf       = optarg; break;

      // --dump-reads is the single runtime knob for all per-read detail
      // outputs. Flip all three flags together here; any call-site
      // reading these fields after parse() sees a consistent state.
      // See SvabaOptions.h for the comment on why these are separate
      // fields instead of a single bool.
      case 1800:
        o.dump_discordant_reads = true;
        o.dump_corrected_reads  = true;
        o.dump_alignments       = true;
        break;

      case 1801:
        o.alwaysRealignCorrected = true;
        break;

      case '?':
      default:
        throw std::runtime_error("Unknown or malformed option; see --help");
    }
  }

  // --non-human: remove human-specific assumptions
  if (o.nonHuman) {
    o.maxMateChrID = -1;  // allow all chromosomes for mate lookup
  }

  // post validation
  if (argc == 1) {
    o.help = true;
  } else if (!o.help) {
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
  //  if (o.bams["t001"] == "-" && !o.rulesJson.empty())
  o.rulesJson = R"json(
{
  "global": {
    "duplicate": false,
    "qcfail":   false
  },
  "": {
    "rules": [
      { "rr":    true },
      { "ff":    true },
      { "rf":    true },
      { "ic":    true },
      { "clip":  5, "length": 30 },
      { "ins":   true },
      { "del":   true },
      { "mapped":      true, "mate_mapped": false },
      { "mate_mapped": true, "mapped":      false }
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
  logger.log(verbose > 1, true, "    Assembler: ", svaba::kAssemblerName);
#if !SVABA_ASSEMBLER_FERMI
  // SGA-only parameters; don't confuse fermi-build operators with
  // values they can't tune without recompiling.
  logger.log(verbose > 1, true, "    ErrorRate (SGA): ",
             (sgaErrorRate < 0.001f ? std::string("EXACT (0)")
                                        : std::to_string(sgaErrorRate)));
  logger.log(verbose > 1, true, "    Num assembly rounds (SGA): ",     sgaNumRounds);
#endif
  logger.log(verbose > 1, true, "    Discordant SD cutoff (tumor): ",  sdDiscCutoff);
  logger.log(verbose > 1, true, "    Discordant SD cutoff (normal): ", sdDiscCutoffNormal);
  logger.log(verbose > 1, true, "    Minimum number of reads for mate lookup: ", mateLookupMin);
  logger.log(verbose > 1, true, "    Min reads to form candidate mate region: ", mateRegionMinCount);
  logger.log(verbose > 1, true, "    Min MAPQ for mate region candidate: ", minMateMAPQ, (minMateMAPQ < 0 ? " (no gate)" : ""));
  logger.log(verbose > 1, true, "    Max ChrID for mate region: ", maxMateChrID, (maxMateChrID < 0 ? " (no limit)" : ""));
  if (nonHuman) logger.log(true, true, "    Non-human genome mode: chromosome gates disabled");
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
