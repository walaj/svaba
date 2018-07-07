#ifndef SVABA_PARAMS_H__
#define SVABA_PARAMS_H__

// moved from AlignmentFragment.h
/////////////////////////////////
#define MAX_CONTIG_SIZE 5000000

// moved from BreakPoint
////////////////////////
#define MAX_ERROR 0.04
#define MIN_ERROR 0.0005

#define T_SPLIT_BUFF 5
#define N_SPLIT_BUFF 5

// if the insertion is this big or larger, don't require splits to span both sides
#define INSERT_SIZE_TOO_BIG_SPAN_READS 16

// if homology is greater than homology / HOMOLOGY_FACTOR, then reject for assembly-only
#define HOMOLOGY_FACTOR 4
#define MIN_SOMATIC_RATIO 15

// when calculating coverage at one base, average over bases (left and right)
#define COVERAGE_AVG_BUFF 10

// moved from DiscordantCluster
///////////////////////////////
#define DISC_PAD 150
#define MIN_PER_CLUST 2
#define DEFAULT_ISIZE_THRESHOLD 800 // shouldn't be hit if isize was learned

// moved from run_svaba
///////////////////////
#define THREAD_READ_LIMIT 20000
#define THREAD_CONTIG_LIMIT 250

// minimum number of reads to support even reporting dscrd cluster 
// (if not assocaited with assembly contig)
#define MIN_DSCRD_READS_DSCRD_ONLY 3 

// moved from svabaAssemblerEngine
//////////////////////////////////
#define MAX_OVERLAPS_PER_ASSEMBLY 20000

#define MIN_CONTIG_MATCH 35
#define MATE_LOOKUP_MIN 3
#define SECONDARY_CAP 10
#define MAX_MATE_ROUNDS 1
#define MATE_REGION_LOOKUP_LIMIT 400
#define MAX_NUM_MATE_WINDOWS 50000000

#define GERMLINE_CNV_PAD 10
#define WINDOW_PAD 500
#define MICROBE_MATCH_MIN 50
#define GET_MATES 1
#define MICROBE 1
#define LARGE_INTRA_LOOKUP_LIMIT 50000
#define SECONDARY_FRAC 0.90

// moved from svabaBamWalker
////////////////////////////
//#define DEBUG_SVABA_BAMWALKER 1
#define MIN_MAPQ_FOR_MATE_LOOKUP 0
//#define TRAIN_READS_FAIL_SAFE 50000
#define MIN_ISIZE_FOR_DISCORDANT_REALIGNMENT 1000
#define DISC_REALIGN_MATE_PAD 100
#define MAX_SECONDARY_HIT_DISC 10
#define MATE_REGION_PAD 250

// trim this many bases from front and back of read when determining coverage
// this should be synced with the split-read buffer in BreakPoint2 for more accurate 
// representation of covearge of INFORMATIVE reads (eg ones that could be split)
#define INFORMATIVE_COVERAGE_BUFFER 0

// moved from vcf
/////////////////
#define VCF_SECONDARY_CAP 200
#define SOMATIC_LOD 1 // just a dummy now. scoring is elsewhere, and output is 0 (germline) or 1 (somatic)
#define DEDUPEPAD 200

#endif
