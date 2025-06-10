// svaba_params.h
#pragma once

#include <cstddef>

inline constexpr int PER_THREAD_BATCH_SIZE = 1000;

// version & date
inline constexpr char SVABA_VERSION[] = "1.3.0";
inline constexpr char SVABA_DATE[]    = "05/2025";

// from AlignmentFragment.h
inline constexpr std::size_t MAX_CONTIG_SIZE = 5'000'000;

// from run_svaba.cpp
inline constexpr int MIN_CLIP_FOR_LOCAL           = 40;
inline constexpr int MAX_NM_FOR_LOCAL             = 10;

// from BreakPoint
inline constexpr double MAX_ERROR                   = 0.04;
inline constexpr double MIN_ERROR                   = 0.0005;
inline constexpr int    T_SPLIT_BUFF                = 5;
inline constexpr int    N_SPLIT_BUFF                = 5;
inline constexpr int    INSERT_SIZE_TOO_BIG_SPAN_READS = 16;
inline constexpr int    HOMOLOGY_FACTOR             = 4;
inline constexpr int    MIN_SOMATIC_RATIO           = 15;
inline constexpr int    COVERAGE_AVG_BUFF           = 10;

// from DiscordantCluster
inline constexpr int DISC_PAD                 = 150;
inline constexpr int MIN_PER_CLUST            = 2;
inline constexpr int DEFAULT_ISIZE_THRESHOLD  = 2000;

// from run_svaba
inline constexpr std::size_t THREAD_READ_LIMIT      = 20'000; 
inline constexpr int         THREAD_CONTIG_LIMIT    =   5'000;
inline constexpr int         MIN_DSCRD_READS_DSCRD_ONLY = 3;

// from svabaAssemblerEngine
inline constexpr std::size_t MAX_OVERLAPS_PER_ASSEMBLY = 20'000;
inline constexpr int         MIN_CONTIG_MATCH           =    35;
inline constexpr int         MATE_LOOKUP_MIN            =     3;
inline constexpr int         SECONDARY_CAP              =    10;
inline constexpr int         MAX_MATE_ROUNDS            =     1;
inline constexpr std::size_t MATE_REGION_LOOKUP_LIMIT  =   400;
inline constexpr std::size_t MAX_NUM_MATE_WINDOWS      = 50'000'000;
inline constexpr int         GERMLINE_CNV_PAD           =    10;
inline constexpr int         GET_MATES                  =     1;
inline constexpr int         LARGE_INTRA_LOOKUP_LIMIT   = 50'000;
inline constexpr double      SECONDARY_FRAC             =  0.90;

// from svabaBamWalker
inline constexpr int MIN_MAPQ_FOR_MATE_LOOKUP            =     0;
inline constexpr int MIN_ISIZE_FOR_DISCORDANT_REALIGNMENT = 1'000;
inline constexpr int DISC_REALIGN_MATE_PAD                =   100;
inline constexpr int MAX_SECONDARY_HIT_DISC               =    10;
inline constexpr int MATE_REGION_PAD                      =   250;

// coverage buffer
inline constexpr int INFORMATIVE_COVERAGE_BUFFER = 0;

// from vcf
inline constexpr int VCF_SECONDARY_CAP = 200;
inline constexpr int SOMATIC_LOD       =   1;
inline constexpr int DEDUPEPAD         = 200;

