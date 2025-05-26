cmake_minimum_required(VERSION 3.10)

# Set the C++ standard required for the project
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set O2 optimization
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

set(CMAKE_BUILD_TYPE Debug)

project(svaba)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/SeqLib)

# Include directories for headers
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/SeqLib
    ${CMAKE_CURRENT_SOURCE_DIR}/SeqLib/bwa
    ${CMAKE_CURRENT_SOURCE_DIR}/SeqLib/fermi-lite    
)

# Look for htslib on the system
find_package(htslib QUIET)

# Find required system level type libraries
find_package(Threads REQUIRED)  
find_package(ZLIB REQUIRED)

## libcurl
# Find libcurl or allow user-specified path
find_package(CURL)

if(CURL_FOUND)
    message(STATUS "Found CURL: ${CURL_LIBRARIES}")
else()
    message(WARNING "CURL not found. You can specify the path with -DCURL_LIBRARY=<path_to_curl> and -DCURL_INCLUDE_DIR=<path_to_headers>.")
    set(CURL_LIBRARY CACHE STRING "Path to libcurl library")
    set(CURL_INCLUDE_DIR CACHE STRING "Path to libcurl include directory")
endif()

## LZMA
find_path(LZMA_INCLUDE_DIR NAMES lzma.h)
find_library(LZMA_LIBRARY NAMES lzma)

if(NOT LZMA_INCLUDE_DIR OR NOT LZMA_LIBRARY)
  message(FATAL_ERROR "LZMA library or headers not found!")
endif()

## BZip2
find_package(BZip2 REQUIRED)

## HTSLIB
if (htslib_FOUND)
    # If htslib was found on the system, use it
    message(STATUS "Using system htslib")
else()
    set(HTSLIB_DIR "" CACHE PATH "Path to HTSLib root directory")

    if (NOT HTSLIB_DIR)
        message(FATAL_ERROR "HTSLIB_DIR not specified. Please specify -DHTSLIB_DIR=/path/to/htslib")
    else()
        # Automatically set include and library paths based on HTSLIB_DIR
        include_directories(${HTSLIB_DIR}/include)
        link_directories(${HTSLIB_DIR}/lib)
    endif()
endif()

# Find all source files
set(SOURCES
src/svaba/run_svaba.cpp
src/svaba/svabaOutputWriter.cpp
src/svaba/BreakPoint.cpp
src/svaba/AlignedContig.cpp
src/svaba/AlignmentFragment.cpp
src/svaba/DiscordantCluster.cpp
src/svaba/DBSnpFilter.cpp
src/svaba/svabaUtils.cpp
src/svaba/svaba.cpp
src/svaba/svabaAssemblerEngine.cpp
src/svaba/tovcf.cpp
src/svaba/vcf.cpp
src/svaba/DiscordantRealigner.cpp
src/svaba/svabaOverlapAlgorithm.cpp
src/svaba/svabaASQG.cpp
src/svaba/svabaAssemble.cpp
src/svaba/KmerFilter.cpp
src/svaba/svabaBamWalker.cpp
src/svaba/refilter.cpp
src/svaba/LearnBamParams.cpp
src/svaba/STCoverage.cpp
src/svaba/Histogram.cpp
src/svaba/BamStats.cpp
src/svaba/svabaRead.cpp

src/SGA/SuffixTools/STCommon.cpp
src/SGA/SuffixTools/Occurrence.cpp
src/SGA/SuffixTools/SuffixArray.cpp
src/SGA/SuffixTools/SuffixCompare.cpp
src/SGA/SuffixTools/InverseSuffixArray.cpp
src/SGA/SuffixTools/SACAInducedCopying.cpp
src/SGA/SuffixTools/BWTAlgorithms.cpp
src/SGA/SuffixTools/BWTReader.cpp
src/SGA/SuffixTools/BWTWriter.cpp
src/SGA/SuffixTools/SAReader.cpp
src/SGA/SuffixTools/SAWriter.cpp
src/SGA/SuffixTools/SBWT.cpp
src/SGA/SuffixTools/RLBWT.cpp
src/SGA/SuffixTools/BWTWriterBinary.cpp
src/SGA/SuffixTools/BWTReaderBinary.cpp
src/SGA/SuffixTools/BWTWriterAscii.cpp
src/SGA/SuffixTools/BWTReaderAscii.cpp
src/SGA/SuffixTools/BWTIntervalCache.cpp
src/SGA/SuffixTools/SampledSuffixArray.cpp

src/SGA/Algorithm/OverlapAlgorithm.cpp
src/SGA/Algorithm/DPAlignment.cpp
src/SGA/Algorithm/SearchSeed.cpp
src/SGA/Algorithm/OverlapBlock.cpp
src/SGA/Algorithm/SearchHistory.cpp
src/SGA/Algorithm/OverlapTools.cpp

src/SGA/Bigraph/Bigraph.cpp
src/SGA/Bigraph/Vertex.cpp
src/SGA/Bigraph/Edge.cpp
src/SGA/Bigraph/EdgeDesc.cpp

src/SGA/SGA/OverlapCommon.cpp

src/SGA/SQG/SQG.cpp
src/SGA/SQG/ASQG.cpp

src/SGA/StringGraph/SGUtil.cpp
src/SGA/StringGraph/SGAlgorithms.cpp
src/SGA/StringGraph/SGVisitors.cpp
src/SGA/StringGraph/CompleteOverlapSet.cpp
src/SGA/StringGraph/RemovalAlgorithm.cpp
src/SGA/StringGraph/SGSearch.cpp
src/SGA/StringGraph/SGWalk.cpp

src/SGA/Util/Util.cpp
src/SGA/Util/stdaln.c
src/SGA/Util/Alphabet.cpp
src/SGA/Util/Contig.cpp
src/SGA/Util/ReadTable.cpp
src/SGA/Util/ReadInfoTable.cpp
src/SGA/Util/SeqReader.cpp
src/SGA/Util/DNAString.cpp
src/SGA/Util/Match.cpp
src/SGA/Util/Pileup.cpp
src/SGA/Util/Interval.cpp
src/SGA/Util/SeqCoord.cpp
src/SGA/Util/QualityVector.cpp
src/SGA/Util/Quality.cpp
src/SGA/Util/PrimerScreen.cpp
src/SGA/Util/CorrectionThresholds.cpp
src/SGA/Util/ClusterReader.cpp
src/SGA/Util/QualityTable.cpp
src/SGA/Util/gzstream.C
src/SGA/Util/BitChar.cpp
src/SGA/Util/MultiOverlap.cpp
)

include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}
    src/SGA/Algorithm
    src/SGA/StringGraph
    src/SGA/SGA
    src/SGA/SuffixTools
    src/SGA/Util
    src/SGA/SQG
    src/SGA/Bigraph
    src/svaba
    src/SGA
)

# Generate the executable
add_executable(svaba ${SOURCES})

# If the submodules already have Makefiles, you can use custom commands to invoke make
# in those directories.
add_custom_target(
    COMMAND make
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# Linking
target_link_libraries(svaba 
    seqlib
    ${CMAKE_CURRENT_SOURCE_DIR}/SeqLib/bwa/libbwa.a
    ${CMAKE_CURRENT_SOURCE_DIR}/SeqLib/fermi-lite/libfml.a
    Threads::Threads
    hts
    ZLIB::ZLIB
    ${CURL_LIBRARY}
    ${LZMA_LIBRARY}
    BZip2::BZip2
)
