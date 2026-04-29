cmake_minimum_required(VERSION 3.10)

# Set the C++ standard required for the project
project(svaba) # LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Default build type" FORCE)
endif()
add_compile_options(-fno-omit-frame-pointer)

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

## LZMA
find_path(LZMA_INCLUDE_DIR NAMES lzma.h)
find_library(LZMA_LIBRARY NAMES lzma)

if(NOT LZMA_INCLUDE_DIR OR NOT LZMA_LIBRARY)
  message(FATAL_ERROR "LZMA library or headers not found!")
endif()

## BZip2
find_package(BZip2 REQUIRED)

## HTSLIB
#
# Detection order:
#   1. find_package(htslib) — works only if htslib ships a
#      htslibConfig.cmake (rare outside conda/vcpkg installs).
#   2. pkg-config — system-installed htslib typically provides a
#      htslib.pc, so `pkg-config --libs htslib` works.
#   3. find_path/find_library for htslib/hts.h and libhts on the
#      default system paths (/usr/include, /usr/local/include,
#      /opt/homebrew/include, etc.).
#   4. Explicit -DHTSLIB_DIR=/path/to/htslib for non-standard installs.
#
# The first three paths auto-detect a system install, so a vanilla
# `cmake .. && make` works when htslib is installed system-wide. Only
# fall back to HTSLIB_DIR when all auto-detections fail, matching the
# pattern used here for ZLIB / BZip2 (also auto-detected).
set(SVABA_HTSLIB_FOUND FALSE)
set(SVABA_HTSLIB_LIBRARIES "")

if (htslib_FOUND)
    message(STATUS "htslib: found via find_package()")
    set(SVABA_HTSLIB_FOUND TRUE)
    set(SVABA_HTSLIB_LIBRARIES hts)
endif()

if (NOT SVABA_HTSLIB_FOUND)
    find_package(PkgConfig QUIET)
    if (PkgConfig_FOUND)
        pkg_check_modules(HTSLIB_PC QUIET htslib)
        if (HTSLIB_PC_FOUND)
            message(STATUS "htslib: found via pkg-config (version ${HTSLIB_PC_VERSION})")
            include_directories(${HTSLIB_PC_INCLUDE_DIRS})
            link_directories(${HTSLIB_PC_LIBRARY_DIRS})
            set(SVABA_HTSLIB_FOUND TRUE)
            set(SVABA_HTSLIB_LIBRARIES ${HTSLIB_PC_LIBRARIES})
        endif()
    endif()
endif()

if (NOT SVABA_HTSLIB_FOUND)
    find_path(HTSLIB_INCLUDE_DIR NAMES htslib/hts.h
              DOC "Directory containing htslib/hts.h")
    find_library(HTSLIB_LIBRARY NAMES hts
                 DOC "Path to the hts library")
    if (HTSLIB_INCLUDE_DIR AND HTSLIB_LIBRARY)
        message(STATUS "htslib: found headers at ${HTSLIB_INCLUDE_DIR}, lib at ${HTSLIB_LIBRARY}")
        include_directories(${HTSLIB_INCLUDE_DIR})
        set(SVABA_HTSLIB_FOUND TRUE)
        set(SVABA_HTSLIB_LIBRARIES ${HTSLIB_LIBRARY})
    endif()
endif()

if (NOT SVABA_HTSLIB_FOUND)
    set(HTSLIB_DIR "" CACHE PATH "Path to HTSLib root directory")
    if (HTSLIB_DIR)
        message(STATUS "htslib: using explicit HTSLIB_DIR=${HTSLIB_DIR}")
        include_directories(${HTSLIB_DIR}/include)
        link_directories(${HTSLIB_DIR}/lib)
        set(SVABA_HTSLIB_FOUND TRUE)
        set(SVABA_HTSLIB_LIBRARIES hts)
    else()
        message(FATAL_ERROR
            "htslib not found. Install htslib system-wide (e.g. "
            "`brew install htslib` or `apt install libhts-dev`) so "
            "pkg-config / find_library locate it automatically, "
            "or pass -DHTSLIB_DIR=/path/to/htslib to point at a "
            "manual build.")
    endif()
endif()

# Find all source files
set(SOURCES
src/svaba/run_svaba.cpp
src/svaba/SvabaOutputWriter.cpp
src/svaba/BreakPoint.cpp
src/svaba/ContigAlignmentScore.cpp
src/svaba/AlignedContig.cpp
src/svaba/AlignmentFragment.cpp
src/svaba/DiscordantCluster.cpp
src/svaba/DBSnpFilter.cpp
src/svaba/SvabaUtils.cpp
src/svaba/svaba.cpp
src/svaba/SvabaAssemblerEngine.cpp
src/svaba/tovcf.cpp
src/svaba/test_svaba.cpp
src/svaba/vcf.cpp
src/svaba/DiscordantRealigner.cpp
src/svaba/SvabaOverlapAlgorithm.cpp
src/svaba/SvabaASQG.cpp
src/svaba/SvabaAssemble.cpp
src/svaba/KmerFilter.cpp
src/svaba/SvabaThreadUnit.cpp
src/svaba/SvabaRegionProcessor.cpp
src/svaba/SvabaBamWalker.cpp
src/svaba/refilter.cpp
src/svaba/SvabaPostprocess.cpp
src/svaba/LearnBamParams.cpp
src/svaba/STCoverage.cpp
src/svaba/SvabaModels.cpp
##src/svaba/Histogram.cpp
##src/svaba/BamStats.cpp
src/svaba/SvabaRead.cpp
src/svaba/SvabaOptions.cpp
src/svaba/SvabaLogger.cpp
src/svaba/SvabaFileLoader.cpp

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

# ---------------------------------------------------------------------------
# Optional compile-time contig tracing (SvabaDebug.h).
#
#   cmake .. -DSVABA_TRACE_CONTIG="c_fermi_chr2_215869501_215894501_11C"
#
# Prints detailed stderr trace for every decision point that contig
# touches. Pass -DSVABA_TRACE_ALL=ON to trace ALL contigs (very noisy).
# ---------------------------------------------------------------------------
set(SVABA_TRACE_CONTIG "" CACHE STRING "Contig name to trace (empty=off)")
option(SVABA_TRACE_ALL "Trace ALL contigs (very noisy)" OFF)
set(SVABA_TRACE_READ "" CACHE STRING "Read QNAME to trace (empty=off)")
option(SVABA_TRACE_ALL_READS "Trace ALL reads (extremely noisy)" OFF)

# Generate the executable
add_executable(svaba ${SOURCES})

if (NOT "${SVABA_TRACE_CONTIG}" STREQUAL "")
    message(STATUS "svaba: tracing contig '${SVABA_TRACE_CONTIG}'")
    target_compile_definitions(svaba PRIVATE
        SVABA_TRACE_CONTIG="${SVABA_TRACE_CONTIG}")
endif()
if (SVABA_TRACE_ALL)
    message(STATUS "svaba: tracing ALL contigs (noisy!)")
    target_compile_definitions(svaba PRIVATE SVABA_TRACE_ALL=1)
endif()
if (NOT "${SVABA_TRACE_READ}" STREQUAL "")
    message(STATUS "svaba: tracing read '${SVABA_TRACE_READ}'")
    target_compile_definitions(svaba PRIVATE
        SVABA_TRACE_READ="${SVABA_TRACE_READ}")
endif()
if (SVABA_TRACE_ALL_READS)
    message(STATUS "svaba: tracing ALL reads (extremely noisy!)")
    target_compile_definitions(svaba PRIVATE SVABA_TRACE_ALL_READS=1)
endif()

# If the submodules already have Makefiles, you can use custom commands to invoke make
# in those directories.
add_custom_target(
    COMMAND make
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# ---------------------------------------------------------------------------
# Optional jemalloc support.
#
# On Linux with high thread counts (-p 16+), jemalloc eliminates glibc
# arena-lock contention that can cost 30-40% wall time. Link it at
# compile time so no LD_PRELOAD is needed at runtime.
#
#   cmake .. -DUSE_JEMALLOC=ON      # link jemalloc (recommended on Linux)
#   cmake .. -DUSE_JEMALLOC=OFF     # use system malloc (default)
#
# macOS users should NOT enable this — Apple's libmalloc with its
# nanomalloc fast path outperforms jemalloc on this workload.
# ---------------------------------------------------------------------------
option(USE_JEMALLOC "Link against jemalloc for reduced allocator contention" OFF)
set(JEMALLOC_LIBRARIES "")
if (USE_JEMALLOC)
    find_library(JEMALLOC_LIB NAMES jemalloc
                 PATHS /usr/lib/x86_64-linux-gnu /usr/lib /usr/local/lib
                 DOC "Path to libjemalloc")
    if (JEMALLOC_LIB)
        message(STATUS "jemalloc: found at ${JEMALLOC_LIB}")
        set(JEMALLOC_LIBRARIES ${JEMALLOC_LIB})
    else()
        message(FATAL_ERROR
            "USE_JEMALLOC=ON but libjemalloc not found. "
            "Install it (apt install libjemalloc-dev) or set "
            "-DJEMALLOC_LIB=/path/to/libjemalloc.so")
    endif()
endif()

# Linking
target_link_libraries(svaba
    seqlib
    ${CMAKE_CURRENT_SOURCE_DIR}/SeqLib/bwa/libbwa.a
    ${CMAKE_CURRENT_SOURCE_DIR}/SeqLib/fermi-lite/libfml.a
    Threads::Threads
    ${SVABA_HTSLIB_LIBRARIES}
    ZLIB::ZLIB
    ${CURL_LIBRARY}
    ${LZMA_LIBRARY}
    BZip2::BZip2
    ${JEMALLOC_LIBRARIES}
)

# ---------------------------------------------------------------------------
# Install support.
#
# `make install` (equivalent: `cmake --install build`) copies the svaba
# binary to ${CMAKE_INSTALL_PREFIX}/bin. GNUInstallDirs sets the usual
# distro-friendly defaults for BINDIR / LIBDIR / etc., so this follows
# standard Unix conventions.
#
# Default CMAKE_INSTALL_PREFIX is /usr/local (may need sudo). To stage
# into a non-system location, pass -DCMAKE_INSTALL_PREFIX=/some/path at
# cmake time, or use `cmake --install build --prefix /some/path` at
# install time. To drop the binary into this repo's `bin/` directory:
#   cmake --install build --prefix ${CMAKE_SOURCE_DIR}
# ---------------------------------------------------------------------------
include(GNUInstallDirs)
install(TARGETS svaba
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(PROGRAMS ${CMAKE_SOURCE_DIR}/scripts/svaba_postprocess.sh
        DESTINATION ${CMAKE_INSTALL_BINDIR})
