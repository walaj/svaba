// SvabaDebug.h — compile-time contig tracing.
//
// Usage:
//   cmake .. -DSVABA_TRACE_CONTIG='"c_fermi_chr2_215869501_215894501_11C"'
//   make -j
//
// This compiles svaba with verbose stderr tracing for every decision
// point that contig touches: assembly filtering, BP identification,
// r2c scoring, split coverage, confidence assignment, and output gating.
//
// To disable (default): just don't pass -DSVABA_TRACE_CONTIG.
// To trace ALL contigs (very noisy): -DSVABA_TRACE_ALL=1
//
#pragma once

#include <iostream>
#include <string>

#ifdef SVABA_TRACE_CONTIG

  // Match by contig name (cname)
  inline bool _svaba_trace_match(const std::string& cname) {
    return cname == SVABA_TRACE_CONTIG;
  }

  #define SVABA_TRACE(cname, msg) \
    do { if (_svaba_trace_match(cname)) { \
      std::cerr << "[TRACE:" << __FILE__ << ":" << __LINE__ << "] " \
                << msg << std::endl; \
    } } while(0)

  #define SVABA_TRACE_IF(cond, cname, msg) \
    do { if ((cond) && _svaba_trace_match(cname)) { \
      std::cerr << "[TRACE:" << __FILE__ << ":" << __LINE__ << "] " \
                << msg << std::endl; \
    } } while(0)

  #define SVABA_TRACING 1

#elif defined(SVABA_TRACE_ALL)

  #define SVABA_TRACE(cname, msg) \
    do { std::cerr << "[TRACE:" << __FILE__ << ":" << __LINE__ << "] " \
                   << "[" << (cname) << "] " << msg << std::endl; \
    } while(0)

  #define SVABA_TRACE_IF(cond, cname, msg) \
    do { if (cond) { std::cerr << "[TRACE:" << __FILE__ << ":" << __LINE__ << "] " \
                               << "[" << (cname) << "] " << msg << std::endl; } \
    } while(0)

  #define SVABA_TRACING 1

#else

  #define SVABA_TRACE(cname, msg) do {} while(0)
  #define SVABA_TRACE_IF(cond, cname, msg) do {} while(0)
  // SVABA_TRACING intentionally not defined

#endif
