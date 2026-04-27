// SvabaDebug.h — compile-time contig and read tracing.
//
// Contig tracing:
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
// Read tracing:
//   cmake .. -DSVABA_TRACE_READ='"LH00306:235:22NHGWLT4:6:1329:45742:15246"'
//   make -j
//
// Traces a specific read (by QNAME) through the BamWalker read-filter
// pipeline: duplicate/QC/blacklist gates, rule_pass, high-NM salvage,
// adapter filter, quality trim, and final buffer admission.
//
// To trace ALL reads (extremely noisy): -DSVABA_TRACE_ALL_READS=1
//
#pragma once

#include <iostream>
#include <string>

// ── Contig tracing ──────────────────────────────────────────────────

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

// ── Read tracing ────────────────────────────────────────────────────

#ifdef SVABA_TRACE_READ

  inline bool _svaba_read_trace_match(const std::string& qname) {
    return qname == SVABA_TRACE_READ;
  }

  #define SVABA_READ_TRACE(qname, msg) \
    do { if (_svaba_read_trace_match(qname)) { \
      std::cerr << "[READ_TRACE:" << __FILE__ << ":" << __LINE__ << "] " \
                << msg << std::endl; \
    } } while(0)

  #define SVABA_READ_TRACE_IF(cond, qname, msg) \
    do { if ((cond) && _svaba_read_trace_match(qname)) { \
      std::cerr << "[READ_TRACE:" << __FILE__ << ":" << __LINE__ << "] " \
                << msg << std::endl; \
    } } while(0)

  #define SVABA_READ_TRACING 1

#elif defined(SVABA_TRACE_ALL_READS)

  #define SVABA_READ_TRACE(qname, msg) \
    do { std::cerr << "[READ_TRACE:" << __FILE__ << ":" << __LINE__ << "] " \
                   << "[" << (qname) << "] " << msg << std::endl; \
    } while(0)

  #define SVABA_READ_TRACE_IF(cond, qname, msg) \
    do { if (cond) { std::cerr << "[READ_TRACE:" << __FILE__ << ":" << __LINE__ << "] " \
                               << "[" << (qname) << "] " << msg << std::endl; } \
    } while(0)

  #define SVABA_READ_TRACING 1

#else

  #define SVABA_READ_TRACE(qname, msg) do {} while(0)
  #define SVABA_READ_TRACE_IF(cond, qname, msg) do {} while(0)
  // SVABA_READ_TRACING intentionally not defined

#endif
