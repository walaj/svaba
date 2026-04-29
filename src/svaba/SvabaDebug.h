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

// ── Kmer restriction ────────────────────────────────────────────────
//
// Compile-time read filter: after BFC error correction, only reads
// whose corrected sequence contains the specified kmer survive into
// assembly/r2c/corrected.bam. All others get to_assemble = false.
//
// Usage:
//   cmake .. -DCMAKE_CXX_FLAGS='-DSVABA_KMER_RESTRICT="\"CCATGCAGAGTGTTGAAGAAAAGGC\""'
//
// Also checks the reverse complement so orientation doesn't matter.
// Zero cost when not compiled in.

#ifdef SVABA_KMER_RESTRICT

  inline std::string _svaba_revcomp(const std::string& s) {
    std::string rc(s.size(), 'N');
    for (size_t i = 0; i < s.size(); ++i) {
      switch (s[s.size() - 1 - i]) {
        case 'A': case 'a': rc[i] = 'T'; break;
        case 'T': case 't': rc[i] = 'A'; break;
        case 'C': case 'c': rc[i] = 'G'; break;
        case 'G': case 'g': rc[i] = 'C'; break;
        default:            rc[i] = 'N'; break;
      }
    }
    return rc;
  }

  inline bool _svaba_kmer_match(const std::string& seq) {
    static const std::string kmer_fwd(SVABA_KMER_RESTRICT);
    static const std::string kmer_rc = _svaba_revcomp(kmer_fwd);
    return seq.find(kmer_fwd) != std::string::npos ||
           seq.find(kmer_rc)  != std::string::npos;
  }

  #define SVABA_KMER_RESTRICTING 1

#endif
