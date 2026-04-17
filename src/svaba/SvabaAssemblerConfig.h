#pragma once

// ---------------------------------------------------------------------------
// svabaAssemblerConfig.h
//
// Compile-time selection of the local-assembly engine used by svaba.
//
//   SVABA_ASSEMBLER_FERMI   (1) -> use fermi2 (ropebwt2-based) assembler
//   SVABA_ASSEMBLER_FERMI   (0) -> use the vendored SGA String Graph Assembler
//
// This used to live as a bare `#define FERMI 1` at the top of
// SvabaRegionProcessor.cpp, which meant nothing else in the build could see
// it (so e.g. run_svaba.cpp could not report which assembler was compiled
// in). Centralizing it here lets any TU `#include "SvabaAssemblerConfig.h"`
// and branch on the same symbol.
//
// To switch assemblers at build time without editing this header, pass
//   -DSVABA_ASSEMBLER_FERMI=0
// (or =1) via CMake / CXXFLAGS.
// ---------------------------------------------------------------------------

#ifndef SVABA_ASSEMBLER_FERMI
#define SVABA_ASSEMBLER_FERMI 1
#endif

// Back-compat alias: existing code uses `#ifndef FERMI` / `#ifdef FERMI`.
#if SVABA_ASSEMBLER_FERMI
  #ifndef FERMI
    #define FERMI 1
  #endif
#else
  #ifdef FERMI
    #undef FERMI
  #endif
#endif

// Human-readable name of the active assembler, for logging.
namespace svaba {
#if SVABA_ASSEMBLER_FERMI
  inline constexpr const char* kAssemblerName = "fermi2";
  inline constexpr const char* kAssemblerTag  = "fermi";
#else
  inline constexpr const char* kAssemblerName = "SGA";
  inline constexpr const char* kAssemblerTag  = "sga";
#endif
}
