#!/usr/bin/env bash
#
# build_asan.sh — build svaba with AddressSanitizer (ASan) to pinpoint the
# teardown heap crash ("final flush" abort / segfault).
#
# Why this script exists: svaba's build has THREE compilation domains and
# ASan only helps where the code is instrumented:
#   1. svaba src + seqlib   -> CMake-built, picks up CMAKE_*_FLAGS below.
#   2. bwa + fermi-lite      -> PREBUILT static libs (SeqLib/bwa/libbwa.a,
#                               SeqLib/fermi-lite/libfml.a). Their Makefiles
#                               hardcode -O2 and ignore CMake, so they must
#                               be rebuilt explicitly with ASan CFLAGS.
#   3. htslib                -> external/precompiled. Left as-is. ASan still
#                               catches corruption that ORIGINATES in the
#                               instrumented code above, which is where the
#                               bug almost certainly lives.
#
# This is a THROWAWAY DEBUG build:
#   * It OVERWRITES the in-tree libbwa.a / libfml.a with ASan versions.
#     Restore production libs with the commands printed at the end.
#   * jemalloc MUST be OFF — ASan provides its own allocator; linking
#     jemalloc alongside ASan crashes at startup. This script forces it off.
#     (With jemalloc off, the malloc_trim guard leaves the trim calls in,
#     but ASan intercepts malloc_trim as a no-op, so that's harmless here.)
#
# Usage:
#   ./scripts/build_asan.sh
#
# Then run on the SAME failing input, same thread count, in a loop until
# it trips (ASan perturbs timing, so the ~10% rate may drop — keep looping):
#
#   ASAN_OPTIONS=detect_leaks=0:abort_on_error=1:halt_on_error=1:fast_unwind_on_malloc=0 \
#     ./build-asan/svaba run <same args, e.g. -t tumor.bam -n normal.bam -G ref.fa -p 16 ...> \
#     2> asan.log
#
# Send back asan.log — the first "ERROR: AddressSanitizer" block names the
# offending file:line plus the allocation stack. That ends the guessing.
#
set -euo pipefail

# Resolve repo root as the parent of this script's directory (portable).
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

JOBS="$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)"

# -O1 keeps ASan fast enough to reproduce while preserving readable stacks;
# -fno-omit-frame-pointer makes the backtraces trustworthy.
SUBMOD_CFLAGS="-g -O1 -fsanitize=address -fno-omit-frame-pointer -Wall -Wno-unused-function"
CMAKE_SAN_FLAGS="-fsanitize=address -fno-omit-frame-pointer"

echo ">>> [1/3] Rebuilding vendored bwa + fermi-lite with ASan (overwrites in-tree .a)"
make -C SeqLib/bwa        clean
make -C SeqLib/fermi-lite clean
make -C SeqLib/bwa        -j"$JOBS" CFLAGS="$SUBMOD_CFLAGS"
make -C SeqLib/fermi-lite -j"$JOBS" CFLAGS="$SUBMOD_CFLAGS"

echo ">>> [2/3] Configuring CMake (svaba + seqlib), ASan on, jemalloc OFF, asserts ON"
# CMAKE_BUILD_TYPE=Debug => -g and (crucially) NO -DNDEBUG, so asserts stay
# live and may catch the bug directly. Do NOT use RelWithDebInfo here (it
# defines NDEBUG and strips the asserts).
cmake -B build-asan \
  -DCMAKE_BUILD_TYPE=Debug \
  -DUSE_JEMALLOC=OFF \
  -DCMAKE_C_FLAGS="$CMAKE_SAN_FLAGS" \
  -DCMAKE_CXX_FLAGS="$CMAKE_SAN_FLAGS" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address"

echo ">>> [3/3] Building"
cmake --build build-asan -j"$JOBS"

echo
echo "================================================================"
echo "Built: ./build-asan/svaba"
echo
echo "Sanity-check ASan is actually wired in:"
echo "    nm build-asan/svaba | grep -c __asan_init     # expect > 0"
echo
echo "Run on the failing input (keep -p and data identical), loop until it trips:"
echo "    ASAN_OPTIONS=detect_leaks=0:abort_on_error=1:halt_on_error=1:fast_unwind_on_malloc=0 \\"
echo "      ./build-asan/svaba run <same args> 2> asan.log"
echo
echo "If stacks aren't symbolized, install a symbolizer (llvm / binutils):"
echo "    apt-get install -y llvm binutils    # then re-run"
echo
echo "RESTORE production submodule libs when done debugging:"
echo "    make -C SeqLib/bwa        clean && make -C SeqLib/bwa        -j$JOBS"
echo "    make -C SeqLib/fermi-lite clean && make -C SeqLib/fermi-lite -j$JOBS"
echo "================================================================"
