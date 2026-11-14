# =====================================================================
# gen_git_version.cmake — stamp the build's git provenance into a header.
#
# Invoked from CMakeLists.txt both at configure time (to guarantee the
# header exists before the first compile) and as a build-time custom
# target (to keep the hash fresh across incremental local builds):
#
#   cmake -DSVABA_SRC_DIR=<repo root> \
#         -DSVABA_IN=<path to SvabaGitVersion.h.in> \
#         -DSVABA_OUT=<path to generated SvabaGitVersion.h> \
#         -P cmake/gen_git_version.cmake
#
# Writes via copy_if_different, so when HEAD/worktree haven't changed the
# output header is left byte-identical and no downstream recompile fires.
# =====================================================================

find_package(Git QUIET)

set(SVABA_GIT_HASH     "unknown")
set(SVABA_GIT_DESCRIBE "unknown")
set(SVABA_GIT_DIRTY    0)

if (Git_FOUND AND EXISTS "${SVABA_SRC_DIR}/.git")

  # full commit hash
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" -C "${SVABA_SRC_DIR}" rev-parse HEAD
    OUTPUT_VARIABLE  _hash
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    RESULT_VARIABLE  _rc_hash)
  if (_rc_hash EQUAL 0 AND _hash)
    set(SVABA_GIT_HASH "${_hash}")
  endif()

  # human-friendly describe (nearest tag + offset, or short hash fallback)
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" -C "${SVABA_SRC_DIR}" describe --tags --always --dirty
    OUTPUT_VARIABLE  _desc
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    RESULT_VARIABLE  _rc_desc)
  if (_rc_desc EQUAL 0 AND _desc)
    set(SVABA_GIT_DESCRIBE "${_desc}")
  endif()

  # uncommitted tracked changes => mark the build dirty
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" -C "${SVABA_SRC_DIR}" status --porcelain --untracked-files=no
    OUTPUT_VARIABLE  _dirty
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
  if (_dirty)
    set(SVABA_GIT_DIRTY 1)
  endif()

endif()

configure_file("${SVABA_IN}" "${SVABA_OUT}.tmp" @ONLY)
execute_process(COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                "${SVABA_OUT}.tmp" "${SVABA_OUT}")
file(REMOVE "${SVABA_OUT}.tmp")
