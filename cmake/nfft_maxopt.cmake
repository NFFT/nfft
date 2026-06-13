# CMake port of m4/ax_cc_maxopt.m4 (gcc/clang branch) plus the -march=native
# selection of m4/ax_gcc_archflag.m4.
#
# When the user does NOT supply their own CFLAGS, the Autotools build compiles the
# C library with aggressive optimization (-O3 -ffast-math -march=native ...). CMake's
# Release/RelWithDebInfo defaults omit -ffast-math and -march=native, leaving the
# *direct* transforms (tight FP loops in kernel/) ~10-25% slower (worst for float).
# We restore those so a plain `cmake` build matches a plain `./configure`.
#
# Opt-out parity (AX_CC_MAXOPT only runs when CFLAGS is unset): maxopt is skipped
# when CMAKE_C_FLAGS is non-empty, so user-supplied flags always win. The
# NFFT_ENABLE_MAXOPT option (default ON) is an additional explicit toggle.
#   plain `cmake ...`                          -> maxopt ON, host-tuned (-march=native)
#   `CFLAGS=... cmake` / -DCMAKE_C_FLAGS=...    -> maxopt OFF, user's flags used
#   `-DNFFT_ENABLE_MAXOPT=OFF`                  -> maxopt OFF, only CMAKE_BUILD_TYPE flags
# Each flag is probed, so unsupported ones (e.g. -march=native on Apple arm64) drop out.
#
# CAVEAT: literal -march=native is more aggressive than Autotools' AX_GCC_ARCHFLAG
# (which emits a named, often conservative -march=<arch>). Under the CodSpeed/Valgrind
# harness it can emit an instruction Valgrind cannot execute and makes instruction
# counts CPU-dependent, so the benchmark CI cells pass their own fixed, Valgrind-safe
# CMAKE_C_FLAGS; see .github/workflows/build-linux.yml.

include(CheckCCompilerFlag)

set(NFFT_MAXOPT_C_FLAGS "")

if(NOT CMAKE_C_FLAGS STREQUAL "")
  message(STATUS "NFFT maxopt skipped: user-supplied CMAKE_C_FLAGS take precedence "
                 "(AX_CC_MAXOPT-style opt-out): '${CMAKE_C_FLAGS}'")
elseif(NFFT_ENABLE_MAXOPT AND CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  foreach(_candidate -O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -march=native)
    string(MAKE_C_IDENTIFIER "NFFT_MAXOPT_HAS_${_candidate}" _ok)
    check_c_compiler_flag("${_candidate}" ${_ok})
    if(${_ok})
      list(APPEND NFFT_MAXOPT_C_FLAGS "${_candidate}")
    endif()
  endforeach()
  message(STATUS "NFFT maxopt C flags (AX_CC_MAXOPT parity): ${NFFT_MAXOPT_C_FLAGS}")
else()
  message(STATUS "NFFT maxopt disabled "
                 "(NFFT_ENABLE_MAXOPT=${NFFT_ENABLE_MAXOPT}, C compiler='${CMAKE_C_COMPILER_ID}'); "
                 "using CMAKE_BUILD_TYPE flags only")
endif()
