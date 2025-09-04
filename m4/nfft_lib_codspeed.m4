# Copyright (c) 2025 Jens Keiner, Stefan Kunis, Daniel Potts
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51
# Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# @synopsis NFFT_LIB_CODSPEED
# @summary Check configure options and assign variables related to the CodSpeed integration library.
# @category C++
#
# @version 2025-09-02
# @license GPLWithACException
# @author Jens Keiner <jens@nfft.org>
#
# If we find the CodSpeed integration library, set the shell variable `nfft_cv_codspeed' to `yes'
# and also set `nfft_codspeed_CPPFLAGS`, `nfft_codspeed_CXXFLAGS`, `nfft_codspeed_LDFLAGS`, and `nfft_codspeed_LIBS`
# to the respective flags nneded to compile and link with the library.
# Otherwise, set `ax_have_codspeed' to `no'.

AC_DEFUN([NFFT_CXX_COMPILE_STDCXX_17], [
    AX_CXX_COMPILE_STDCXX([17], [noext], [mandatory])
])

AC_DEFUN([NFFT_LIB_CODSPEED],
[
  AC_LANG_PUSH([C++])

  AC_REQUIRE([AC_PROG_CXX])
  AC_REQUIRE([AX_CXX_COMPILE_STDCXX_17])

  AC_ARG_WITH(codspeed, [AS_HELP_STRING([--with-codspeed=DIR],
  [use CodSpeed C++ library project (https://github.com/CodSpeedHQ/codspeed-cpp) in DIR])], with_codspeed=$withval, with_codspeed="no")
  
  if test "x$with_codspeed" != "xno"; then
    nfft_codspeed_CPPFLAGS="-I${with_codspeed}/google_benchmark/include -I${with_codspeed}/core/include -I${with_codspeed}/core/instrument-hooks/includes"
    nfft_codspeed_LDFLAGS="-L${with_codspeed}/google_benchmark/src -L${with_codspeed}/google_benchmark/codspeed"
  else
    nfft_codspeed_CPPFLAGS=""
    nfft_codspeed_LDFLAGS=""
  fi

  nfft_codspeed_CPPFLAGS="-DBENCHMARK_STATIC_DEFINE -DCODSPEED_ENABLED -DCODSPEED_INSTRUMENTATION $nfft_codspeed_CPPFLAGS"
  nfft_codspeed_CXXFLAGS=""
  nfft_codspeed_LDFLAGS="$nfft_codspeed_LDFLAGS"
  nfft_codspeed_LIBS="-lbenchmark -lcodspeed -linstrument_hooks"
  
  saved_CPPFLAGS="$CPPFLAGS"
  saved_CXXFLAGS="$CXXFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"
  
  CPPFLAGS="$nfft_codspeed_CPPFLAGS $CPPFLAGS"
  CXXFLAGS="$nfft_codspeed_CXXFLAGS $CXXFLAGS"
  LDFLAGS="$nfft_codspeed_LDFLAGS $LDFLAGS"
  LIBS="$nfft_codspeed_LIBS $LIBS"

  AC_CHECK_HEADER([benchmark/benchmark.h], [nfft_codspeed_header=yes], [nfft_codspeed_header=no])

  AC_MSG_CHECKING([whether CodSpeed C++ integration library test code compiles])
  AC_COMPILE_IFELSE([
    AC_LANG_SOURCE([[
      #include <benchmark/benchmark.h>
      #include <cstring>
      
      static void BM_StringCreation(benchmark::State& state) {
        for (auto _ : state)
          std::string empty_string;
      }
      BENCHMARK(BM_StringCreation);
      
      BENCHMARK_MAIN();
    ]])
  ], [
    AC_MSG_RESULT([yes])
    nfft_codspeed_compiles=yes
  ], [
    AC_MSG_RESULT([no])
    nfft_codspeed_compiles=no
  ])

  AC_MSG_CHECKING([whether CodSpeed C++ integration library test code links])
  AC_LINK_IFELSE([
    AC_LANG_SOURCE([[
      #include <benchmark/benchmark.h>
      #include <cstring>
      
      static void BM_StringCreation(benchmark::State& state) {
        for (auto _ : state)
          std::string empty_string;
      }
      BENCHMARK(BM_StringCreation);
      
      BENCHMARK_MAIN();
    ]])
  ], [
    AC_MSG_RESULT([yes])
    nfft_codspeed_links=yes
  ], [
    AC_MSG_RESULT([no])
    nfft_codspeed_links=no
  ])

  # Restore saved flags.
  CPPFLAGS="$saved_CPPFLAGS"
  CXXFLAGS="$saved_CXXFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  if test "x$nfft_codspeed_header" == "xyes" -a "x$nfft_codspeed_compiles" == "xyes" -a "x$nfft_codspeed_links" == "xyes"; then
    nfft_cv_codspeed=yes
  else
    nfft_cv_codspeed=no
    nfft_codspeed_CPPFLAGS=""
    nfft_codspeed_CXXFLAGS=""
    nfft_codspeed_LDFLAGS=""
    nfft_codspeed_LIBS=""
  fi

  AC_MSG_CHECKING([for CodSpeed C++ integration library])
  AC_MSG_RESULT([$nfft_cv_codspeed])

  AC_LANG_POP([C++])
])
