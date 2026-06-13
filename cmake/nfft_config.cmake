include(CheckIncludeFile)
include(CheckFunctionExists)
include(CheckSymbolExists)
include(CheckTypeSize)

# Window function
set(KAISER_BESSEL OFF)
set(GAUSSIAN OFF)
set(B_SPLINE OFF)
set(SINC_POWER OFF)
set(DIRAC_DELTA OFF)
if(NFFT_WINDOW STREQUAL "kaiserbessel")
  set(KAISER_BESSEL ON)
  set(WINDOW_NAME "Kaiser Bessel")
elseif(NFFT_WINDOW STREQUAL "gaussian")
  set(GAUSSIAN ON)
  set(WINDOW_NAME "Gaussian")
elseif(NFFT_WINDOW STREQUAL "bspline")
  set(B_SPLINE ON)
  set(WINDOW_NAME "B Spline")
elseif(NFFT_WINDOW STREQUAL "sinc")
  set(SINC_POWER ON)
  set(WINDOW_NAME "Sinc Power")
elseif(NFFT_WINDOW STREQUAL "dirac")
  set(DIRAC_DELTA ON)
  set(WINDOW_NAME "Dirac Delta")
else()
  message(FATAL_ERROR "Unknown NFFT_WINDOW='${NFFT_WINDOW}'")
endif()

# Precision
if(NFFT_ENABLE_FLOAT)
  set(NFFT_SINGLE ON)
elseif(NFFT_ENABLE_LONG_DOUBLE)
  set(NFFT_LDOUBLE ON)
endif()

# Benchmarks prefix
if(NOT DEFINED BENCHMARKS_PREFIX)
  set(BENCHMARKS_PREFIX "")
endif()

# Header checks
foreach(h math stdio stdlib time sys/time complex string float limits stdarg
          stddef sys/types stdint inttypes stdbool malloc c_asm intrinsics
          mach/mach_time)
  string(TOUPPER "${h}" H)
  string(REPLACE "/" "_" H "${H}")
  string(REPLACE "." "_" H "${H}")
  check_include_file("${h}.h" HAVE_${H}_H)
endforeach()

# Type checks
check_type_size("long double" SIZEOF_LONG_DOUBLE LANGUAGE C)
if(HAVE_SIZEOF_LONG_DOUBLE)
  set(HAVE_LONG_DOUBLE ON)
endif()
set(CMAKE_EXTRA_INCLUDE_FILES "sys/time.h")
check_type_size("hrtime_t" SIZEOF_HRTIME_T LANGUAGE C)
unset(CMAKE_EXTRA_INCLUDE_FILES)
if(HAVE_SIZEOF_HRTIME_T)
  set(HAVE_HRTIME_T ON)
endif()

# Function checks
foreach(fn gethrtime read_real_time time_base_to_time clock_gettime
           mach_absolute_time memset posix_memalign memalign sysctl abort
           snprintf sqrt sleep usleep nanosleep drand48 srand48 gethostname)
  string(TOUPPER "${fn}" FN)
  check_function_exists("${fn}" HAVE_${FN})
endforeach()

# Declaration checks
# AC_CHECK_DECLS sets HAVE_DECL_X to 0 or 1 (always defined). We emit the same
# via #define lines accumulated into NFFT_HAVE_DECLS.
set(_nfft_decl_defs "")
set(CMAKE_REQUIRED_QUIET ON)

# Base (double) math family
set(_math_base
  copysign nextafter nan ceil floor nearbyint rint round lrint lround llrint
  llround trunc fmod remainder remquo fdim fmax fmin fma fabs sqrt cbrt hypot
  exp exp2 expm1 log log2 log10 log1p logb ilogb modf frexp ldexp scalbn scalbln
  pow cos sin tan cosh sinh tanh acos asin atan atan2 acosh asinh atanh tgamma
  lgamma j0 j1 jn y0 y1 yn erf erfc creal cimag cabs carg conj cproj csqrt cexp
  clog cpow csin ccos ctan casin cacos catan csinh ccosh ctanh casinh cacosh
  catanh)

# Build double/float/long-double variants
set(_decl_names memalign posix_memalign sleep nanosleep drand48 srand48)
foreach(m IN LISTS _math_base)
  list(APPEND _decl_names "${m}" "${m}f" "${m}l")
endforeach()

foreach(name IN LISTS _decl_names)
  string(TOUPPER "${name}" NAME)
  check_symbol_exists("${name}" "math.h;complex.h;stdlib.h;unistd.h;time.h" _decl_${name})
  if(_decl_${name})
    string(APPEND _nfft_decl_defs "#define HAVE_DECL_${NAME} 1\n")
  else()
    string(APPEND _nfft_decl_defs "#define HAVE_DECL_${NAME} 0\n")
  endif()
  unset(_decl_${name} CACHE)
endforeach()
unset(CMAKE_REQUIRED_QUIET)
set(NFFT_HAVE_DECLS "${_nfft_decl_defs}")

# Additional headers
foreach(h alloca dlfcn strings sys/stat unistd)
  string(TOUPPER "${h}" H)
  string(REPLACE "/" "_" H "${H}")
  string(REPLACE "." "_" H "${H}")
  check_include_file("${h}.h" HAVE_${H}_H)
endforeach()

# alloca() function + STDC_HEADERS (these are effectively always true on the
# POSIX/C99 toolchains NFFT targets; mirror the Autotools AC_FUNC_ALLOCA /
# STDC_HEADERS results).
# alloca is typically a compiler builtin/macro, so probe via a compile test
# (AC_FUNC_ALLOCA parity): available if <alloca.h> exists or the builtin works.
check_symbol_exists(alloca "alloca.h" _nfft_alloca_in_header)
if(_nfft_alloca_in_header)
  set(HAVE_ALLOCA 1)
else()
  include(CheckCSourceCompiles)
  check_c_source_compiles(
    "int main(void){ char *p = (char*)__builtin_alloca(8); return p != 0; }"
    _nfft_alloca_builtin)
  if(_nfft_alloca_builtin)
    set(HAVE_ALLOCA 1)
  endif()
endif()
check_function_exists(vprintf HAVE_VPRINTF)
if(HAVE_STDLIB_H AND HAVE_STDARG_H AND HAVE_STRING_H AND HAVE_FLOAT_H)
  set(STDC_HEADERS 1)
endif()

# libm
set(HAVE_LIBM 1)

# Module enable flags
set(HAVE_NFCT 1)
set(HAVE_NFST 1)
if(NOT NFFT_ENABLE_FLOAT AND NOT NFFT_ENABLE_LONG_DOUBLE)
  set(HAVE_NNFFT 1)
  set(HAVE_NSFFT 1)
  set(HAVE_MRI 1)
  set(HAVE_FPT 1)
  set(HAVE_NFSFT 1)
  set(HAVE_NFSOFT 1)
endif()

# FFTW threads default when OpenMP is enabled
if(NFFT_ENABLE_OPENMP)
  if(NOT DEFINED HAVE_FFTW_THREADS)
    set(HAVE_FFTW_THREADS 1)
  endif()
endif()

# MATLAB argument checks
set(MATLAB_ARGCHECKS 1)

# Exhaustive unit tests option
if(NFFT_ENABLE_EXHAUSTIVE_UNIT_TESTS)
  set(NFFT_EXHAUSTIVE_UNIT_TESTS 1)
endif()

# Sizes of standard types
check_type_size("double"             SIZEOF_DOUBLE             LANGUAGE C)
check_type_size("float"              SIZEOF_FLOAT              LANGUAGE C)
check_type_size("int"                SIZEOF_INT                LANGUAGE C)
check_type_size("long"               SIZEOF_LONG               LANGUAGE C)
check_type_size("long long"          SIZEOF_LONG_LONG          LANGUAGE C)
check_type_size("ptrdiff_t"          SIZEOF_PTRDIFF_T          LANGUAGE C)
check_type_size("size_t"             SIZEOF_SIZE_T             LANGUAGE C)
check_type_size("unsigned int"       SIZEOF_UNSIGNED_INT       LANGUAGE C)
check_type_size("unsigned long"      SIZEOF_UNSIGNED_LONG      LANGUAGE C)
check_type_size("unsigned long long" SIZEOF_UNSIGNED_LONG_LONG LANGUAGE C)

# uintptr_t presence (AC_TYPE_UINTPTR_T).
check_type_size("uintptr_t" SIZEOF_UINTPTR_T LANGUAGE C)
if(HAVE_SIZEOF_UINTPTR_T)
  set(HAVE_UINTPTR_T 1)
endif()

# Package identification (autoconf AC_INIT parity)
set(NFFT_PACKAGE         "nfft")
set(NFFT_PACKAGE_NAME    "NFFT")
set(NFFT_PACKAGE_TARNAME "nfft")
set(NFFT_PACKAGE_BUGREPORT "mail@nfft.org")
set(NFFT_PACKAGE_URL     "")
set(NFFT_PACKAGE_VERSION_STR "${PROJECT_VERSION}${NFFT_VERSION_TYPE}")
set(NFFT_PACKAGE_STRING  "NFFT ${PROJECT_VERSION}${NFFT_VERSION_TYPE}")
set(LT_OBJDIR ".libs/")

# 'restrict' keyword spelling (AC_C_RESTRICT). GCC/Clang accept __restrict__.
set(NFFT_RESTRICT "__restrict__")
