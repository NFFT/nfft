# CUnit test suite.

# Shared source list.
set(NFFT_TEST_SOURCES
  check.c
  util.c    util.h
  bench_emit.c bench_emit.h
  reflect.c reflect.h
  bspline.c bspline.h
  bessel.c  bessel.h
  nfft.c    nfft.h
  nfct.c    nfct.h
  nfst.c    nfst.h)

  # Non-threaded test program.
add_executable(checkall ${NFFT_TEST_SOURCES})

# Generated includes first (build/), then source includes (include/).
target_include_directories(checkall PRIVATE
  ${NFFT_GENERATED_INCLUDE_DIR}      # config.h, ticks.h
  ${PROJECT_SOURCE_DIR}/include)     # nfft3.h, infft.h
target_link_libraries(checkall PRIVATE nfft3 FFTW3::fftw3 CUnit::CUnit)

if(UNIX)
# Add libm.
  target_link_libraries(checkall PRIVATE m)
endif()

add_test(NAME checkall COMMAND checkall)

# Optional threaded test program.
if(NFFT_ENABLE_OPENMP)
  add_executable(checkall_threads ${NFFT_TEST_SOURCES})
  target_include_directories(checkall_threads PRIVATE
    ${NFFT_GENERATED_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/include)
  target_link_libraries(checkall_threads PRIVATE
    nfft3_omp OpenMP::OpenMP_C FFTW3::fftw3 ${NFFT_OPENMP_ATOMIC_LIBS} CUnit::CUnit)
  if(FFTW3_OMP_FOUND)
    target_link_libraries(checkall_threads PRIVATE FFTW3::fftw3_omp)
  elseif(FFTW3_THREADS_FOUND)
    target_link_libraries(checkall_threads PRIVATE FFTW3::fftw3_threads)
  endif()
  if(UNIX)
    target_link_libraries(checkall_threads PRIVATE m)
  endif()
  add_test(NAME checkall_threads COMMAND checkall_threads)
endif()
