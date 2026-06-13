# Helpers for the example/application programs. In Autotools these are
# noinst_PROGRAMS (built by `make`, never installed); here they are ordinary
# targets, part of all. Two link lines:
#   serial : links the precision-suffixed serial library  (nfft3)
#   omp    : links _omp + OpenMP + the FFTW omp/threads lib + atomics
#
# All example/application targets live in a single CMakeLists per area. Many 
# programs share the name "simple_test", so the <subdir> argument routes each
# binary to build/<area>/<subdir>/, avoiding collision. Global target names must still be unique; OUTPUT_NAME restores the original name.
#
# Usage:
#   nfft_add_serial_program(<target> <subdir> <output_name> <source>...)
#   nfft_add_omp_program(<target> <subdir> <output_name> <source>...)

function(nfft_add_serial_program target subdir output)
  add_executable(${target} ${ARGN})
  set_target_properties(${target} PROPERTIES
    OUTPUT_NAME ${output}
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${subdir})
  target_include_directories(${target} PRIVATE
    ${NFFT_GENERATED_INCLUDE_DIR}        # build/ : config.h, ticks.h (must win)
    ${PROJECT_SOURCE_DIR}/include)
  target_link_libraries(${target} PRIVATE nfft3 FFTW3::fftw3)
  if(UNIX)
    target_link_libraries(${target} PRIVATE m)
  endif()
endfunction()

function(nfft_add_omp_program target subdir output)
  add_executable(${target} ${ARGN})
  set_target_properties(${target} PROPERTIES
    OUTPUT_NAME ${output}
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${subdir})
  target_include_directories(${target} PRIVATE
    ${NFFT_GENERATED_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/include)
  target_link_libraries(${target} PRIVATE
    nfft3_omp OpenMP::OpenMP_C FFTW3::fftw3 ${NFFT_OPENMP_ATOMIC_LIBS})
  if(FFTW3_OMP_FOUND)
    target_link_libraries(${target} PRIVATE FFTW3::fftw3_omp)
  elseif(FFTW3_THREADS_FOUND)
    target_link_libraries(${target} PRIVATE FFTW3::fftw3_threads)
  endif()
  if(UNIX)
    target_link_libraries(${target} PRIVATE m)
  endif()
endfunction()
