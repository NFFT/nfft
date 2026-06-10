# Core kernel sources, always built.
set(NFFT_KERNEL_SOURCES
  util/assert.c util/bessel_i0.c util/bspline.c util/damp.c util/error.c
  util/float.c util/int.c util/lambda.c util/malloc.c util/print.c util/rand.c
  util/sinc.c util/sort.c util/thread.c util/time.c util/vector1.c
  util/vector2.c util/vector3.c util/version.c util/window.c
  nfft/nfft.c
  nfct/nfct.c
  nfst/nfst.c
  solver/solver.c)

# Add double-precision-only modules, maybe.
if(NOT NFFT_ENABLE_FLOAT AND NOT NFFT_ENABLE_LONG_DOUBLE)
  list(APPEND NFFT_KERNEL_SOURCES
    nnfft/nnfft.c
    nsfft/nsfft.c
    mri/mri.c
    fpt/fpt.c
    nfsft/legendre.c nfsft/nfsft.c
    nfsoft/nfsoft.c nfsoft/wigner.c)
endif()

# Apply the include dirs, install-interface, output/export name, and (Windows)
# DLL export macros shared by the nfft3 / nfft3_omp targets. `outname` is the
# bare library name (e.g. nfft3, nfft3f, nfft3_omp) — passed explicitly rather
# than computed, since the precision suffix is a configure-time constant.
function(nfft_configure_library tgt outname)
  set_target_properties(${tgt} PROPERTIES
    OUTPUT_NAME "${outname}"
    VERSION 4.0.3        # libtool -version-info 4:3:0 -> current.age.revision
    SOVERSION 0
    EXPORT_NAME "${outname}")
  target_include_directories(${tgt}
    PUBLIC
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
      ${NFFT_GENERATED_INCLUDE_DIR})   # build/: config.h + ticks.h
  if(UNIX)
    target_link_libraries(${tgt} PUBLIC m)  # Add math library.
  endif()
  if(WIN32 AND BUILD_SHARED_LIBS)
    # Library TUs export (COMPILING_NFFT) ; both library and consumers need NFFT_DLL.
    target_compile_definitions(${tgt} PUBLIC NFFT_DLL)
    target_compile_definitions(${tgt} PRIVATE COMPILING_NFFT)
  endif()
  if(MSVC)
    target_compile_options(${tgt} PRIVATE /bigobj)
  endif()
  # Aggressive optimization compiler flags. Empty when NFFT_ENABLE_MAXOPT is OFF.
  target_compile_options(${tgt} PRIVATE ${NFFT_MAXOPT_C_FLAGS})
endfunction()

add_library(nfft3 ${NFFT_KERNEL_SOURCES})
nfft_configure_library(nfft3 "nfft3${NFFT_PREC_SUFFIX}")
target_link_libraries(nfft3 PUBLIC FFTW3::fftw3)
add_library(NFFT3::nfft3${NFFT_PREC_SUFFIX} ALIAS nfft3)

# Add threaded library, maybe.
if(NFFT_ENABLE_OPENMP)
  add_library(nfft3_omp ${NFFT_KERNEL_SOURCES})
  nfft_configure_library(nfft3_omp "nfft3${NFFT_PREC_SUFFIX}_omp")
  target_link_libraries(nfft3_omp PUBLIC OpenMP::OpenMP_C FFTW3::fftw3 ${NFFT_OPENMP_ATOMIC_LIBS})
  if(FFTW3_OMP_FOUND)
    target_link_libraries(nfft3_omp PUBLIC FFTW3::fftw3_omp)
  elseif(FFTW3_THREADS_FOUND)
    target_link_libraries(nfft3_omp PUBLIC FFTW3::fftw3_threads)
  endif()
endif()

# Expose kernel source list to the interface subdirs. Interfaces recompile 
# these into a PIC static lib that links no FFTW of its own (FFTW flavour 
# chosen at the final shared-object link). Export ABSOLUTE paths.
set(_iface_abs "")
foreach(_s IN LISTS NFFT_KERNEL_SOURCES)
  list(APPEND _iface_abs "${CMAKE_CURRENT_SOURCE_DIR}/${_s}")
endforeach()
set(NFFT_KERNEL_SOURCES_ABS "${_iface_abs}" PARENT_SCOPE)
