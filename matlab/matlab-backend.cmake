# MATLAB mex backend. FindMatlab has no Octave support, hence a separate file 
# from the Octave helper. FindMatlab does NOT locate MATLAB's bundled libmwfftw3, 
# so we find it manually under ${Matlab_ROOT_DIR}/bin/<arch>.
set(Matlab_ROOT_DIR "${NFFT_WITH_MATLAB}")

# Request only MX_LIBRARY: matlab_add_mex compiles with CMake's own compiler, so
# MEX_COMPILER (the `mex` script) is not needed.
find_package(Matlab REQUIRED COMPONENTS MX_LIBRARY)

# Locate the bundled FFTW (libmwfftw3) — arch keyed off the mex extension.
set(_mwarch_mexa64    glnxa64)
set(_mwarch_mexmaca64 maca64)
set(_mwarch_mexmaci64 maci64)
set(_mwarch_mexw64    win64)
set(_arch "${_mwarch_${Matlab_MEX_EXTENSION}}")
if(NOT _arch)
  message(FATAL_ERROR "Unsupported Matlab_MEX_EXTENSION='${Matlab_MEX_EXTENSION}'; "
                      "extend the arch map in matlab/matlab-backend.cmake.")
endif()

# MATLAB ships the bundled FFTW as the VERSIONED libmwfftw3.so.3 (not an
# unversioned libmwfftw3.so), which find_library will not find AND which CMake
# decomposes to a broken `-lmwfftw3`. Try the unversioned name first (some
# installs have it); else glob the versioned lib and create a build-local
# `libmwfftw3.so` symlink to it, then find that. No fall-back to system fftw3.
set(_mlbin "${Matlab_ROOT_DIR}/bin/${_arch}")
find_library(Matlab_FFTW_LIBRARY NAMES mwfftw3 HINTS "${_mlbin}" NO_DEFAULT_PATH)
if(NOT Matlab_FFTW_LIBRARY)
  file(GLOB _mwfftw_versioned "${_mlbin}/libmwfftw3.so.*" "${_mlbin}/libmwfftw3.*.dylib")
  if(_mwfftw_versioned)
    list(GET _mwfftw_versioned 0 _mwfftw_real)
    file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/matlab-fftw")
    file(CREATE_LINK "${_mwfftw_real}"
         "${PROJECT_BINARY_DIR}/matlab-fftw/libmwfftw3.so" SYMBOLIC)
    find_library(Matlab_FFTW_LIBRARY NAMES mwfftw3
      HINTS "${PROJECT_BINARY_DIR}/matlab-fftw" NO_DEFAULT_PATH)
  endif()
endif()
if(NOT Matlab_FFTW_LIBRARY)
  message(FATAL_ERROR "MATLAB bundled FFTW (libmwfftw3[.so.N]) not found under ${_mlbin}.")
endif()
message(STATUS "MATLAB ${Matlab_VERSION} (ext=${Matlab_MEX_EXTENSION}); "
               "libmwfftw3=${Matlab_FFTW_LIBRARY}")

# This backend is provisional: The CMake wiring is verified against the stub fixture, 
# but it has not been built+run against real MATLAB.
message(WARNING "Building the MATLAB mex backend with CMake is experimental and not well tested.")

if(NFFT_ENABLE_MATLAB_THREADS)
  message(WARNING "Threaded MATLAB mex links libgomp, which can clash with MATLAB's "
                  "bundled OpenMP runtime (libiomp5/MKL) at load time. "
                  "verify on the MATLAB host or set NFFT_ENABLE_MATLAB_THREADS=OFF.")
endif()

# nfft_add_matlab_mex(<module> <source.c> [EXTRA_INC ...] [EXTRA_LINK ...])
# Builds build/matlab/<module>/<module>mex.<mexext> via matlab_add_mex, links the
# interface kernel + MATLAB's bundled FFTW + the matlab mex helper, stages .m
# wrappers (+ @f_hat) beside it, and installs lib -> libdir, .m -> datadir.
function(nfft_add_matlab_mex module source)
  cmake_parse_arguments(MX "" "" "EXTRA_INC;EXTRA_LINK" ${ARGN})
  set(tgt "matlab_${module}mex")
  set(outdir "${CMAKE_CURRENT_BINARY_DIR}/${module}")
  nfft_mex_helper_library(_helper matlab ${Matlab_INCLUDE_DIRS})

  if(NFFT_ENABLE_MATLAB_THREADS)
    set(_kernel nfft3_iface_omp)
  else()
    set(_kernel nfft3_iface)
  endif()

  matlab_add_mex(NAME ${tgt}
    SRC "${PROJECT_SOURCE_DIR}/matlab/${module}/${source}"
    OUTPUT_NAME ${module}mex
    LINK_TO ${MX_EXTRA_LINK} ${_kernel} ${_helper} ${Matlab_FFTW_LIBRARY}
    R2017b)                      # separate complex C API (code uses mxGetPr/mxGetPi)
  target_include_directories(${tgt} PRIVATE
    ${NFFT_GENERATED_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/matlab ${Matlab_INCLUDE_DIRS} ${MX_EXTRA_INC})
  set_target_properties(${tgt} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${outdir}")

  file(GLOB _m "${PROJECT_SOURCE_DIR}/matlab/${module}/*.m")
  add_custom_command(TARGET ${tgt} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory "${outdir}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_m} "${outdir}/"
    COMMENT "matlab/${module}: stage .m wrappers (matlab)")
  if(IS_DIRECTORY "${PROJECT_SOURCE_DIR}/matlab/${module}/@f_hat")
    add_custom_command(TARGET ${tgt} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_directory
              "${PROJECT_SOURCE_DIR}/matlab/${module}/@f_hat" "${outdir}/@f_hat")
  endif()

  set(_mdest "${CMAKE_INSTALL_DATADIR}/nfft/matlab/${module}")
  install(TARGETS ${tgt} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(FILES ${_m} DESTINATION "${_mdest}")
  if(IS_DIRECTORY "${PROJECT_SOURCE_DIR}/matlab/${module}/@f_hat")
    install(DIRECTORY "${PROJECT_SOURCE_DIR}/matlab/${module}/@f_hat" DESTINATION "${_mdest}")
  endif()
endfunction()
