# Shared helpers for the language-interface subdirs.
# nfft_add_julia_module(<module> <source.c> [EXTRA_INC <dir>...] [EXTRA_LINK <tgt>...])
# Builds build/julia/<module>/lib<module>julia.<ext> (SHARED, PREFIX "" + the "lib"
# in OUTPUT_NAME so the name is Windows-correct too), no version/soname, links the
# interface kernel (omp under OpenMP) + system FFTW, stages the .jl wrappers beside
# it, and CO-INSTALLS lib + .jl into datadir/nfft/julia/<module>/.
function(nfft_add_julia_module module source)
  cmake_parse_arguments(JL "" "" "EXTRA_INC;EXTRA_LINK" ${ARGN})
  set(tgt "julia_${module}")
  set(outdir "${CMAKE_CURRENT_BINARY_DIR}/${module}")

  add_library(${tgt} SHARED "${CMAKE_CURRENT_SOURCE_DIR}/${module}/${source}")
  set_target_properties(${tgt} PROPERTIES
    OUTPUT_NAME "lib${module}julia" PREFIX ""
    LIBRARY_OUTPUT_DIRECTORY "${outdir}" RUNTIME_OUTPUT_DIRECTORY "${outdir}"
    NO_SONAME ON VERSION "" SOVERSION "")
  target_include_directories(${tgt}
    PRIVATE ${NFFT_GENERATED_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/include ${JL_EXTRA_INC})
  # The .jl wrappers resolve library symbols by name at runtime (ccall/dlsym), not
  # at link time. The interface kernel is a static archive, so a normal link prunes
  # any object no module source references at link time. Link the archive whole so 
  # every kernel symbol survives and is exported — matching libtool's convenience-lib
  # behavior in the Autotools build.
  # (LINK_LIBRARY:WHOLE_ARCHIVE needs CMake >= 3.24)
  if(NFFT_ENABLE_OPENMP)
    set(_iface nfft3_iface_omp)
  else()
    set(_iface nfft3_iface)
  endif()
  target_link_libraries(${tgt} PRIVATE
    ${JL_EXTRA_LINK} "$<LINK_LIBRARY:WHOLE_ARCHIVE,${_iface}>")
  if(NFFT_ENABLE_OPENMP)
    if(FFTW3_OMP_FOUND)
      target_link_libraries(${tgt} PRIVATE FFTW3::fftw3_omp)
    elseif(FFTW3_THREADS_FOUND)
      target_link_libraries(${tgt} PRIVATE FFTW3::fftw3_threads)
    endif()
  endif()
  target_link_libraries(${tgt} PRIVATE FFTW3::fftw3)

  file(GLOB _jl "${CMAKE_CURRENT_SOURCE_DIR}/${module}/*.jl")
  add_custom_command(TARGET ${tgt} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_jl} "${outdir}/"
    COMMENT "julia/${module}: stage .jl wrappers")

  # Co-install lib + .jl so @__DIR__ resolves post-install.
  set(_jdest "${CMAKE_INSTALL_DATADIR}/nfft/julia/${module}")
  install(TARGETS ${tgt} LIBRARY DESTINATION "${_jdest}" RUNTIME DESTINATION "${_jdest}")
  install(FILES ${_jl} DESTINATION "${_jdest}")
endfunction()
# nfft_mex_helper_library(<outvar> <backend> <inc_dirs>...)
# Shared mex helper (matlab/args.c + malloc.c) as a PIC static lib, once per
# backend ("octave" or "matlab").
function(nfft_mex_helper_library outvar backend)
  set(tgt "nfft_mex_${backend}")
  if(NOT TARGET ${tgt})
    add_library(${tgt} STATIC
      ${PROJECT_SOURCE_DIR}/matlab/args.c ${PROJECT_SOURCE_DIR}/matlab/malloc.c)
    set_target_properties(${tgt} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_include_directories(${tgt} PRIVATE
      ${NFFT_GENERATED_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/include
      ${PROJECT_SOURCE_DIR}/matlab ${ARGN})
    if(NFFT_ENABLE_MATLAB_THREADS)
      target_link_libraries(${tgt} PRIVATE OpenMP::OpenMP_C)
    endif()
  endif()
  set(${outvar} ${tgt} PARENT_SCOPE)
endfunction()

# nfft_add_octave_mex(<module> <source.c> [EXTRA_INC ...] [EXTRA_LINK ...])
# Builds build/matlab/<module>/<module>mex.mex (MODULE, PREFIX "", SUFFIX .mex),
# links interface kernel + system FFTW + Octave libs + the octave mex helper,
# stages the .m wrappers (+ @f_hat) beside it, and installs lib + .m to datadir.
function(nfft_add_octave_mex module source)
  cmake_parse_arguments(MX "" "" "EXTRA_INC;EXTRA_LINK" ${ARGN})
  set(tgt "octave_${module}mex")
  set(outdir "${CMAKE_CURRENT_BINARY_DIR}/${module}")
  nfft_mex_helper_library(_helper octave ${Octave_INCLUDE_DIRS})

  add_library(${tgt} MODULE "${PROJECT_SOURCE_DIR}/matlab/${module}/${source}")
  set_target_properties(${tgt} PROPERTIES
    OUTPUT_NAME "${module}mex" PREFIX "" SUFFIX ".${Octave_MEX_EXTENSION}"
    LIBRARY_OUTPUT_DIRECTORY "${outdir}")
  target_include_directories(${tgt} PRIVATE
    ${NFFT_GENERATED_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/matlab ${Octave_INCLUDE_DIRS} ${MX_EXTRA_INC})
  # Do NOT define HAVE_OCTAVE — <mex.h> defines it (octave/mex.h:56).

  if(NFFT_ENABLE_MATLAB_THREADS)
    target_link_libraries(${tgt} PRIVATE ${MX_EXTRA_LINK} nfft3_iface_omp)
    if(FFTW3_OMP_FOUND)
      target_link_libraries(${tgt} PRIVATE FFTW3::fftw3_omp)
    elseif(FFTW3_THREADS_FOUND)
      target_link_libraries(${tgt} PRIVATE FFTW3::fftw3_threads)
    endif()
  else()
    target_link_libraries(${tgt} PRIVATE ${MX_EXTRA_LINK} nfft3_iface)
  endif()
  target_link_libraries(${tgt} PRIVATE FFTW3::fftw3 ${_helper})
  target_link_directories(${tgt} PRIVATE ${Octave_LIBRARY_DIRS})
  target_link_libraries(${tgt} PRIVATE ${Octave_LIBRARIES})

  file(GLOB _m "${PROJECT_SOURCE_DIR}/matlab/${module}/*.m")
  add_custom_command(TARGET ${tgt} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory "${outdir}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_m} "${outdir}/"
    COMMENT "matlab/${module}: stage .m wrappers (octave)")
  if(IS_DIRECTORY "${PROJECT_SOURCE_DIR}/matlab/${module}/@f_hat")
    add_custom_command(TARGET ${tgt} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_directory
              "${PROJECT_SOURCE_DIR}/matlab/${module}/@f_hat" "${outdir}/@f_hat")
  endif()

  # Install: lib -> libdir (final name <mod>mex.mex), .m + @f_hat -> datadir.
  set(_mdest "${CMAKE_INSTALL_DATADIR}/nfft/matlab/${module}")
  install(TARGETS ${tgt} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(FILES ${_m} DESTINATION "${_mdest}")
  if(IS_DIRECTORY "${PROJECT_SOURCE_DIR}/matlab/${module}/@f_hat")
    install(DIRECTORY "${PROJECT_SOURCE_DIR}/matlab/${module}/@f_hat" DESTINATION "${_mdest}")
  endif()
endfunction()
