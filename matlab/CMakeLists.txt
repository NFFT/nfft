# MATLAB/Octave mex interface.
if(NFFT_WITH_OCTAVE)
  find_package(Octave REQUIRED)
  set(_mex nfft_add_octave_mex)
elseif(NFFT_WITH_MATLAB)
  include("${CMAKE_CURRENT_SOURCE_DIR}/matlab-backend.cmake")   # defines nfft_add_matlab_mex
  set(_mex nfft_add_matlab_mex)
endif()

cmake_language(CALL ${_mex} nfft nfftmex.c)
cmake_language(CALL ${_mex} nfct nfctmex.c)
cmake_language(CALL ${_mex} nfst nfstmex.c)
if(NFFT_ENABLE_APPLICATIONS)
  if(NFFT_ENABLE_MATLAB_THREADS)
    cmake_language(CALL ${_mex} fastsum fastsummex.c
      EXTRA_INC ${PROJECT_SOURCE_DIR}/applications/fastsum
      EXTRA_LINK fastsum_threads_lib fastsum_kernels)
  else()
    cmake_language(CALL ${_mex} fastsum fastsummex.c
      EXTRA_INC ${PROJECT_SOURCE_DIR}/applications/fastsum
      EXTRA_LINK fastsum_lib fastsum_kernels)
  endif()
endif()
if(HAVE_NFSFT)
  cmake_language(CALL ${_mex} nfsft nfsftmex.c)
endif()
if(HAVE_NFSOFT)
  cmake_language(CALL ${_mex} nfsoft nfsoftmex.c)
endif()
if(HAVE_NNFFT)
  cmake_language(CALL ${_mex} nnfft nnfftmex.c)
endif()
if(HAVE_FPT)
  cmake_language(CALL ${_mex} fpt fptmex.c)
endif()

# ---- Data-only installs --------------------------------------------------
install(FILES
  ${PROJECT_SOURCE_DIR}/matlab/infft1d/coeffs_2048.mat
  ${PROJECT_SOURCE_DIR}/matlab/infft1d/infft.m
  ${PROJECT_SOURCE_DIR}/matlab/infft1d/README
  ${PROJECT_SOURCE_DIR}/matlab/infft1d/simple_test.m
  ${PROJECT_SOURCE_DIR}/matlab/infft1d/test_adjoint.m
  ${PROJECT_SOURCE_DIR}/matlab/infft1d/test_overdetermined.m
  ${PROJECT_SOURCE_DIR}/matlab/infft1d/test_underdetermined.m
  DESTINATION ${CMAKE_INSTALL_DATADIR}/nfft/matlab/infft1d)

# matlab/tests dist_matlabtests_DATA (delegates + runners; generated *UnitTests.m excluded)
file(GLOB _testdata
  "${PROJECT_SOURCE_DIR}/matlab/tests/*TestcaseDelegate*.m"
  "${PROJECT_SOURCE_DIR}/matlab/tests/*TestcaseCheckDelegate.m"
  "${PROJECT_SOURCE_DIR}/matlab/tests/*TestcaseInitDelegate.m"
  "${PROJECT_SOURCE_DIR}/matlab/tests/*TestcaseTrafoDelegate.m"
  "${PROJECT_SOURCE_DIR}/matlab/tests/*UnitTestsRunAndExit.m"
  "${PROJECT_SOURCE_DIR}/matlab/tests/nfsftUnitTests.m"
  "${PROJECT_SOURCE_DIR}/matlab/tests/nfsoftUnitTests.m")
install(FILES ${_testdata} DESTINATION ${CMAKE_INSTALL_DATADIR}/nfft/matlab/tests)

# Interface tests (only nfftUnitTests.m / nfsftUnitTests.m are .m.in) ----
if(NFFT_WITH_OCTAVE AND NFFT_ENABLE_TESTS)
  if(NFFT_ENABLE_EXHAUSTIVE_UNIT_TESTS)
    set(matlab_exhaustive_unit_tests_flag 1)
  else()
    set(matlab_exhaustive_unit_tests_flag 0)
  endif()
  set(ABS_SRCDIR "${PROJECT_SOURCE_DIR}")
  set(_tbin "${CMAKE_CURRENT_BINARY_DIR}/tests")
  file(MAKE_DIRECTORY "${_tbin}")
  configure_file(${PROJECT_SOURCE_DIR}/matlab/tests/nfftUnitTests.m.in  "${_tbin}/nfftUnitTests.m" @ONLY)
  configure_file(${PROJECT_SOURCE_DIR}/matlab/tests/nfsftUnitTests.m.in "${_tbin}/nfsftUnitTests.m" @ONLY)
  file(GLOB _testm "${PROJECT_SOURCE_DIR}/matlab/tests/*.m")
  file(COPY ${_testm} DESTINATION "${_tbin}")        # delegates + runners (source .m)

  # The runner needs each listed module dir on the Octave path. The nfsoft
  # tests additionally need nfft/nfsft's lgwt.m (matlab/nfsft/lgwt.m), mirroring
  # check_nfsoft_octave.sh.in which addpath()s matlab/{tests,nfsft,nfsoft}.
  function(nfft_add_octave_interface_test name runner)
    set(_pathargs --path "${_tbin}")
    foreach(_moddir IN LISTS ARGN)
      list(APPEND _pathargs --path "${CMAKE_CURRENT_BINARY_DIR}/${_moddir}")
    endforeach()
    add_test(NAME ${name}
      COMMAND ${Octave_CLI} --norc --no-window-system
              ${_pathargs}
              --eval "${runner}")
    set_tests_properties(${name} PROPERTIES LABELS "interface;octave" TIMEOUT 1800)
  endfunction()

  nfft_add_octave_interface_test(octave_nfft_unittests nfftUnitTestsRunAndExit nfft)
  if(HAVE_NFSFT)
    nfft_add_octave_interface_test(octave_nfsft_unittests nfsftUnitTestsRunAndExit nfsft)
  endif()
  if(HAVE_NFSOFT)
    nfft_add_octave_interface_test(octave_nfsoft_unittests nfsoftUnitTestsRunAndExit nfsft nfsoft)
  endif()
endif()
