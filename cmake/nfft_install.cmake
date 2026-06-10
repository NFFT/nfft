# Install/export layer. Included from the root CMakeLists after the targets
# are defined.

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(NFFT_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/nfft3${NFFT_PREC_SUFFIX}"
    CACHE STRING "Install dir for NFFT3 CMake package files")

# Targets
# Only the serial library is in the export set, so find_package(NFFT3) yields just 
# NFFT3::nfft3 and imposes no OpenMP dependency on consumers. The _omp variant is 
# installed as a plain library file (usable via -lnfft3_omp / pkg-config) but is 
# not an exported CMake target.
install(TARGETS nfft3
  EXPORT NFFT3Targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
# No INCLUDES DESTINATION. The install include dir is already carried by the
# target's $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>.

if(NFFT_ENABLE_OPENMP)
  install(TARGETS nfft3_omp  # file only, no EXPORT
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

# ---- Public headers (precision-agnostic; installed once) ----------------
install(FILES
  ${PROJECT_SOURCE_DIR}/include/nfft3.h
  ${PROJECT_SOURCE_DIR}/include/nfft3mp.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# ---- pkg-config (reuse the Autotools template) --------------------------
set(prefix      "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}")
set(libdir      "${CMAKE_INSTALL_FULL_LIBDIR}")
set(VERSION     "${PROJECT_VERSION}")
set(PREC_SUFFIX "${NFFT_PREC_SUFFIX}")
configure_file(${PROJECT_SOURCE_DIR}/nfft3.pc.in
               ${PROJECT_BINARY_DIR}/nfft3${NFFT_PREC_SUFFIX}.pc @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/nfft3${NFFT_PREC_SUFFIX}.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

# ---- CMake package config ----------------------------------------------
configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/NFFT3Config.cmake.in
  ${PROJECT_BINARY_DIR}/NFFT3${NFFT_PREC_SUFFIX}Config.cmake
  INSTALL_DESTINATION ${NFFT_INSTALL_CMAKEDIR})

write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/NFFT3${NFFT_PREC_SUFFIX}ConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

install(FILES
  ${PROJECT_BINARY_DIR}/NFFT3${NFFT_PREC_SUFFIX}Config.cmake
  ${PROJECT_BINARY_DIR}/NFFT3${NFFT_PREC_SUFFIX}ConfigVersion.cmake
  ${PROJECT_SOURCE_DIR}/cmake/FindFFTW3.cmake          # shipped for find_dependency(FFTW3)
  DESTINATION ${NFFT_INSTALL_CMAKEDIR})

install(EXPORT NFFT3Targets
  FILE NFFT3Targets.cmake
  NAMESPACE NFFT3::
  DESTINATION ${NFFT_INSTALL_CMAKEDIR})
