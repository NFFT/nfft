# Applications. All programs are noinst. Sources are qualified by their subdir; 
# the helper's <subdir> arg routes each binary to build/applications/<subdir>/.
# applications/iterS2 is currently omitted; the source is kept for later re-integration.
set(B ${CMAKE_CURRENT_BINARY_DIR})

# fastgauss / polarFFT / radon (all precisions)
configure_file(fastgauss/fastgauss.c.in ${B}/fastgauss/fastgauss.c @ONLY)
nfft_add_serial_program(app_fastgauss fastgauss fastgauss ${B}/fastgauss/fastgauss.c)

configure_file(polarFFT/polar_fft_test.c.in    ${B}/polarFFT/polar_fft_test.c    @ONLY)
configure_file(polarFFT/mpolar_fft_test.c.in   ${B}/polarFFT/mpolar_fft_test.c   @ONLY)
configure_file(polarFFT/linogram_fft_test.c.in ${B}/polarFFT/linogram_fft_test.c @ONLY)
nfft_add_serial_program(app_polar_fft_test    polarFFT polar_fft_test    ${B}/polarFFT/polar_fft_test.c)
nfft_add_serial_program(app_mpolar_fft_test   polarFFT mpolar_fft_test   ${B}/polarFFT/mpolar_fft_test.c)
nfft_add_serial_program(app_linogram_fft_test polarFFT linogram_fft_test ${B}/polarFFT/linogram_fft_test.c)

configure_file(radon/radon.c.in         ${B}/radon/radon.c         @ONLY)
configure_file(radon/inverse_radon.c.in ${B}/radon/inverse_radon.c @ONLY)
nfft_add_serial_program(app_radon         radon radon         ${B}/radon/radon.c)
nfft_add_serial_program(app_inverse_radon radon inverse_radon ${B}/radon/inverse_radon.c)

# fastsum (all precisions)
# libfastsum.la (serial), libkernels.la, and under OpenMP libfastsum_threads.la 
# (fastsum.c recompiled with OpenMP) as static libs
# fastsum.c/fastsum_test.c include local fastsum.h/kernels.h, so add the fastsum 
# source dir to the include path.
set(_fastsum_src ${PROJECT_SOURCE_DIR}/applications/fastsum)
set(_fastsum_inc
  ${NFFT_GENERATED_INCLUDE_DIR}        # build/ : config.h, ticks.h (must win)
  ${PROJECT_SOURCE_DIR}/include
  ${_fastsum_src})                     # fastsum.h, kernels.h

# PIC on the fastsum archives: besides the app executables, they are linked whole
# into the shared Julia fastsum module (julia/CMakeLists.txt). A shared object needs
# -fPIC objects (the default -fPIE is not enough). Same treatment as the mex
# helper (cmake/nfft_interfaces.cmake).
add_library(fastsum_lib STATIC ${_fastsum_src}/fastsum.c)
set_target_properties(fastsum_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(fastsum_lib PRIVATE ${_fastsum_inc})
target_link_libraries(fastsum_lib PRIVATE FFTW3::fftw3)

add_library(fastsum_kernels STATIC ${_fastsum_src}/kernels.c)
set_target_properties(fastsum_kernels PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(fastsum_kernels PRIVATE ${_fastsum_inc})
# kernels.c -> kernels.h -> nfft3.h -> <fftw3.h>, so it needs FFTW's include dir.
# Link FFTW3::fftw3 (as fastsum_lib does) to pick it up portably; ${_fastsum_inc}
# alone omits it, which only breaks where FFTW is non-system (macOS/Homebrew at
# /opt/homebrew/include) — Linux finds <fftw3.h> on the default system path.
target_link_libraries(fastsum_kernels PRIVATE FFTW3::fftw3)
# kernels.c calls libm (exp, …); as a STATIC archive its libm reference must be
# resolved after the archive on the link line. Propagate -lm as a usage
# requirement so consumers (and the helper's earlier -lm) link it in the right
# order — Autotools gets this transitively via @fftw3_LIBS@/libtool reordering.
if(UNIX)
  target_link_libraries(fastsum_kernels INTERFACE m)
endif()

# fastsum_test (serial)
nfft_add_serial_program(app_fastsum_test fastsum fastsum_test ${_fastsum_src}/fastsum_test.c)
target_include_directories(app_fastsum_test PRIVATE ${_fastsum_src})
target_link_libraries(app_fastsum_test PRIVATE fastsum_lib fastsum_kernels)

if(NFFT_ENABLE_OPENMP)
  add_library(fastsum_threads_lib STATIC ${_fastsum_src}/fastsum.c)
  set_target_properties(fastsum_threads_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_include_directories(fastsum_threads_lib PRIVATE ${_fastsum_inc})
  target_link_libraries(fastsum_threads_lib PRIVATE OpenMP::OpenMP_C FFTW3::fftw3)
  if(FFTW3_OMP_FOUND)
    target_link_libraries(fastsum_threads_lib PRIVATE FFTW3::fftw3_omp)
  elseif(FFTW3_THREADS_FOUND)
    target_link_libraries(fastsum_threads_lib PRIVATE FFTW3::fftw3_threads)
  endif()

  # fastsum_matlab (OpenMP stack when OpenMP on).
  nfft_add_omp_program(app_fastsum_matlab fastsum fastsum_matlab ${_fastsum_src}/fastsum_matlab.c)
  target_include_directories(app_fastsum_matlab PRIVATE ${_fastsum_src})
  target_link_libraries(app_fastsum_matlab PRIVATE fastsum_threads_lib fastsum_kernels)

  # fastsum_test_threads (OpenMP stack; same source as fastsum_test).
  nfft_add_omp_program(app_fastsum_test_threads fastsum fastsum_test_threads ${_fastsum_src}/fastsum_test.c)
  target_include_directories(app_fastsum_test_threads PRIVATE ${_fastsum_src})
  target_link_libraries(app_fastsum_test_threads PRIVATE fastsum_threads_lib fastsum_kernels)

  # benchomp: full-OpenMP (benchomp, detail_threads) link the threads stack.
  nfft_add_omp_program(app_fastsum_benchomp fastsum fastsum_benchomp ${_fastsum_src}/fastsum_benchomp.c)
  target_include_directories(app_fastsum_benchomp PRIVATE ${_fastsum_src})
  target_link_libraries(app_fastsum_benchomp PRIVATE fastsum_threads_lib fastsum_kernels)

  nfft_add_omp_program(app_fastsum_benchomp_detail_threads fastsum
    fastsum_benchomp_detail_threads ${_fastsum_src}/fastsum_benchomp_detail.c)
  target_include_directories(app_fastsum_benchomp_detail_threads PRIVATE ${_fastsum_src})
  target_link_libraries(app_fastsum_benchomp_detail_threads PRIVATE fastsum_threads_lib fastsum_kernels)

  # Serial benchomp helpers (serial library).
  nfft_add_serial_program(app_fastsum_benchomp_createdataset fastsum
    fastsum_benchomp_createdataset ${_fastsum_src}/fastsum_benchomp_createdataset.c)
  target_include_directories(app_fastsum_benchomp_createdataset PRIVATE ${_fastsum_src})
  target_link_libraries(app_fastsum_benchomp_createdataset PRIVATE fastsum_lib fastsum_kernels)

  nfft_add_serial_program(app_fastsum_benchomp_detail_single fastsum
    fastsum_benchomp_detail_single ${_fastsum_src}/fastsum_benchomp_detail.c)
  target_include_directories(app_fastsum_benchomp_detail_single PRIVATE ${_fastsum_src})
  target_link_libraries(app_fastsum_benchomp_detail_single PRIVATE fastsum_lib fastsum_kernels)
else()
  # fastsum_matlab (serial stack when OpenMP off).
  nfft_add_serial_program(app_fastsum_matlab fastsum fastsum_matlab ${_fastsum_src}/fastsum_matlab.c)
  target_include_directories(app_fastsum_matlab PRIVATE ${_fastsum_src})
  target_link_libraries(app_fastsum_matlab PRIVATE fastsum_lib fastsum_kernels)
endif()

# Double-only applications
if(HAVE_NFSFT)
  nfft_add_serial_program(app_fastsumS2    fastsumS2    fastsumS2    fastsumS2/fastsumS2.c)
  nfft_add_serial_program(app_quadratureS2 quadratureS2 quadratureS2 quadratureS2/quadratureS2.c)
endif()

if(HAVE_MRI)
  # mri2d
  nfft_add_serial_program(app_mri2d_construct_data_2d        mri/mri2d construct_data_2d        mri/mri2d/construct_data_2d.c)
  nfft_add_serial_program(app_mri2d_construct_data_inh_2d1d  mri/mri2d construct_data_inh_2d1d  mri/mri2d/construct_data_inh_2d1d.c)
  nfft_add_serial_program(app_mri2d_construct_data_inh_3d    mri/mri2d construct_data_inh_3d    mri/mri2d/construct_data_inh_3d.c)
  nfft_add_serial_program(app_mri2d_reconstruct_data_gridding mri/mri2d reconstruct_data_gridding mri/mri2d/reconstruct_data_gridding.c)
  nfft_add_serial_program(app_mri2d_reconstruct_data_2d      mri/mri2d reconstruct_data_2d      mri/mri2d/reconstruct_data_2d.c)
  nfft_add_serial_program(app_mri2d_reconstruct_data_inh_3d  mri/mri2d reconstruct_data_inh_3d  mri/mri2d/reconstruct_data_inh_3d.c)
  nfft_add_serial_program(app_mri2d_reconstruct_data_inh_2d1d mri/mri2d reconstruct_data_inh_2d1d mri/mri2d/reconstruct_data_inh_2d1d.c)
  nfft_add_serial_program(app_mri2d_reconstruct_data_inh_nnfft mri/mri2d reconstruct_data_inh_nnfft mri/mri2d/reconstruct_data_inh_nnfft.c)

  # mri3d
  nfft_add_serial_program(app_mri3d_construct_data_2d1d      mri/mri3d construct_data_2d1d      mri/mri3d/construct_data_2d1d.c)
  nfft_add_serial_program(app_mri3d_construct_data_3d        mri/mri3d construct_data_3d        mri/mri3d/construct_data_3d.c)
  nfft_add_serial_program(app_mri3d_reconstruct_data_gridding mri/mri3d reconstruct_data_gridding mri/mri3d/reconstruct_data_gridding.c)
  nfft_add_serial_program(app_mri3d_reconstruct_data_2d1d    mri/mri3d reconstruct_data_2d1d    mri/mri3d/reconstruct_data_2d1d.c)
  nfft_add_serial_program(app_mri3d_reconstruct_data_3d      mri/mri3d reconstruct_data_3d      mri/mri3d/reconstruct_data_3d.c)
endif()
