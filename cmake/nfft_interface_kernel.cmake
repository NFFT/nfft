# Interface-kernel static libraries.
# nfft3_iface     : kernel/*.c, compiled PIC, serial, NO FFTW linked.
# nfft3_iface_omp : same, compiled with OpenMP (only when NFFT_ENABLE_OPENMP).
# The consuming shared object (a Julia module or a mex file) links the FFTW it
# needs. One pair serves Julia, Octave, and MATLAB. Invoke from root scope after
# add_subdirectory(kernel) and only when an interface is enabled.
function(nfft_define_interface_kernels)
  if(TARGET nfft3_iface)
    return()
  endif()
  add_library(nfft3_iface STATIC ${NFFT_KERNEL_SOURCES_ABS})
  set_target_properties(nfft3_iface PROPERTIES POSITION_INDEPENDENT_CODE ON)
  # FFTW's header dir only. nfft3.h pulls in <fftw3.h>; on Linux it sits on the default 
  # include path, but on macOS (Homebrew /opt/homebrew/include) the kernel sources fail 
  # to compile without this, since CMake does not consume CPPFLAGS.
  target_include_directories(nfft3_iface
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    PRIVATE ${NFFT_GENERATED_INCLUDE_DIR} ${FFTW3_INCLUDE_DIR})
  target_compile_options(nfft3_iface PRIVATE ${NFFT_MAXOPT_C_FLAGS})
  if(UNIX)
    target_link_libraries(nfft3_iface PUBLIC m)
  endif()

  if(NFFT_ENABLE_OPENMP)
    add_library(nfft3_iface_omp STATIC ${NFFT_KERNEL_SOURCES_ABS})
    set_target_properties(nfft3_iface_omp PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_include_directories(nfft3_iface_omp
      PUBLIC ${PROJECT_SOURCE_DIR}/include
      PRIVATE ${NFFT_GENERATED_INCLUDE_DIR} ${FFTW3_INCLUDE_DIR})
    target_compile_options(nfft3_iface_omp PRIVATE ${NFFT_MAXOPT_C_FLAGS})
    target_link_libraries(nfft3_iface_omp
      PUBLIC OpenMP::OpenMP_C ${NFFT_OPENMP_ATOMIC_LIBS})
    if(UNIX)
      target_link_libraries(nfft3_iface_omp PUBLIC m)
    endif()
  endif()
endfunction()
