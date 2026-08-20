# Examples. All programs are noinst. Sources are qualified by their subdir; 
# the helper's <subdir> arg routes each binary to build/examples/<subdir>/.
set(B ${CMAKE_CURRENT_BINARY_DIR})

# nfft (all precisions)
configure_file(nfft/simple_test.c.in ${B}/nfft/simple_test.c @ONLY)
nfft_add_serial_program(example_nfft_simple_test  nfft simple_test  ${B}/nfft/simple_test.c)
nfft_add_serial_program(example_nfft_ndft_fast    nfft ndft_fast    nfft/ndft_fast.c)
nfft_add_serial_program(example_nfft_taylor_nfft  nfft taylor_nfft  nfft/taylor_nfft.c)
nfft_add_serial_program(example_nfft_flags        nfft flags        nfft/flags.c)
nfft_add_serial_program(example_nfft_nfft_times   nfft nfft_times   nfft/nfft_times.c)
nfft_add_serial_program(example_nfft_nfast_native nfft nfast_native nfft/nfast_native.c)
target_compile_definitions(example_nfft_nfast_native PRIVATE
  NFAST_NATIVE_DATA="${CMAKE_CURRENT_SOURCE_DIR}/nfft/data/nfft_1d_8192_128.txt")

# nfft (OpenMP)
# nfft_benchomp_detail.c is built TWICE: serial (_single) and OpenMP (_threads),
# so the two can be compared — do not collapse them. The serial helpers link the
# serial library but still need the OpenMP RUNTIME lib (they call omp_* without
# -fopenmp): add ${OpenMP_C_LIBRARIES} only (no flags), mirroring $(OPENMP_LIBS).
if(NFFT_ENABLE_OPENMP)
  configure_file(nfft/simple_test_threads.c.in ${B}/nfft/simple_test_threads.c @ONLY)
  nfft_add_omp_program(example_nfft_simple_test_threads nfft simple_test_threads
    ${B}/nfft/simple_test_threads.c)
  nfft_add_omp_program(example_nfft_benchomp nfft nfft_benchomp nfft/nfft_benchomp.c)
  nfft_add_omp_program(example_nfft_benchomp_detail_threads nfft
    nfft_benchomp_detail_threads nfft/nfft_benchomp_detail.c)
  nfft_add_serial_program(example_nfft_benchomp_createdataset nfft
    nfft_benchomp_createdataset nfft/nfft_benchomp_createdataset.c)
  target_link_libraries(example_nfft_benchomp_createdataset PRIVATE ${OpenMP_C_LIBRARIES})
  nfft_add_serial_program(example_nfft_benchomp_detail_single nfft
    nfft_benchomp_detail_single nfft/nfft_benchomp_detail.c)
  target_link_libraries(example_nfft_benchomp_detail_single PRIVATE ${OpenMP_C_LIBRARIES})
endif()

# nfct/nfst (all precisions)
if(HAVE_NFCT)
  configure_file(nfct/simple_test.c.in ${B}/nfct/simple_test.c @ONLY)
  nfft_add_serial_program(example_nfct_simple_test nfct simple_test ${B}/nfct/simple_test.c)
endif()
if(HAVE_NFST)
  configure_file(nfst/simple_test.c.in ${B}/nfst/simple_test.c @ONLY)
  nfft_add_serial_program(example_nfst_simple_test nfst simple_test ${B}/nfst/simple_test.c)
endif()

# solver (all precisions)
configure_file(solver/simple_test.c.in ${B}/solver/simple_test.c @ONLY)
configure_file(solver/glacier.c.in     ${B}/solver/glacier.c     @ONLY)
nfft_add_serial_program(example_solver_simple_test solver simple_test ${B}/solver/simple_test.c)
nfft_add_serial_program(example_solver_glacier     solver glacier     ${B}/solver/glacier.c)

# nnfft (double-only)
if(HAVE_NNFFT)
  nfft_add_serial_program(example_nnfft_simple_test nnfft simple_test nnfft/simple_test.c)
  nfft_add_serial_program(example_nnfft_accuracy    nnfft accuracy    nnfft/accuracy.c)
endif()

# nsfft (double-only)
if(HAVE_NSFFT)
  nfft_add_serial_program(example_nsfft_simple_test nsfft simple_test nsfft/simple_test.c)
  nfft_add_serial_program(example_nsfft_nsfft_test  nsfft nsfft_test  nsfft/nsfft_test.c)
endif()

# fpt (double-only)
if(HAVE_FPT)
  nfft_add_serial_program(example_fpt_simple_test fpt simple_test fpt/simple_test.c)
endif()

# nfsoft (double-only)
if(HAVE_NFSOFT)
  nfft_add_serial_program(example_nfsoft_simple_test nfsoft simple_test nfsoft/simple_test.c)
endif()

# nfsft (double-only)
if(HAVE_NFSFT)
  # All plain .c here (no precision template). Serial always; OpenMP adds the
  # threads + benchomp set. The nfsft benchomp_createdataset/detail_single
  # are pure serial (unlike nfft's).
  nfft_add_serial_program(example_nfsft_simple_test nfsft simple_test nfsft/simple_test.c)
  if(NFFT_ENABLE_OPENMP)
    nfft_add_omp_program(example_nfsft_simple_test_threads nfsft simple_test_threads
      nfsft/simple_test_threads.c)
    nfft_add_omp_program(example_nfsft_benchomp nfsft nfsft_benchomp nfsft/nfsft_benchomp.c)
    nfft_add_omp_program(example_nfsft_benchomp_detail_threads nfsft
      nfsft_benchomp_detail_threads nfsft/nfsft_benchomp_detail.c)
    nfft_add_serial_program(example_nfsft_benchomp_createdataset nfsft
      nfsft_benchomp_createdataset nfsft/nfsft_benchomp_createdataset.c)
    nfft_add_serial_program(example_nfsft_benchomp_detail_single nfsft
      nfsft_benchomp_detail_single nfsft/nfsft_benchomp_detail.c)
  endif()
endif()
