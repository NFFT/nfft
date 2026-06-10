# Julia interface (double + shared only).
nfft_add_julia_module(nfft libnfftjulia.c)
if(HAVE_NFCT)
  nfft_add_julia_module(nfct libnfctjulia.c)
endif()
if(HAVE_NFST)
  nfft_add_julia_module(nfst libnfstjulia.c)
endif()
if(HAVE_NFSFT)
  nfft_add_julia_module(nfsft libnfsftjulia.c)
endif()
if(NFFT_ENABLE_APPLICATIONS)
  if(NFFT_ENABLE_OPENMP)
    nfft_add_julia_module(fastsum libfastsumjulia.c
      EXTRA_INC ${PROJECT_SOURCE_DIR}/applications/fastsum
      EXTRA_LINK fastsum_threads_lib fastsum_kernels)
  else()
    nfft_add_julia_module(fastsum libfastsumjulia.c
      EXTRA_INC ${PROJECT_SOURCE_DIR}/applications/fastsum
      EXTRA_LINK fastsum_lib fastsum_kernels)
  endif()
endif()
