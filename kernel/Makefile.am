if HAVE_NFCT
  LIB_NFCT=nfct/libnfct.la
  DIR_NFCT=nfct
else
  LIB_NFCT=
endif

if HAVE_NFST
  LIB_NFST=nfst/libnfst.la
  DIR_NFST=nfst
else
  LIB_NFST=
  DIR_NFST=
endif

if HAVE_NFSFT
  LIB_NFSFT=nfsft/libnfsft.la
  DIR_NFSFT=nfsft
if HAVE_THREADS
  LIB_NFSFT_THREADS=nfsft/libnfsft_threads.la
else
  LIB_NFSFT_THREADS=
endif
else
  LIB_NFSFT=
  DIR_NFSFT=
  LIB_NFSFT_THREADS=
endif

if HAVE_NFSOFT
  LIB_NFSOFT=nfsoft/libnfsoft.la
  DIR_NFSOFT=nfsoft
if HAVE_THREADS
  LIB_NFSOFT_THREADS=nfsoft/libnfsoft_threads.la
else
  LIB_NFSOFT_THREADS=
endif
else
  LIB_NFSOFT=
  DIR_NFSOFT=
  LIB_NFSOFT_THREADS=
endif

if HAVE_NNFFT
  LIB_NNFFT=nnfft/libnnfft.la
  DIR_NNFFT=nnfft
else
  LIB_NNFFT=
  DIR_NNFFT=
endif

if HAVE_NSFFT
  LIB_NSFFT=nsfft/libnsfft.la
  DIR_NSFFT=nsfft
else
  LIB_NSFFT=
  DIR_NSFFT=
endif

if HAVE_MRI
  LIB_MRI=mri/libmri.la
  DIR_MRI=mri
else
  LIB_MRI=
  DIR_MRI=
endif

if HAVE_FPT
  LIB_FPT=fpt/libfpt.la
  DIR_FPT=fpt
if HAVE_THREADS
  LIB_FPT_THREADS=fpt/libfpt_threads.la
else
  LIB_FPT_THREADS=
endif
else
  LIB_FPT=
  DIR_FPT=
  LIB_FPT_THREADS=
endif

SUBDIRS = util nfft $(DIR_NFCT) $(DIR_NFST) $(DIR_NNFFT) $(DIR_NSFFT) \
  $(DIR_MRI) $(DIR_FPT) $(DIR_NFSFT) $(DIR_NFSOFT) solver .

if HAVE_THREADS
  LIBKERNEL_THREADS_LA = libkernel_threads.la
else
  LIBKERNEL_THREADS_LA =
endif


noinst_LTLIBRARIES = libkernel.la $(LIBKERNEL_THREADS_LA)

libkernel_la_SOURCES =

libkernel_la_LIBADD = util/libutil.la nfft/libnfft.la $(LIB_NFCT) $(LIB_NFST) \
  $(LIB_NNFFT) $(LIB_NSFFT) $(LIB_MRI) $(LIB_FPT) $(LIB_NFSFT) $(LIB_NFSOFT) \
  solver/libsolver.la

if HAVE_THREADS
  libkernel_threads_la_SOURCES =

  libkernel_threads_la_LIBADD = util/libutil_threads.la nfft/libnfft_threads.la $(LIB_NFCT) $(LIB_NFST) \
    $(LIB_NNFFT) $(LIB_NSFFT) $(LIB_MRI) $(LIB_FPT_THREADS) $(LIB_NFSFT_THREADS) $(LIB_NFSOFT_THREADS) \
    solver/libsolver.la

if HAVE_OPENMP
  libkernel_threads_la_CFLAGS = $(OPENMP_CFLAGS)
endif
endif
