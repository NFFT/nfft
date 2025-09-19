AM_CPPFLAGS = -I$(top_srcdir)/include

if ENABLE_OPENMP
  SIMPLE_TEST_THREADS=simple_test_threads
else
  SIMPLE_TEST_THREADS=
endif

if ENABLE_OPENMP
NFSFT_BENCHOMP_PROGS=nfsft_benchomp nfsft_benchomp_createdataset nfsft_benchomp_detail_single nfsft_benchomp_detail_threads
else
  NFSFT_BENCHOMP_PROGS=
endif

## The program name
noinst_PROGRAMS = simple_test $(SIMPLE_TEST_THREADS) $(NFSFT_BENCHOMP_PROGS)

## The source files.
simple_test_SOURCES = simple_test.c
simple_test_LDFLAGS = @fftw3_LDFLAGS@
simple_test_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@

if ENABLE_OPENMP
  simple_test_threads_SOURCES = simple_test_threads.c
  simple_test_threads_CFLAGS = $(OPENMP_CFLAGS)
  simple_test_threads_LDFLAGS = $(OPENMP_CFLAGS) @fftw3_LDFLAGS@
  simple_test_threads_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_omp.la @fftw3_LIBS_omp@ @fftw3_LIBS@ 
endif

if ENABLE_OPENMP
nfsft_benchomp_SOURCES = nfsft_benchomp.c
nfsft_benchomp_CFLAGS = $(OPENMP_CFLAGS)
nfsft_benchomp_LDFLAGS = $(OPENMP_LDFLAGS) @fftw3_LDFLAGS@
nfsft_benchomp_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_omp.la @fftw3_LIBS_omp@ @fftw3_LIBS@ 

nfsft_benchomp_createdataset_SOURCES = nfsft_benchomp_createdataset.c
nfsft_benchomp_createdataset_LDFLAGS = @fftw3_LDFLAGS@
nfsft_benchomp_createdataset_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@

nfsft_benchomp_detail_single_SOURCES = nfsft_benchomp_detail.c
nfsft_benchomp_detail_single_LDFLAGS = @fftw3_LDFLAGS@
nfsft_benchomp_detail_single_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@

nfsft_benchomp_detail_threads_SOURCES = nfsft_benchomp_detail.c
nfsft_benchomp_detail_threads_CFLAGS = $(OPENMP_CFLAGS)
nfsft_benchomp_detail_threads_LDFLAGS = $(OPENMP_CFLAGS) @fftw3_LDFLAGS@
nfsft_benchomp_detail_threads_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_omp.la @fftw3_LIBS_omp@ @fftw3_LIBS@ 
endif
