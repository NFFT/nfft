AM_CPPFLAGS = -I$(top_srcdir)/include

if ENABLE_OPENMP
  SIMPLE_TEST_THREADS=simple_test_threads
else
  SIMPLE_TEST_THREADS=
endif

if ENABLE_OPENMP
NFFT_BENCHOMP_PROGS=nfft_benchomp nfft_benchomp_createdataset nfft_benchomp_detail_single nfft_benchomp_detail_threads
else
  NFFT_BENCHOMP_PROGS=
endif

noinst_PROGRAMS = simple_test $(SIMPLE_TEST_THREADS) ndft_fast taylor_nfft flags nfft_times $(NFFT_BENCHOMP_PROGS) nfast_native

if ENABLE_OPENMP
  simple_test_threads_SOURCES = simple_test_threads.c
  simple_test_threads_LDFLAGS = @fftw3_LDFLAGS@
  simple_test_threads_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_omp.la @fftw3_LIBS_omp@ @fftw3_LIBS@ $(OPENMP_LIBS)
endif

simple_test_SOURCES = simple_test.c
simple_test_LDFLAGS = @fftw3_LDFLAGS@
simple_test_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@

ndft_fast_SOURCES = ndft_fast.c
ndft_fast_LDFLAGS = @fftw3_LDFLAGS@
ndft_fast_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@

taylor_nfft_SOURCES = taylor_nfft.c
taylor_nfft_LDFLAGS = @fftw3_LDFLAGS@
taylor_nfft_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@

flags_SOURCES = flags.c
flags_LDFLAGS = @fftw3_LDFLAGS@
flags_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@

nfft_times_SOURCES = nfft_times.c
nfft_times_LDFLAGS = @fftw3_LDFLAGS@
nfft_times_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@

nfast_native_SOURCES = nfast_native.c
nfast_native_CPPFLAGS = $(AM_CPPFLAGS) -DNFAST_NATIVE_DATA=\"$(abs_srcdir)/data/nfft_1d_8192_128.txt\"
nfast_native_LDFLAGS = @fftw3_LDFLAGS@
nfast_native_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@

if ENABLE_OPENMP
nfft_benchomp_SOURCES = nfft_benchomp.c
nfft_benchomp_CFLAGS = $(OPENMP_CFLAGS)
nfft_benchomp_LDFLAGS = $(OPENMP_CFLAGS) @fftw3_LDFLAGS@
nfft_benchomp_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_omp.la @fftw3_LIBS_omp@ @fftw3_LIBS@ $(OPENMP_LIBS)

nfft_benchomp_createdataset_SOURCES = nfft_benchomp_createdataset.c
nfft_benchomp_createdataset_LDFLAGS = @fftw3_LDFLAGS@
nfft_benchomp_createdataset_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@ $(OPENMP_LIBS) 

nfft_benchomp_detail_single_SOURCES = nfft_benchomp_detail.c
nfft_benchomp_detail_single_LDFLAGS = @fftw3_LDFLAGS@
nfft_benchomp_detail_single_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@ $(OPENMP_LIBS) 

nfft_benchomp_detail_threads_SOURCES = nfft_benchomp_detail.c
nfft_benchomp_detail_threads_CFLAGS = $(OPENMP_CFLAGS)
nfft_benchomp_detail_threads_LDFLAGS = $(OPENMP_CFLAGS) @fftw3_LDFLAGS@
nfft_benchomp_detail_threads_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_omp.la @fftw3_LIBS_omp@ @fftw3_LIBS@ $(OPENMP_LIBS) 
endif

EXTRA_DIST = flags.m ndft_fast.m taylor_nfft.m README data/nfft_1d_8192_128.txt
