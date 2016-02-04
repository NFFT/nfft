AM_CPPFLAGS = -I$(top_srcdir)/include

if HAVE_THREADS
  SIMPLE_TEST_THREADS=simple_test_threads
else
  SIMPLE_TEST_THREADS=
endif

if HAVE_THREADS
if HAVE_OPENMP
  NFFT_BENCHOMP_PROGS=nfft_benchomp nfft_benchomp_createdataset nfft_benchomp_detail_single nfft_benchomp_detail_threads
else
  NFFT_BENCHOMP_PROGS=
endif
else
  NFFT_BENCHOMP_PROGS=
endif

noinst_PROGRAMS = simple_test $(SIMPLE_TEST_THREADS) ndft_fast taylor_nfft flags nfft_times $(NFFT_BENCHOMP_PROGS)

if HAVE_THREADS
  simple_test_threads_SOURCES = simple_test_threads.c
  simple_test_threads_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_threads.la @fftw3_LDFLAGS@ @fftw3_threads_LIBS@
endif

simple_test_SOURCES = simple_test.c
simple_test_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LDFLAGS@ @fftw3_LIBS@

ndft_fast_SOURCES = ndft_fast.c
ndft_fast_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LDFLAGS@ @fftw3_LIBS@

taylor_nfft_SOURCES = taylor_nfft.c
taylor_nfft_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LDFLAGS@ @fftw3_LIBS@

flags_SOURCES = flags.c
flags_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LDFLAGS@ @fftw3_LIBS@

nfft_times_SOURCES = nfft_times.c
nfft_times_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LDFLAGS@ @fftw3_LIBS@

if HAVE_THREADS
if HAVE_OPENMP
  nfft_benchomp_SOURCES = nfft_benchomp.c
  nfft_benchomp_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_threads.la @fftw3_LDFLAGS@ @fftw3_threads_LIBS@
  nfft_benchomp_CFLAGS = $(OPENMP_CFLAGS)

  nfft_benchomp_createdataset_SOURCES = nfft_benchomp_createdataset.c
  nfft_benchomp_createdataset_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LDFLAGS@ @fftw3_LIBS@

  nfft_benchomp_detail_single_SOURCES = nfft_benchomp_detail.c
  nfft_benchomp_detail_single_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LDFLAGS@ @fftw3_LIBS@

  nfft_benchomp_detail_threads_SOURCES = nfft_benchomp_detail.c
  nfft_benchomp_detail_threads_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_threads.la @fftw3_LDFLAGS@ @fftw3_threads_LIBS@
  nfft_benchomp_detail_threads_CFLAGS = $(OPENMP_CFLAGS)
endif
endif

EXTRA_DIST = flags.m ndft_fast.m taylor_nfft.m README
