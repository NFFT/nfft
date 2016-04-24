AM_CPPFLAGS = -I$(top_srcdir)/include

if HAVE_THREADS
  SIMPLE_TEST_THREADS=simple_test_threads
else
  SIMPLE_TEST_THREADS=
endif

if HAVE_THREADS
if HAVE_OPENMP
  NFSFT_BENCHOMP_PROGS=nfsft_benchomp nfsft_benchomp_createdataset nfsft_benchomp_detail_single nfsft_benchomp_detail_threads
else
  NFSFT_BENCHOMP_PROGS=
endif
else
  NFSFT_BENCHOMP_PROGS=
endif

## The program name
noinst_PROGRAMS = simple_test $(SIMPLE_TEST_THREADS) $(NFSFT_BENCHOMP_PROGS)

## The source files.
simple_test_SOURCES = simple_test.c
simple_test_LDADD = $(top_builddir)/libnfft3.la @fftw3_LDFLAGS@ @fftw3_LIBS@

if HAVE_THREADS
  simple_test_threads_SOURCES = simple_test_threads.c
  simple_test_threads_LDADD = $(top_builddir)/libnfft3_threads.la @fftw3_LDFLAGS@ @fftw3_threads_LIBS@
if HAVE_OPENMP
  simple_test_threads_CFLAGS = $(OPENMP_CFLAGS)
endif
endif

if HAVE_THREADS
if HAVE_OPENMP
  nfsft_benchomp_SOURCES = nfsft_benchomp.c
  nfsft_benchomp_LDADD = $(top_builddir)/libnfft3_threads.la @fftw3_LDFLAGS@ @fftw3_threads_LIBS@
  nfsft_benchomp_CFLAGS = $(OPENMP_CFLAGS)

  nfsft_benchomp_createdataset_SOURCES = nfsft_benchomp_createdataset.c
  nfsft_benchomp_createdataset_LDADD = $(top_builddir)/libnfft3.la @fftw3_LDFLAGS@ @fftw3_LIBS@

  nfsft_benchomp_detail_single_SOURCES = nfsft_benchomp_detail.c
  nfsft_benchomp_detail_single_LDADD = $(top_builddir)/libnfft3.la @fftw3_LDFLAGS@ @fftw3_LIBS@

  nfsft_benchomp_detail_threads_SOURCES = nfsft_benchomp_detail.c
  nfsft_benchomp_detail_threads_LDADD = $(top_builddir)/libnfft3_threads.la @fftw3_LDFLAGS@ @fftw3_threads_LIBS@
  nfsft_benchomp_detail_threads_CFLAGS = $(OPENMP_CFLAGS)
endif
endif

