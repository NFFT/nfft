# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/include

if HAVE_THREADS
  SIMPLE_TEST_THREADS=simple_test_threads
else
  SIMPLE_TEST_THREADS=
endif

noinst_PROGRAMS = simple_test $(SIMPLE_TEST_THREADS) ndft_fast taylor_nfft flags nfft_times

if HAVE_THREADS
  simple_test_threads_SOURCES = simple_test_threads.c
  simple_test_threads_LDADD = $(top_builddir)/libnfft3_threads.la @fftw3_LDFLAGS@ @fftw3_threads_LIBS@
endif

simple_test_SOURCES = simple_test.c
simple_test_LDADD = $(top_builddir)/libnfft3.la @fftw3_LDFLAGS@ @fftw3_LIBS@

ndft_fast_SOURCES = ndft_fast.c
ndft_fast_LDADD = $(top_builddir)/libnfft3.la @fftw3_LDFLAGS@ @fftw3_LIBS@

taylor_nfft_SOURCES = taylor_nfft.c
taylor_nfft_LDADD = $(top_builddir)/libnfft3.la @fftw3_LDFLAGS@ @fftw3_LIBS@

flags_SOURCES = flags.c
flags_LDADD = $(top_builddir)/libnfft3.la @fftw3_LDFLAGS@ @fftw3_LIBS@

nfft_times_SOURCES = nfft_times.c
nfft_times_LDADD = $(top_builddir)/libnfft3.la @fftw3_LDFLAGS@ @fftw3_LIBS@

EXTRA_DIST = flags.m ndft_fast.m taylor_nfft.m README
