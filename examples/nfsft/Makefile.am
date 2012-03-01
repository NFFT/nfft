# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/include

if HAVE_THREADS
  SIMPLE_TEST_THREADS=simple_test_threads
else
  SIMPLE_TEST_THREADS=
endif

## The program name
noinst_PROGRAMS = simple_test $(SIMPLE_TEST_THREADS)

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
