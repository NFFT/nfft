# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/include

if HAVE_THREADS
  FASTSUM_TEST_THREADS=fastsum_test_threads
  LIBFASTSUM_THREADS=libfastsum_threads.la
else
  FASTSUM_TEST_THREADS=
  LIBFASTSUM_THREADS=
endif

noinst_LTLIBRARIES = libfastsum.la libkernels.la $(LIBFASTSUM_THREADS)
noinst_PROGRAMS = fastsum_test fastsum_matlab $(FASTSUM_TEST_THREADS)

libfastsum_la_SOURCES = fastsum.c fastsum.h kernels.h
libfastsum_la_LIBADD = @fftw3_LIBS@
libfastsum_la_LDFLAGS = @fftw3_LDFLAGS@

libkernels_la_SOURCES = kernels.c

if HAVE_THREADS
  libfastsum_threads_la_SOURCES = fastsum.c fastsum.h kernels.h
  libfastsum_threads_la_LIBADD = @fftw3_threads_LIBS@
  libfastsum_threads_la_LDFLAGS = @fftw3_threads_LDFLAGS@
  libfastsum_threads_la_CFLAGS = $(OPENMP_CFLAGS)
endif

fastsum_test_SOURCES = fastsum_test.c
fastsum_test_LDADD = libfastsum.la libkernels.la $(top_builddir)/libnfft3.la

fastsum_matlab_SOURCES = fastsum_matlab.c
fastsum_matlab_LDADD = libfastsum.la libkernels.la $(top_builddir)/libnfft3.la

if HAVE_THREADS
  fastsum_test_threads_SOURCES = fastsum_test_threads.c
  fastsum_test_threads_LDADD = libfastsum_threads.la libkernels.la $(top_builddir)/libnfft3_threads.la
  fastsum_test_threads_CFLAGS = $(OPENMP_CFLAGS)
endif

EXTRA_DIST = fastsum.m fastsum_test.m README
