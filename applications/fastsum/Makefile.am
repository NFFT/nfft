# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/include

noinst_LTLIBRARIES = libfastsum.la libkernels.la
noinst_PROGRAMS = fastsum_test fastsum_matlab

libfastsum_la_SOURCES = fastsum.c fastsum.h kernels.h
libfastsum_la_LIBADD = @fftw3_LIBS@
libfastsum_la_LDFLAGS = @fftw3_LDFLAGS@

libkernels_la_SOURCES = kernels.c

fastsum_test_SOURCES = fastsum_test.c
fastsum_test_LDADD = libfastsum.la libkernels.la $(top_builddir)/libnfft3.la

fastsum_matlab_SOURCES = fastsum_matlab.c
fastsum_matlab_LDADD = libfastsum.la libkernels.la $(top_builddir)/libnfft3.la

EXTRA_DIST = fastsum.m fastsum_test.m README
