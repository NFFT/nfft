# $Id$

AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = 

applicationdir = $(pkgdatadir)/applications/fastsum

noinst_LTLIBRARIES = libfastsum.la libkernels.la
application_PROGRAMS = fastsum_test fastsum_matlab

pkginclude_HEADERS = fastsum.h kernels.h

libfastsum_la_SOURCES = fastsum.c

libkernels_la_SOURCES = kernels.c

fastsum_test_SOURCES = fastsum_test.c
fastsum_test_LDADD = libfastsum.la libkernels.la $(top_builddir)/libnfft3.la -lfftw3 -lm

fastsum_matlab_SOURCES = fastsum_matlab.c
fastsum_matlab_LDADD = libfastsum.la libkernels.la $(top_builddir)/libnfft3.la -lfftw3 -lm

EXTRA_DIST = fastsum.m fastsum_test.m README
application_DATA = $(EXTRA_DIST)