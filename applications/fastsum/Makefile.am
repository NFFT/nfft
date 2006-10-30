AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

applicationdir = $(pkgdatadir)/applications/fastsum

## The library and program names.
noinst_LTLIBRARIES = libfastsum.la libkernels.la
application_PROGRAMS = fastsum_test fastsum_matlab

## The library header.
pkginclude_HEADERS = fastsum.h kernels.h

## The source files.
libfastsum_la_SOURCES = fastsum.c

libkernels_la_SOURCES = kernels.c

fastsum_test_SOURCES = fastsum_test.c
fastsum_test_LDADD = .libs/libfastsum.a .libs/libkernels.a $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

fastsum_matlab_SOURCES = fastsum_matlab.c
fastsum_matlab_LDADD = .libs/libfastsum.a .libs/libkernels.a $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

## Include these extra files in the distribution.
EXTRA_DIST = fastsum.m fastsum_test.m README
application_DATA = $(EXTRA_DIST)