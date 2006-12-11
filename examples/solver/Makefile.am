AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

exampledir = $(pkgdatadir)/examples/solver

## The program name
example_PROGRAMS = simple_test glacier

## The source files.
simple_test_SOURCES = simple_test.c
simple_test_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

glacier_SOURCES = glacier.c
glacier_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

## Include these extra files in the distribution.
EXTRA_DIST = glacier.m glacier_cv.m vol87.dat README
example_DATA = $(EXTRA_DIST)
