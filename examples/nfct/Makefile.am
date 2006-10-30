AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

exampledir = $(pkgdatadir)/examples/nfct

## The program name
example_PROGRAMS = simple_test

## The source files.
simple_test_SOURCES = simple_test.c
simple_test_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

## Include these extra files in the distribution.
EXTRA_DIST = README
example_DATA = $(EXTRA_DIST)
