AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_builddir)/.libs

exampledir = $(pkgdatadir)/examples/nsfft

## The program name
example_PROGRAMS = simple_test nsfft_test

EXTRA_DIST = README
example_DATA = $(EXTRA_DIST)

## The source files.
simple_test_SOURCES = simple_test.c
simple_test_LDADD = -lm -lfftw3 $(top_builddir)/.libs/libnfft3.a

nsfft_test_SOURCES = nsfft_test.c
nsfft_test_LDADD = -lm -lfftw3 $(top_builddir)/.libs/libnfft3.a
