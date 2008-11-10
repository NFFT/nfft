# $Id$

AM_CFLAGS = -I$(top_srcdir)/include

exampledir = $(pkgdatadir)/examples/nsfft

example_PROGRAMS = simple_test nsfft_test

simple_test_SOURCES = simple_test.c
simple_test_LDADD = $(top_builddir)/libnfft3.la -lfftw3

nsfft_test_SOURCES = nsfft_test.c
nsfft_test_LDADD = $(top_builddir)/libnfft3.la -lfftw3

EXTRA_DIST = README
example_DATA = $(EXTRA_DIST)
