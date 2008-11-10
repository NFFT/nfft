# $Id$

AM_CFLAGS = -I$(top_srcdir)/include

exampledir = $(pkgdatadir)/examples/nfst

example_PROGRAMS = simple_test

simple_test_SOURCES = simple_test.c
simple_test_LDADD = $(top_builddir)/libnfft3.la -lfftw3

EXTRA_DIST = README
example_DATA = $(EXTRA_DIST)
