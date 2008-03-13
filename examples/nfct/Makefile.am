# $Id$

AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS =

exampledir = $(pkgdatadir)/examples/nfct

example_PROGRAMS = simple_test

simple_test_SOURCES = simple_test.c
simple_test_LDADD = $(top_builddir)/libnfft3.la -lfftw3 -lm

EXTRA_DIST = README
example_DATA = $(EXTRA_DIST)
