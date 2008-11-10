# $Id$

AM_CFLAGS = -I$(top_srcdir)/include

exampledir = $(pkgdatadir)/examples/nfsft

## The program name
example_PROGRAMS = simple_test

## The source files.
simple_test_SOURCES = simple_test.c
simple_test_LDADD = $(top_builddir)/libnfft3.la -lfftw3
