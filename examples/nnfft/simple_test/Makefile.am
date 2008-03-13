# $Id$

AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS =

exampledir = $(pkgdatadir)/examples/nnfft/simple_test

example_PROGRAMS = simple_test

simple_test_SOURCES = simple_test.c
simple_test_LDADD = $(top_builddir)/libnfft3.la -lfftw3 -lm
