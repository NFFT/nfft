# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/include

exampledir = $(pkgdatadir)/examples/nnfft

example_PROGRAMS = simple_test accuracy

simple_test_SOURCES = simple_test.c
simple_test_LDADD = $(top_builddir)/libnfft3.la

accuracy_SOURCES = accuracy.c
accuracy_LDADD = $(top_builddir)/libnfft3.la

EXTRA_DIST = accuracy.m README
example_DATA = $(EXTRA_DIST)
