# $Id$

AM_CFLAGS = -I$(top_srcdir)/include

exampledir = $(pkgdatadir)/examples/nnfft/accuracy

example_PROGRAMS = accuracy

accuracy_SOURCES = accuracy.c
accuracy_LDADD = $(top_builddir)/libnfft3.la -lfftw3
accuracy_LDFLAGS = -all-static
