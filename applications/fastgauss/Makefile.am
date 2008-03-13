# $Id$

AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = 

applicationdir = $(pkgdatadir)/applications/fastgauss

application_PROGRAMS = fastgauss

fastgauss_SOURCES = fastgauss.c
fastgauss_LDADD = $(top_builddir)/libnfft3.la -lfftw3 -lm

EXTRA_DIST = levelplots.m  show_results.m README
application_DATA = $(EXTRA_DIST)
