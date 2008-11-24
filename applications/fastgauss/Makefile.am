# $Id$

AM_CFLAGS = -I$(top_srcdir)/include

applicationdir = $(pkgdatadir)/applications/fastgauss

application_PROGRAMS = fastgauss

fastgauss_SOURCES = fastgauss.c
fastgauss_LDADD = $(top_builddir)/libnfft3.la @FFTW3_LIBS@ @FFTW3_LIBCMD@

EXTRA_DIST = levelplots.m  show_results.m README
application_DATA = $(EXTRA_DIST)
