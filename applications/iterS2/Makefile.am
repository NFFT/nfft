# $Id$

AM_CFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/kernel/nfsft

applicationdir = $(pkgdatadir)/applications/iterS2

application_PROGRAMS = iterS2

iterS2_SOURCES = iterS2.c
iterS2_LDADD = $(top_builddir)/libnfft3.la @FFTW3_LIBS@ @FFTW3_LIBCMD@

EXTRA_DIST = iterS2.m writeTestcase.m readTestcase.m README
application_DATA = $(EXTRA_DIST)