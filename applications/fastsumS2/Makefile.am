# $Id$

AM_CFLAGS = -I$(top_srcdir)/include

applicationdir = $(pkgdatadir)/applications/fastsumS2

application_PROGRAMS = fastsumS2

fastsumS2_SOURCES = fastsumS2.c
fastsumS2_LDADD = $(top_builddir)/libnfft3.la -lfftw3
fastsumS2_LDFLAGS = -all-static

EXTRA_DIST = fastsumS2.m writeTestcase.m readTestcase.m README
application_DATA = $(EXTRA_DIST)