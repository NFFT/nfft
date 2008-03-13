# $Id$

AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = 

applicationdir = $(pkgdatadir)/applications/fastsumS2

application_PROGRAMS = fastsumS2

fastsumS2_SOURCES = fastsumS2.c
fastsumS2_LDADD = $(top_builddir)/libnfft3.la -lfftw3 -lm

EXTRA_DIST = fastsumS2.m writeTestcase.m readTestcase.m README
application_DATA = $(EXTRA_DIST)