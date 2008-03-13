# $Id$

AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS =

applicationdir = $(pkgdatadir)/applications/quadratureS2

application_PROGRAMS = quadratureS2

quadratureS2_SOURCES = quadratureS2.c
quadratureS2_LDADD = $(top_builddir)/libnfft3.la -lfftw3 -lm

EXTRA_DIST = quadratureS2.m writeTestcase.m writeWeights.m readTestcase.m lgwt.m plotGrid.m
application_DATA = $(EXTRA_DIST)
