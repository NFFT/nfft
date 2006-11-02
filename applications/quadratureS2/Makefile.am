AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

applicationdir = $(pkgdatadir)/applications/quadratureS2

## The program name
application_PROGRAMS = quadratureS2

## The source files.
quadratureS2_SOURCES = quadratureS2.c
quadratureS2_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

EXTRA_DIST = quadratureS2.m writeTestcase.m writeWeights.m lgwt.m plotGrid.m
application_DATA = $(EXTRA_DIST)
