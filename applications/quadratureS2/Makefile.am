AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)/.libs -all-static

applicationdir = $(pkgdatadir)/applications/quadratureS2

## The program name
application_PROGRAMS = quadratureS2

## The source files.
quadratureS2_SOURCES = quadratureS2.c
quadratureS2_LDADD = -lnfft3 -lfftw3 -lm

EXTRA_DIST = quadratureS2.m writeTestcase.m writeWeights.m readTestcase.m lgwt.m plotGrid.m
application_DATA = $(EXTRA_DIST)
