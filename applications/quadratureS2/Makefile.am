AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = quadratureS2

## The source files.
quadratureS2_SOURCES = quadratureS2.c
quadratureS2_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

EXTRA_DIST = quadratureS2.m writeTestcase.m writeWeights.m
