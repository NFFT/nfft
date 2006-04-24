AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = fastsumS2

## The source files.
fastsumS2_SOURCES = fastsumS2.c
fastsumS2_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3 -lgsl -lgslcblas

EXTRA_DIST = fastsumS2.m writeTestcase.m readTestcase.m
