## Add top-level include directory to compiler search path
AM_CFLAGS = -I$(top_srcdir)/include

applicationdir = $(pkgdatadir)/applications/fastsumS2

## The program name
application_PROGRAMS = fastsumS2

## The source files.
fastsumS2_SOURCES = fastsumS2.c

## The libraries to link against
fastsumS2_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

## Include these extra files in the distribution.
EXTRA_DIST = fastsumS2.m writeTestcase.m readTestcase.m README
application_DATA = $(EXTRA_DIST)