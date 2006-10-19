## Add top-level include directory to compiler search path
AM_CFLAGS = -I$(top_srcdir)/include

## The program name
noinst_PROGRAMS = fastsumS2

## The source files.
fastsumS2_SOURCES = fastsumS2.c

## The libraries to link against
fastsumS2_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

## Include these extra files in the distribution.
EXTRA_DIST = fastsumS2.m writeTestcase.m readTestcase.m README doxygen.c
