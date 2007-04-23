# $Id$

## Add top-level include directory to compiler search path
AM_CFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/kernel/nfsft

applicationdir = $(pkgdatadir)/applications/iterS2

## The program name
application_PROGRAMS = iterS2

## The source files.
iterS2_SOURCES = iterS2.c

## The libraries to link against
iterS2_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

## Include these extra files in the distribution.
EXTRA_DIST = iterS2.m writeTestcase.m readTestcase.m README
application_DATA = $(EXTRA_DIST)