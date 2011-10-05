# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/kernel/nfsft

noinst_PROGRAMS = iterS2

iterS2_SOURCES = iterS2.c doxygen.h
iterS2_LDADD = $(top_builddir)/libnfft3.la @fftw3_LDFLAGS@ @fftw3_LIBS@

EXTRA_DIST = iterS2.m writeTestcase.m readTestcase.m README
