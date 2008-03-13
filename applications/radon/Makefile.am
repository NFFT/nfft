# $Id$

AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS =

applicationdir = $(pkgdatadir)/applications/radon

application_PROGRAMS = radon inverse_radon

radon_SOURCES = radon.c
radon_LDADD = $(top_builddir)/libnfft3.la -lfftw3 -lm

inverse_radon_SOURCES = inverse_radon.c
inverse_radon_LDADD = $(top_builddir)/libnfft3.la -lfftw3 -lm

EXTRA_DIST = phantom.m radon.m ridgelet.m README
application_DATA = $(EXTRA_DIST)