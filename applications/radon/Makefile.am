# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/include

applicationdir = $(bindir)/nfft/applications/radon

application_PROGRAMS = radon inverse_radon

radon_SOURCES = radon.c
radon_LDADD = $(top_builddir)/libnfft3.la @fftw3_LDFLAGS@ @fftw3_LIBS@

inverse_radon_SOURCES = inverse_radon.c
inverse_radon_LDADD = $(top_builddir)/libnfft3.la @fftw3_LDFLAGS@ @fftw3_LIBS@

EXTRA_DIST = phantom.m radon.m ridgelet.m README
application_DATA = $(EXTRA_DIST)