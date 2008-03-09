# $Id$

AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -all-static

applicationdir = $(pkgdatadir)/applications/polarFFT

application_PROGRAMS = polar_fft_test mpolar_fft_test linogram_fft_test

polar_fft_test_SOURCES = polar_fft_test.c
polar_fft_test_LDADD = $(top_builddir)/libnfft3.la -lfftw3 -lm

mpolar_fft_test_SOURCES = mpolar_fft_test.c
mpolar_fft_test_LDADD = $(top_builddir)/libnfft3.la -lfftw3 -lm

linogram_fft_test_SOURCES = linogram_fft_test.c
linogram_fft_test_LDADD = $(top_builddir)/libnfft3.la -lfftw3 -lm

EXTRA_DIST = phantom.m fft_test.m README
application_DATA = $(EXTRA_DIST)
