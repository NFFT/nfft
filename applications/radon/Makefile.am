AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

applicationdir = $(pkgdatadir)/applications/radon

## The program names
application_PROGRAMS = radon inverse_radon

## The source files
radon_SOURCES = radon.c
inverse_radon_SOURCES = inverse_radon.c

radon_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3
inverse_radon_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

## extra files to include in distribution
EXTRA_DIST = phantom.m radon.m ridgelet.m README
application_DATA = $(EXTRA_DIST)