AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = radon

## The source files.
radon_SOURCES = radon.c

radon_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3
