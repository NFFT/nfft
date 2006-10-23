AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The program name.
noinst_PROGRAMS = fastgauss

## The source files.
fastgauss_SOURCES = fastgauss.c
fastgauss_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3

## Include these extra files in the distribution.
EXTRA_DIST = levelplots.m  show_results.m README