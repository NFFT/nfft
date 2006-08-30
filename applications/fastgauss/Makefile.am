AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

EXTRA_DIST = levelplots.m  output_error.m  output_error_p.m  show_results.m README

## The program name
noinst_PROGRAMS = fastgauss

## The source files
fastgauss_SOURCES = fastgauss.c
fastgauss_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3