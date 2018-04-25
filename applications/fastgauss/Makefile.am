AM_CPPFLAGS = -I$(top_srcdir)/include

noinst_PROGRAMS = fastgauss

fastgauss_SOURCES = fastgauss.c
fastgauss_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LDFLAGS@ @fftw3_LIBS@

EXTRA_DIST = levelplots.m  show_results.m README output_error.m output_error_p.m
