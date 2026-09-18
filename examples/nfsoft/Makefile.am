AM_CPPFLAGS = -I$(top_srcdir)/include

## The program name
noinst_PROGRAMS = simple_test

## The source files.
simple_test_SOURCES = simple_test.c
simple_test_LDFLAGS = @fftw3_LDFLAGS@
simple_test_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@

EXTRA_DIST = README
