AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = simple_itest

## The source files.
simple_test_SOURCES = simple_itest.c

simple_test_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3
