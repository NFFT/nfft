AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = simple_test

## The source files.
simple_test_SOURCES = simple_test.c
simple_test_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3
