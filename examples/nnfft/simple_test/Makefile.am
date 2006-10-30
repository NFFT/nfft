AM_CFLAGS = -I$(top_srcdir)/include

exampledir = $(pkgdatadir)/examples/nnfft/simple_test

## The program name
example_PROGRAMS = simple_test

## The source files.
simple_test_SOURCES = simple_test.c

simple_test_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3
