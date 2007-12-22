AM_CFLAGS = -I$(top_srcdir)/include

exampledir = $(pkgdatadir)/examples/nfct

## The program name
example_PROGRAMS = simple_test

## The source files.
simple_test_SOURCES = simple_test.c
simple_test_LDADD = $(top_builddir)/libnfft3.la -lfftw3 -lm 

## Include these extra files in the distribution.
EXTRA_DIST = README
example_DATA = $(EXTRA_DIST)
