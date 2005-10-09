AM_CFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/lowlevel/nfsft/obsolete
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = simple_test
##simple_test

## The source files.
simple_test_SOURCES = simple_test.c
##simple_test.c

simple_test_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3
