AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The program names
noinst_PROGRAMS = fastsum_test

## The source files
fastsum_test_SOURCES = fastsum_test.c

fastsum_test_LDADD = $(top_builddir)/highlevel/fastsum/libfastsum.la $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

## extra files to include in distribution
# Jens: File not in repository!
# EXTRA_DIST = fastsum_test.m
