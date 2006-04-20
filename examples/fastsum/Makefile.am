AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The program names
noinst_PROGRAMS = fastsum_test fastsum_matlab

## The source files
dist_fastsum_test_SOURCES = kernels.h kernels.c fastsum_test.c
fastsum_test_LDADD = $(top_builddir)/highlevel/fastsum/libfastsum.la $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

dist_fastsum_matlab_SOURCES = kernels.h kernels.c fastsum_matlab.c
fastsum_matlab_LDADD = $(top_builddir)/highlevel/fastsum/libfastsum.la $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

## extra files to include in distribution
EXTRA_DIST = fastsum_test.m fastsum.m
