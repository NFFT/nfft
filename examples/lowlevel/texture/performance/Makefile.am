AM_CFLAGS = -Wall -I$(top_srcdir)/include -I../texture_util

## The program name
noinst_PROGRAMS = performance

## The source files for the library
performance_SOURCES = performance.c

## Additional libraries for the library
performance_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3 ../texture_util/libtexture_util.la
