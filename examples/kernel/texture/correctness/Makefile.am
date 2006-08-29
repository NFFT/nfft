AM_CFLAGS = -Wall -I$(top_srcdir)/include -I$(top_srcdir)/lowlevel/nfsft/obsolete -I../texture_util

## The program name
noinst_PROGRAMS = textureTest

## The source files for the library
textureTest_SOURCES = textureTest.c

## Additional libraries for the library
textureTest_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3 ../texture_util/libtexture_util.la
