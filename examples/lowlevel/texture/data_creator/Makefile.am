AM_CFLAGS = "-Wall -I$(top_srcdir)/include -I../texture_util"

## The program name
noinst_PROGRAMS = data_creator

## The source files for the library
data_creator_SOURCES = data_creator.c

## Additional libraries for the library
data_creator_LDADD = ../texture_util/libtexture_util.la $(top_builddir)/libnfft3.la -lm -lfftw3
