AM_CFLAGS = -Wall -I$(top_srcdir)/include

## The program name
noinst_LTLIBRARIES = libtexture_util.la
noinst_PROGRAMS = compare

## The source files for the library
libtexture_util_la_SOURCES = texture_util.c
compare_SOURCES = compare.c

## Additional libraries for the library
libtexture_util_la_LIBADD = $(top_builddir)/libnfft3.la -lm -lfftw3
compare_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3 libtexture_util.la
