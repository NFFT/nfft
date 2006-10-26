AM_CFLAGS = -Wall -I$(top_srcdir)/include -I../..

## The program name
nodist_LTLIBRARIES = libtexture_util.la
nodist_PROGRAMS = compare

## The source files for the library
libtexture_util_la_SOURCES = texture_util.c texture_util.h
compare_SOURCES = compare.c

## Additional libraries for the library
#libtexture_util_la_LIBADD = $(top_srcdir)/libnfft3.la ../../libtexture.la -lm -lfftw3
libtexture_util_la_LIBADD = ../../libtexture.la -lm -lfftw3
compare_LDADD = libtexture_util.la

nodistdir = .  
