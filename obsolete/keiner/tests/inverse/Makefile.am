AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = inverse

## The source files.
inverse_SOURCES = inverse.c

inverse_LDADD = $(top_builddir)/libnfsft.la -lm -lfftw3

EXTRA_DIST =  