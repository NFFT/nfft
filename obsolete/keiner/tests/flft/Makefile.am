AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = flft

## The source files.
flft_SOURCES = flft.c

flft_LDADD = $(top_builddir)/libnfsft.la -lm -lfftw3

EXTRA_DIST =  