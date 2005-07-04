AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = reuter

## The source files.
reuter_SOURCES = reuter.c

reuter_LDADD = $(top_builddir)/libnfsft.la -lm -lfftw3