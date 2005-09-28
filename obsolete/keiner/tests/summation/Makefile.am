AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = summation

## The source files.
summation_SOURCES = summation.c

summation_LDADD = $(top_builddir)/libnfsft.la -lm -lfftw3 -lgsl