AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = egm96

## The source files.
egm96_SOURCES = egm96.c

egm96_LDADD = $(top_builddir)/libnfsft.la -lm -lfftw3