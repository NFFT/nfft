AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = cgne

## The source files.
cgne_SOURCES = cgne.c

cgne_LDADD = $(top_builddir)/libnfsft.la -lm -lfftw3