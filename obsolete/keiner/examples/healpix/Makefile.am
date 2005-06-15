AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = healpix

## The source files.
healpix_SOURCES = healpix.c

healpix_LDADD = $(top_builddir)/libnfsft.la -lm -lfftw3