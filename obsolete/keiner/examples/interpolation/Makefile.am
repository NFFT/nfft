AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = interpolation

## The source files.
interpolation_SOURCES = interpolation.c

interpolation_LDADD = $(top_builddir)/libnfsft.la -lm -lfftw3