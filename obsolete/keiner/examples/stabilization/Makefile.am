AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = stabilization

## The source files.
stabilization_SOURCES = stabilization.c

template_LDADD = $(top_builddir)/libnfsft.la -lm -lfftw3