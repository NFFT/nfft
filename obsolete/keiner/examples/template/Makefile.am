AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = template

## The source files.
template_SOURCES = template.c

template_LDADD = $(top_builddir)/libnfsft.la -lm -lfftw3