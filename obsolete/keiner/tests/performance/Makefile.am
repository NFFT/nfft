AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = performance

## The source files.
performance_SOURCES = performance.c

performance_LDADD = $(top_builddir)/libnfsft.la -lm -lfftw3