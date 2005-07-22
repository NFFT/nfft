AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = accuracy2

## The source files.
accuracy2_SOURCES = accuracy2.c

accuracy2_LDADD = $(top_builddir)/libnfsft.la -lm -lfftw3

EXTRA_DIST =  