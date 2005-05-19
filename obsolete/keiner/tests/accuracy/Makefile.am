AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = accuracy

## The source files.
accuracy_SOURCES = accuracy.c

accuracy_LDADD = $(top_builddir)/libnfsft.la -lm -lfftw3

EXTRA_DIST =  