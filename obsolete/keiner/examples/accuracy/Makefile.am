AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = accuracy

## The source files.
accuracy_SOURCES = accuracy.c

accuracy_LDADD = $(top_builddir)/libnfsft.la -lm -lfftw3

EXTRA_DIST = gl1.dat gl2.dat gl4.dat gl8.dat gl16.dat gl32.dat gl64.dat gl128.dat gl256.dat gl512.dat gl1024.dat