AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = accuracy

## The source files.
accuracy_SOURCES = accuracy.c

accuracy_LDADD = $(top_builddir)/libnfft.la -lm -lfftw3
