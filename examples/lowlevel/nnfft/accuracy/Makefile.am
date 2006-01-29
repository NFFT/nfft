AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = accuracy

## The source files.
accuracy_SOURCES = accuracy.c

accuracy_LDADD = -lm -lfftw3 -lnfft3
