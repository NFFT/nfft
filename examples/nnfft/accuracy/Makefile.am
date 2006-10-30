AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

exampledir = $(pkgdatadir)/examples/nnfft/accuracy

## The program name
example_PROGRAMS = accuracy

## The source files.
accuracy_SOURCES = accuracy.c

accuracy_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3
