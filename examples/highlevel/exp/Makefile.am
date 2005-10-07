AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = exp

## The source files.
exp_SOURCES = exp.c

exp_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3