AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = func_rec_2d

## The source files.
func_rec_2d_SOURCES = func_rec_2d.c

func_rec_2d_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3
