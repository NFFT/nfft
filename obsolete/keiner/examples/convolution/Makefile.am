AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = convolution

## The source files.
convolution_SOURCES = convolution.c

convolution_LDADD = $(top_builddir)/libnfsft.la -lm -lfftw3