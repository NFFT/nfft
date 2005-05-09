AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = wrapper

## The source files.
wrapper_SOURCES = wrapper.c

wrapper_LDADD = $(top_builddir)/libnfsft.la -lm -lfftw3