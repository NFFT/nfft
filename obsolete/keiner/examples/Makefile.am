AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = nfsft

## The source files.
nfsft_SOURCES = main.c

#  add " -lm" to the following if sfft requires the math library
nfsft_LDADD = -lfftw3 -lm $(top_builddir)/libnfsft.la

