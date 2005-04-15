AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = nfsft

## The source files.
nfsft_SOURCES = main.c

## 
nfsft_LDADD = -lm -lfftw3 -lnfft $(top_builddir)/libnfsft.la

