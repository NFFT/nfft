AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = nfsft

## The source files.
nfsft_SOURCES = main.c

## 
nfsft_LDADD = -lnfsft -lnfft -lfftw3 -lm
