if HAVE_NFSFT
  DIR_NFSFT=nfsft
else
  DIR_NFSFT=
endif

if HAVE_NFSOFT
  DIR_NFSOFT=nfsoft
else
  DIR_NFSOFT=
endif

DIR_NFFT=nfft
DIR_NNFFT=nnfft

SUBDIRS = . $(DIR_NFFT) $(DIR_NFSFT) $(DIR_NFSOFT) $(DIR_NNFFT)

AM_CPPFLAGS = -I$(top_srcdir)/include $(matlab_CPPFLAGS)

noinst_LTLIBRARIES = libmatlab.la

libmatlab_la_SOURCES = args.c malloc.c imex.h

if HAVE_MATLAB_THREADS
  libmatlab_la_CFLAGS = $(OPENMP_CFLAGS)
endif

EXTRA_DIST = README
