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

if HAVE_NFCT
  DIR_NFCT=nfct
else
  DIR_NFCT=
endif

if HAVE_NFST
  DIR_NFST=nfst
else
  DIR_NFST=
endif

if HAVE_NNFFT
  DIR_NNFFT=nnfft
else
  DIR_NNFFT=
endif

if HAVE_APPLICATIONS
  DIR_FASTSUM=fastsum
else
  DIR_FASTSUM=
endif


SUBDIRS = . $(DIR_NFFT) $(DIR_NFCT) $(DIR_NFST) $(DIR_NFSFT) $(DIR_NFSOFT) $(DIR_NNFFT) $(DIR_FASTSUM)

AM_CPPFLAGS = -I$(top_srcdir)/include $(matlab_CPPFLAGS)

noinst_LTLIBRARIES = libmatlab.la

libmatlab_la_SOURCES = args.c malloc.c imex.h

if HAVE_MATLAB_THREADS
  libmatlab_la_CFLAGS = $(OPENMP_CFLAGS)
endif

EXTRA_DIST = README
