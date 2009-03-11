# $Id$

if HAVE_NFSFT
  DIR_NFSFT=nfsft
else
  DIR_NFSFT=
endif

DIR_NFFT=nfft

SUBDIRS = . $(DIR_NFFT) $(DIR_NFSFT)

AM_CPPFLAGS = -I$(top_srcdir)/include $(matlab_CPPFLAGS)

noinst_LTLIBRARIES = libmatlab.la

libmatlab_la_SOURCES = args.c malloc.c imex.h

EXTRA_DIST = README
