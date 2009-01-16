# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/include

exampledir = $(pkgdatadir)/examples/bug_nfsft

example_PROGRAMS = bug_nfsft bug2_nfsft

bug_nfsft_SOURCES = bug_nfsft.c
bug2_nfsft_SOURCES = bug2_nfsft.c

bug_nfsft_LDADD = $(top_builddir)/libnfft3.la
bug2_nfsft_LDADD = $(top_builddir)/libnfft3.la

EXTRA_DIST = 
example_DATA = $(EXTRA_DIST)
