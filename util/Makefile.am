# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/3rdparty/cstripack

if HAVE_THREADS
  LIBUTIL_THREADS_LA = libutil_threads.la
else
  LIBUTIL_THREADS_LA =
endif

noinst_LTLIBRARIES = libutil.la $(LIBUTIL_THREADS_LA)

libutil_la_SOURCES = util.c

if HAVE_THREADS
  libutil_threads_la_SOURCES = util.c
if HAVE_OPENMP
  libutil_threads_la_CFLAGS = $(OPENMP_CFLAGS)
endif
endif

