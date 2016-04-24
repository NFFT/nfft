AM_CPPFLAGS = -I$(top_srcdir)/include

if HAVE_THREADS
  LIBFPT_THREADS_LA = libfpt_threads.la
else
  LIBFPT_THREADS_LA =
endif

noinst_LTLIBRARIES = libfpt.la $(LIBFPT_THREADS_LA)

libfpt_la_SOURCES = fpt.c

noinst_HEADERS = fpt.h

if HAVE_THREADS
  libfpt_threads_la_SOURCES = fpt.c
if HAVE_OPENMP
  libfpt_threads_la_CFLAGS = $(OPENMP_CFLAGS)
endif
endif

