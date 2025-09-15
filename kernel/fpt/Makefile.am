AM_CPPFLAGS = -I$(top_srcdir)/include

if ENABLE_OPENMP
  LIBFPT_THREADS_LA = libfpt_threads.la
else
  LIBFPT_THREADS_LA =
endif

noinst_LTLIBRARIES = libfpt.la $(LIBFPT_THREADS_LA)

libfpt_la_SOURCES = fpt.c

noinst_HEADERS = fpt.h

if ENABLE_OPENMP
libfpt_threads_la_SOURCES = fpt.c
libfpt_threads_la_CFLAGS = $(OPENMP_CFLAGS)
libfpt_threads_la_LDLAGS = $(OPENMP_CFLAGS)
endif

