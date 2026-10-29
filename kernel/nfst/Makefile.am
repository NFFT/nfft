AM_CPPFLAGS = -I$(top_srcdir)/include

if ENABLE_OPENMP
  LIBNFST_THREADS_LA = libnfst_threads.la
else
  LIBNFST_THREADS_LA =
endif

noinst_LTLIBRARIES = libnfst.la $(LIBNFST_THREADS_LA)

libnfst_la_SOURCES = nfst.c

if ENABLE_OPENMP
libnfst_threads_la_SOURCES = nfst.c
libnfst_threads_la_CFLAGS = $(OPENMP_CFLAGS)
endif
