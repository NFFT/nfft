AM_CPPFLAGS = -I$(top_srcdir)/include

if ENABLE_OPENMP
  LIBNFCT_THREADS_LA = libnfct_threads.la
else
  LIBNFCT_THREADS_LA =
endif

noinst_LTLIBRARIES = libnfct.la $(LIBNFCT_THREADS_LA)

libnfct_la_SOURCES = nfct.c

if ENABLE_OPENMP
libnfct_threads_la_SOURCES = nfct.c
libnfct_threads_la_CFLAGS = $(OPENMP_CFLAGS)
endif
