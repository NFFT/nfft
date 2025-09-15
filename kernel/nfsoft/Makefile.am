AM_CPPFLAGS = -I$(top_srcdir)/include

if ENABLE_OPENMP
  LIBNFSOFT_THREADS_LA = libnfsoft_threads.la
else
  LIBNFSOFT_THREADS_LA =
endif

noinst_LTLIBRARIES = libnfsoft.la $(LIBNFSOFT_THREADS_LA)
libnfsoft_la_SOURCES = nfsoft.c wigner.h wigner.c

if ENABLE_OPENMP
libnfsoft_threads_la_SOURCES = nfsoft.c wigner.h wigner.c
libnfsoft_threads_la_CFLAGS = $(OPENMP_CFLAGS)
libnfsoft_threads_la_LDLAGS = $(OPENMP_CFLAGS)
endif