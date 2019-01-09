AM_CPPFLAGS = -I$(top_srcdir)/include

if HAVE_THREADS
  LIBNFSOFT_THREADS_LA = libnfsoft_threads.la
else
  LIBNFSOFT_THREADS_LA =
endif

noinst_LTLIBRARIES = libnfsoft.la $(LIBNFSOFT_THREADS_LA)
libnfsoft_la_SOURCES = nfsoft.c wigner.h wigner.c

if HAVE_THREADS
  libnfsoft_threads_la_SOURCES = nfsoft.c wigner.h wigner.c
if HAVE_OPENMP
  libnfsoft_threads_la_CFLAGS = $(OPENMP_CFLAGS)
endif
endif