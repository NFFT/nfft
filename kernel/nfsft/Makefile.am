AM_CPPFLAGS = -I$(top_srcdir)/include

if HAVE_THREADS
  LIBNFSFT_THREADS_LA = libnfsft_threads.la
else
  LIBNFSFT_THREADS_LA =
endif

noinst_LTLIBRARIES = libnfsft.la $(LIBNFSFT_THREADS_LA)
libnfsft_la_SOURCES = nfsft.c legendre.h legendre.c api.h

if HAVE_THREADS
  libnfsft_threads_la_SOURCES = nfsft.c legendre.h legendre.c api.h
if HAVE_OPENMP
  libnfsft_threads_la_CFLAGS = $(OPENMP_CFLAGS)
endif
endif

