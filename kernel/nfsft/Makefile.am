AM_CPPFLAGS = -I$(top_srcdir)/include

if ENABLE_OPENMP
  LIBNFSFT_THREADS_LA = libnfsft_threads.la
else
  LIBNFSFT_THREADS_LA =
endif

noinst_LTLIBRARIES = libnfsft.la $(LIBNFSFT_THREADS_LA)
libnfsft_la_SOURCES = nfsft.c legendre.h legendre.c api.h

if ENABLE_OPENMP
libnfsft_threads_la_SOURCES = nfsft.c legendre.h legendre.c api.h
libnfsft_threads_la_CFLAGS = $(OPENMP_CFLAGS)
libnfsft_threads_la_LDLAGS = $(OPENMP_CFLAGS)
endif

