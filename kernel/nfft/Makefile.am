AM_CPPFLAGS = -I$(top_srcdir)/include

if ENABLE_OPENMP
  LIBNFFT_THREADS_LA = libnfft_threads.la
else
  LIBNFFT_THREADS_LA =
endif

noinst_LTLIBRARIES = libnfft.la $(LIBNFFT_THREADS_LA)
libnfft_la_SOURCES = nfft.c 

if ENABLE_OPENMP
libnfft_threads_la_SOURCES = nfft.c
libnfft_threads_la_CFLAGS = $(OPENMP_CFLAGS)
libnfft_threads_la_LDLAGS = $(OPENMP_CFLAGS)
endif
