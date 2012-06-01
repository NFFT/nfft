AM_CPPFLAGS = -I$(top_srcdir)/include

if HAVE_THREADS
  LIBNFFT_THREADS_LA = libnfft_threads.la
else
  LIBNFFT_THREADS_LA =
endif

noinst_LTLIBRARIES = libnfft.la $(LIBNFFT_THREADS_LA)
libnfft_la_SOURCES = nfft.c 

if HAVE_THREADS
  libnfft_threads_la_SOURCES = nfft.c
if HAVE_OPENMP
  libnfft_threads_la_CFLAGS = $(OPENMP_CFLAGS)
endif
endif
