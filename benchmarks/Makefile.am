AM_CPPFLAGS = -I$(top_srcdir)/include @nfft_benchmarks_CPPFLAGS@ @CPPFLAGS@

if HAVE_BENCHMARKS
  NFFT_BENCHMARKS = bench_nfft_direct
else
  NFFT_BENCHMARKS =
endif

noinst_PROGRAMS = $(NFFT_BENCHMARKS)

# NFFT direct benchmarks
bench_nfft_direct_SOURCES = bench_nfft_direct.cpp util.h
bench_nfft_direct_CXXFLAGS = @nfft_benchmarks_CXXFLAGS@ @CXXFLAGS@
bench_nfft_direct_LDFLAGS = @nfft_benchmarks_LDFLAGS@ @fftw3_LDFLAGS@ @LDFLAGS@
bench_nfft_direct_LDADD = @nfft_benchmarks_LIBS@ $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@ @LIBS@
