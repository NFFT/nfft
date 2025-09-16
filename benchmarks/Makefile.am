AM_CPPFLAGS = -I$(top_srcdir)/include @nfft_benchmarks_CPPFLAGS@ @CPPFLAGS@

if HAVE_BENCHMARKS
if BENCHMARK_AGNOSTIC_WINDOW
NFFT_BENCHMARKS_DIRECT = bench_nfft_direct
bench_nfft_direct_SOURCES = bench_nfft_direct.cpp util.h
bench_nfft_direct_CXXFLAGS = @nfft_benchmarks_CXXFLAGS@ @CXXFLAGS@
bench_nfft_direct_LDFLAGS = @nfft_benchmarks_LDFLAGS@ @fftw3_LDFLAGS@ @LDFLAGS@
bench_nfft_direct_LDADD = @nfft_benchmarks_LIBS@ $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@ @LIBS@
else
  NFFT_BENCHMARKS_DIRECT =
endif
else
  NFFT_BENCHMARKS_DIRECT =
endif

# Aggregate all benchmarks
NFFT_BENCHMARKS = $(NFFT_BENCHMARKS_DIRECT)

noinst_PROGRAMS = $(NFFT_BENCHMARKS)

# Run all built benchmarks
bench: $(NFFT_BENCHMARKS)
	for bench in $(NFFT_BENCHMARKS); do \
		echo "Running $$bench..."; \
		./$$bench; \
	done
