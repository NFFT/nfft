AM_CPPFLAGS = -I$(top_srcdir)/include @nfft_benchmarks_CPPFLAGS@ @CPPFLAGS@

if HAVE_BENCHMARKS
if BENCHMARK_AGNOSTIC_WINDOW
NFFT_BENCHMARKS_DIRECT = bench_nfft_direct
bench_nfft_direct_SOURCES = bench_nfft_direct.cpp util.h
bench_nfft_direct_CXXFLAGS = @nfft_benchmarks_CXXFLAGS@ @CXXFLAGS@
bench_nfft_direct_LDFLAGS = @nfft_benchmarks_LDFLAGS@ @fftw3_LDFLAGS@ @LDFLAGS@
bench_nfft_direct_LDADD = @nfft_benchmarks_LIBS@ $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@ @LIBS@
if ENABLE_OPENMP
NFFT_BENCHMARKS_DIRECT_OMP = bench_nfft_direct_omp
bench_nfft_direct_omp_SOURCES = bench_nfft_direct.cpp util.h
bench_nfft_direct_omp_CXXFLAGS = @nfft_benchmarks_CXXFLAGS@ $(OPENMP_CFLAGS) @CXXFLAGS@
bench_nfft_direct_omp_LDFLAGS = @nfft_benchmarks_LDFLAGS@ @fftw3_LDFLAGS@ $(OPENMP_CFLAGS) @LDFLAGS@
bench_nfft_direct_omp_LDADD = @nfft_benchmarks_LIBS@ $(top_builddir)/libnfft3@PREC_SUFFIX@_omp.la @fftw3_LIBS_omp@ @fftw3_LIBS@ $(OPENMP_LIBS)
endif
else
  NFFT_BENCHMARKS_DIRECT =
  NFFT_BENCHMARKS_DIRECT_OMP =
endif
else
  NFFT_BENCHMARKS_DIRECT =
  NFFT_BENCHMARKS_DIRECT_OMP =
endif

# Aggregate all benchmarks
NFFT_BENCHMARKS = $(NFFT_BENCHMARKS_DIRECT)

if ENABLE_OPENMP
NFFT_BENCHMARKS_OMP = $(NFFT_BENCHMARKS_DIRECT_OMP)
endif

noinst_PROGRAMS = $(NFFT_BENCHMARKS) $(NFFT_BENCHMARKS_OMP)

# Run all built benchmarks
bench: $(NFFT_BENCHMARKS) $(NFFT_BENCHMARKS_OMP)
	for bench in $(NFFT_BENCHMARKS); do \
		echo "Running $$bench..."; \
		ulimit -c unlimited; \
		./$$bench; \
	done
	for bench in $(NFFT_BENCHMARKS_OMP); do \
		echo "Running $$bench..."; \
		ulimit -c unlimited; \
		./$$bench; \
	done
