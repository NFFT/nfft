AM_CPPFLAGS = -I$(top_srcdir)/include @nfft_benchmarks_CPPFLAGS@ @CPPFLAGS@

if HAVE_BENCHMARKS
if BENCHMARK_AGNOSTIC_WINDOW
NFFT_BENCHMARKS_DIRECT = bench_nfft_direct
bench_nfft_direct_SOURCES = bench_nfft_direct.cpp util.h
bench_nfft_direct_CXXFLAGS = @nfft_benchmarks_CXXFLAGS@ @CXXFLAGS@
bench_nfft_direct_LDFLAGS = @nfft_benchmarks_LDFLAGS@ @fftw3_LDFLAGS@ @LDFLAGS@
bench_nfft_direct_LDADD = @nfft_benchmarks_LIBS@ $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@ @LIBS@

NFFT_BENCHMARKS_FAST = bench_nfft_fast
bench_nfft_fast_SOURCES = bench_nfft_fast.cpp util.h
bench_nfft_fast_CXXFLAGS = @nfft_benchmarks_CXXFLAGS@ @CXXFLAGS@
bench_nfft_fast_LDFLAGS = @nfft_benchmarks_LDFLAGS@ @fftw3_LDFLAGS@ @LDFLAGS@
bench_nfft_fast_LDADD = @nfft_benchmarks_LIBS@ $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@ @LIBS@

NFCT_BENCHMARKS_DIRECT = bench_nfct_direct
bench_nfct_direct_SOURCES = bench_nfct_direct.cpp util.h
bench_nfct_direct_CXXFLAGS = @nfft_benchmarks_CXXFLAGS@ @CXXFLAGS@
bench_nfct_direct_LDFLAGS = @nfft_benchmarks_LDFLAGS@ @fftw3_LDFLAGS@ @LDFLAGS@
bench_nfct_direct_LDADD = @nfft_benchmarks_LIBS@ $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@ @LIBS@

NFST_BENCHMARKS_DIRECT = bench_nfst_direct
bench_nfst_direct_SOURCES = bench_nfst_direct.cpp util.h
bench_nfst_direct_CXXFLAGS = @nfft_benchmarks_CXXFLAGS@ @CXXFLAGS@
bench_nfst_direct_LDFLAGS = @nfft_benchmarks_LDFLAGS@ @fftw3_LDFLAGS@ @LDFLAGS@
bench_nfst_direct_LDADD = @nfft_benchmarks_LIBS@ $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@ @LIBS@
if ENABLE_OPENMP
NFFT_BENCHMARKS_DIRECT_OMP = bench_nfft_direct_omp
bench_nfft_direct_omp_SOURCES = bench_nfft_direct.cpp util.h
bench_nfft_direct_omp_CXXFLAGS = @nfft_benchmarks_CXXFLAGS@ $(OPENMP_CFLAGS) @CXXFLAGS@
bench_nfft_direct_omp_LDFLAGS = @nfft_benchmarks_LDFLAGS@ @fftw3_LDFLAGS@ $(OPENMP_CFLAGS) @LDFLAGS@
bench_nfft_direct_omp_LDADD = @nfft_benchmarks_LIBS@ $(top_builddir)/libnfft3@PREC_SUFFIX@_omp.la @fftw3_LIBS_omp@ @fftw3_LIBS@ $(OPENMP_LIBS)

NFFT_BENCHMARKS_FAST_OMP = bench_nfft_fast_omp
bench_nfft_fast_omp_SOURCES = bench_nfft_fast.cpp util.h
bench_nfft_fast_omp_CXXFLAGS = @nfft_benchmarks_CXXFLAGS@ $(OPENMP_CFLAGS) @CXXFLAGS@
bench_nfft_fast_omp_LDFLAGS = @nfft_benchmarks_LDFLAGS@ @fftw3_LDFLAGS@ $(OPENMP_CFLAGS) @LDFLAGS@
bench_nfft_fast_omp_LDADD = @nfft_benchmarks_LIBS@ $(top_builddir)/libnfft3@PREC_SUFFIX@_omp.la @fftw3_LIBS_omp@ @fftw3_LIBS@ $(OPENMP_LIBS)

NFCT_BENCHMARKS_DIRECT_OMP = bench_nfct_direct_omp
bench_nfct_direct_omp_SOURCES = bench_nfct_direct.cpp util.h
bench_nfct_direct_omp_CXXFLAGS = @nfft_benchmarks_CXXFLAGS@ $(OPENMP_CFLAGS) @CXXFLAGS@
bench_nfct_direct_omp_LDFLAGS = @nfft_benchmarks_LDFLAGS@ @fftw3_LDFLAGS@ $(OPENMP_CFLAGS) @LDFLAGS@
bench_nfct_direct_omp_LDADD = @nfft_benchmarks_LIBS@ $(top_builddir)/libnfft3@PREC_SUFFIX@_omp.la @fftw3_LIBS_omp@ @fftw3_LIBS@ $(OPENMP_LIBS)

NFST_BENCHMARKS_DIRECT_OMP = bench_nfst_direct_omp
bench_nfst_direct_omp_SOURCES = bench_nfst_direct.cpp util.h
bench_nfst_direct_omp_CXXFLAGS = @nfft_benchmarks_CXXFLAGS@ $(OPENMP_CFLAGS) @CXXFLAGS@
bench_nfst_direct_omp_LDFLAGS = @nfft_benchmarks_LDFLAGS@ @fftw3_LDFLAGS@ $(OPENMP_CFLAGS) @LDFLAGS@
bench_nfst_direct_omp_LDADD = @nfft_benchmarks_LIBS@ $(top_builddir)/libnfft3@PREC_SUFFIX@_omp.la @fftw3_LIBS_omp@ @fftw3_LIBS@ $(OPENMP_LIBS)
if HAVE_NNFFT
NNFFT_BENCHMARKS = bench_nnfft_direct bench_nnfft_fast
bench_nnfft_direct_SOURCES = bench_nnfft_direct.cpp util.h
bench_nnfft_direct_CXXFLAGS = @nfft_benchmarks_CXXFLAGS@ @CXXFLAGS@
bench_nnfft_direct_LDFLAGS = @nfft_benchmarks_LDFLAGS@ @fftw3_LDFLAGS@ @LDFLAGS@
bench_nnfft_direct_LDADD = @nfft_benchmarks_LIBS@ $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@ @LIBS@
bench_nnfft_fast_SOURCES = bench_nnfft_fast.cpp util.h
bench_nnfft_fast_CXXFLAGS = @nfft_benchmarks_CXXFLAGS@ @CXXFLAGS@
bench_nnfft_fast_LDFLAGS = @nfft_benchmarks_LDFLAGS@ @fftw3_LDFLAGS@ @LDFLAGS@
bench_nnfft_fast_LDADD = @nfft_benchmarks_LIBS@ $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@ @LIBS@
if ENABLE_OPENMP
NNFFT_BENCHMARKS_OMP = bench_nnfft_direct_omp bench_nnfft_fast_omp
bench_nnfft_direct_omp_SOURCES = bench_nnfft_direct.cpp util.h
bench_nnfft_direct_omp_CXXFLAGS = @nfft_benchmarks_CXXFLAGS@ $(OPENMP_CFLAGS) @CXXFLAGS@
bench_nnfft_direct_omp_LDFLAGS = @nfft_benchmarks_LDFLAGS@ @fftw3_LDFLAGS@ $(OPENMP_CFLAGS) @LDFLAGS@
bench_nnfft_direct_omp_LDADD = @nfft_benchmarks_LIBS@ $(top_builddir)/libnfft3@PREC_SUFFIX@_omp.la @fftw3_LIBS_omp@ @fftw3_LIBS@ $(OPENMP_LIBS)
bench_nnfft_fast_omp_SOURCES = bench_nnfft_fast.cpp util.h
bench_nnfft_fast_omp_CXXFLAGS = @nfft_benchmarks_CXXFLAGS@ $(OPENMP_CFLAGS) @CXXFLAGS@
bench_nnfft_fast_omp_LDFLAGS = @nfft_benchmarks_LDFLAGS@ @fftw3_LDFLAGS@ $(OPENMP_CFLAGS) @LDFLAGS@
bench_nnfft_fast_omp_LDADD = @nfft_benchmarks_LIBS@ $(top_builddir)/libnfft3@PREC_SUFFIX@_omp.la @fftw3_LIBS_omp@ @fftw3_LIBS@ $(OPENMP_LIBS)
endif
endif
endif
else
  NFFT_BENCHMARKS_DIRECT =
  NFFT_BENCHMARKS_DIRECT_OMP =
  NFFT_BENCHMARKS_FAST =
  NFFT_BENCHMARKS_FAST_OMP =
  NFCT_BENCHMARKS_DIRECT =
  NFCT_BENCHMARKS_DIRECT_OMP =
  NFST_BENCHMARKS_DIRECT =
  NFST_BENCHMARKS_DIRECT_OMP =
endif
else
  NFFT_BENCHMARKS_DIRECT =
  NFFT_BENCHMARKS_DIRECT_OMP =
  NFFT_BENCHMARKS_FAST =
  NFFT_BENCHMARKS_FAST_OMP =
  NFCT_BENCHMARKS_DIRECT =
  NFCT_BENCHMARKS_DIRECT_OMP =
  NFST_BENCHMARKS_DIRECT =
  NFST_BENCHMARKS_DIRECT_OMP =
endif

# Aggregate all benchmarks
NFFT_BENCHMARKS = $(NFFT_BENCHMARKS_DIRECT) $(NFFT_BENCHMARKS_FAST) $(NFCT_BENCHMARKS_DIRECT) $(NFST_BENCHMARKS_DIRECT) $(NNFFT_BENCHMARKS)

if ENABLE_OPENMP
NFFT_BENCHMARKS_OMP = $(NFFT_BENCHMARKS_DIRECT_OMP) $(NFFT_BENCHMARKS_FAST_OMP) $(NFCT_BENCHMARKS_DIRECT_OMP) $(NFST_BENCHMARKS_DIRECT_OMP) $(NNFFT_BENCHMARKS_OMP)
endif

noinst_PROGRAMS = $(NFFT_BENCHMARKS) $(NFFT_BENCHMARKS_OMP)

# Run all built benchmarks
bench: $(NFFT_BENCHMARKS) $(NFFT_BENCHMARKS_OMP)
	for bench in $(NFFT_BENCHMARKS); do \
		echo "Running $$bench..."; \
		./$$bench; \
	done
	for bench in $(NFFT_BENCHMARKS_OMP); do \
		echo "Running $$bench..."; \
		./$$bench; \
	done
