.PHONY: libnfstjulia-link clean-libnfstjulia-link

# compiler flags
AM_CPPFLAGS = -I$(top_srcdir)/include

# library
lib_LTLIBRARIES = libnfstjulia.la
libnfstjulia_la_SOURCES = libnfstjulia.c

if ENABLE_OPENMP
  libadd_for_fftw_threads=@fftw3_LIBS_omp@ @fftw3_LIBS@ 
else
  libadd_for_fftw_threads=
endif

libnfstjulia_la_LIBADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_julia.la $(libadd_for_fftw_threads) @fftw3_LIBS@
libnfstjulia_la_LDFLAGS = -no-undefined -shared -avoid-version @fftw3_LDFLAGS@

EXTRA_DIST = NFST.jl simple_test_1d.jl simple_test_2d.jl simple_test_3d.jl

libnfstjulia-link: all-am
	soname=`$(EGREP) "^dlname=" libnfstjulia.la | $(SED) -e "s|^dlname='\(.*\)'|\1|"` ; \
	rm -f "$$soname"; \
	$(LN_S) ".libs/$$soname" "./$$soname"

clean-libnfstjulia-link:
	soname=`$(EGREP) "^dlname=" libnfstjulia.la | $(SED) -e "s|^dlname='\(.*\)'|\1|"`; \
	rm -f "$$soname"

all: all-am libnfstjulia-link

clean: clean-libnfstjulia-link clean-am
