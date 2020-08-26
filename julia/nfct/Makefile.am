.PHONY: libnfctjulia-link clean-libnfctjulia-link

# compiler flags
AM_CPPFLAGS = -I$(top_srcdir)/include

# library
lib_LTLIBRARIES = libnfctjulia.la
libnfctjulia_la_SOURCES = libnfctjulia.c

if HAVE_THREADS
  libadd_for_fftw_threads=@fftw3_threads_LIBS@
else
  libadd_for_fftw_threads=
endif

libnfctjulia_la_LIBADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_julia.la @fftw3_LDFLAGS@ $(libadd_for_fftw_threads) @fftw3_LIBS@

libnfctjulia_la_LDFLAGS = -no-undefined -shared -avoid-version @fftw3_LDFLAGS@

EXTRA_DIST = NFCT.jl simple_test_1d.jl simple_test_2d.jl simple_test_3d.jl

libnfctjulia-link: all-am
	soname=`$(EGREP) "^dlname=" libnfctjulia.la | $(SED) -e "s|^dlname='\(.*\)'|\1|"` ; \
	rm -f "$$soname"; \
	$(LN_S) ".libs/$$soname" "./$$soname"

clean-libnfctjulia-link:
	soname=`$(EGREP) "^dlname=" libnfctjulia.la | $(SED) -e "s|^dlname='\(.*\)'|\1|"`; \
	rm -f "$$soname"

all: all-am libnfctjulia-link

clean: clean-libnfctjulia-link clean-am
