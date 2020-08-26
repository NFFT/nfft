.PHONY: libnfftjulia-link clean-libnfftjulia-link

# compiler flags
AM_CPPFLAGS = -I$(top_srcdir)/include

# library
lib_LTLIBRARIES = libnfftjulia.la
libnfftjulia_la_SOURCES = libnfftjulia.c

if HAVE_THREADS
  libadd_for_fftw_threads=@fftw3_threads_LIBS@
else
  libadd_for_fftw_threads=
endif

libnfftjulia_la_LIBADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_julia.la @fftw3_LDFLAGS@ $(libadd_for_fftw_threads) @fftw3_LIBS@

libnfftjulia_la_LDFLAGS = -no-undefined -shared -avoid-version @fftw3_LDFLAGS@

EXTRA_DIST = NFFT.jl simple_test_1d.jl simple_test_2d.jl simple_test_3d.jl README.md

libnfftjulia-link: all-am
	soname=`$(EGREP) "^dlname=" libnfftjulia.la | $(SED) -e "s|^dlname='\(.*\)'|\1|"` ; \
	rm -f "$$soname"; \
	$(LN_S) ".libs/$$soname" "./$$soname"

clean-libnfftjulia-link:
	soname=`$(EGREP) "^dlname=" libnfftjulia.la | $(SED) -e "s|^dlname='\(.*\)'|\1|"`; \
	rm -f "$$soname"

all: all-am libnfftjulia-link

clean: clean-libnfftjulia-link clean-am

#all-local:: libnfftjulia-link
	
#clean-local:: clean-libnfftjulia-link
