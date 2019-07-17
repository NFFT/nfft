.PHONY: libnfsftjulia-link clean-libnfsftjulia-link

# compiler flags
AM_CPPFLAGS = -I$(top_srcdir)/include

# library
lib_LTLIBRARIES = libnfsftjulia.la
libnfsftjulia_la_SOURCES = libnfsftjulia.c

if HAVE_THREADS
  libadd_for_fftw_threads=@fftw3_threads_LIBS@
else
  libadd_for_fftw_threads=
endif

libnfsftjulia_la_LIBADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_julia.la @fftw3_LDFLAGS@ $(libadd_for_fftw_threads) @fftw3_LIBS@

libnfsftjulia_la_LDFLAGS = -no-undefined -module -shared -avoid-version @fftw3_LDFLAGS@

EXTRA_DIST = NFSFT.jl simple_test.jl

libnfsftjulia-link: all-am
	soname=`$(EGREP) "^dlname=" libnfsftjulia.la | $(SED) -e "s|^dlname='\(.*\)'|\1|"` ; \
	rm -f "$$soname"; \
	$(LN_S) ".libs/$$soname" "./$$soname"

clean-libnfsftjulia-link:
	soname=`$(EGREP) "^dlname=" libnfsftjulia.la | $(SED) -e "s|^dlname='\(.*\)'|\1|"`; \
	rm -f "$$soname"

all: all-am libnfsftjulia-link

clean: clean-libnfsftjulia-link clean-am
