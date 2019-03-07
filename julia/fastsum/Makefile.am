.PHONY: libfastsumjulia-link clean-libfastsumjulia-link

# compiler flags
AM_CPPFLAGS = -I$(top_srcdir)/include

# library
lib_LTLIBRARIES = libfastsumjulia.la
libfastsumjulia_la_SOURCES = libfastsumjulia.c

if HAVE_THREADS
  libadd_for_fftw_threads=@fftw3_threads_LIBS@
else
  libadd_for_fftw_threads=
endif

libfastsumjulia_la_LIBADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_julia.la @fftw3_LDFLAGS@ $(libadd_for_fftw_threads) @fftw3_LIBS@

libfastsumjulia_la_LDFLAGS = -no-undefined -module -shared -avoid-version @fftw3_LDFLAGS@

EXTRA_DIST = NFFT.jl 

libfastsumjulia-link: all-am
	soname=`$(EGREP) "^dlname=" libfastsumjulia.la | $(SED) -e "s|^dlname='\(.*\)'|\1|"` ; \
	rm -f "$$soname"; \
	$(LN_S) ".libs/$$soname" "./$$soname"

clean-libfastsumjulia-link:
	soname=`$(EGREP) "^dlname=" libfastsumjulia.la | $(SED) -e "s|^dlname='\(.*\)'|\1|"`; \
	rm -f "$$soname"

all: all-am libfastsumjulia-link

clean: clean-libfastsumjulia-link clean-am

#all-local:: libfastsumjulia-link

#clean-local:: clean-libfastsumjulia-link
