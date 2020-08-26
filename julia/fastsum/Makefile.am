.PHONY: libfastsumjulia-link clean-libfastsumjulia-link

# compiler flags
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/applications/fastsum

# library
lib_LTLIBRARIES = libfastsumjulia.la
libfastsumjulia_la_SOURCES = libfastsumjulia.c

if HAVE_THREADS
  libadd_for_fftw_threads=@fftw3_threads_LIBS@
  libadd_libfastsum=$(top_builddir)/applications/fastsum/libfastsum_threads.la
else
  libadd_for_fftw_threads=
  libadd_libfastsum=$(top_builddir)/applications/fastsum/libfastsum.la
endif
# fastsum_test_LDADD = libfastsum.la libkernels.la $(top_builddir)/libnfft3@PREC_SUFFIX@.la
libfastsumjulia_la_LIBADD = $(libadd_libfastsum) $(top_builddir)/applications/fastsum/libkernels.la $(top_builddir)/libnfft3@PREC_SUFFIX@_julia.la @fftw3_LDFLAGS@ $(libadd_for_fftw_threads) @fftw3_LIBS@


libfastsumjulia_la_LDFLAGS = -no-undefined -shared -avoid-version @fftw3_LDFLAGS@

EXTRA_DIST = fastsum.jl simple_test.jl

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
