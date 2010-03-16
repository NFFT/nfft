# $Id$

LIBTOOL_DEPS = @LIBTOOL_DEPS@

libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck

# Subdirectories
DIST_SUBDIRS= 3rdparty include util kernel . examples applications matlab

if MAINTAINER_MODE
  EXTRA_DIRS= tests
else
  EXTRA_DIRS=
endif

if HAVE_MATLAB
  MATLAB_DIRS= matlab
else
  MATLAB_DIRS=
endif

if HAVE_THREADS
  LIBNFFT3_THREADS_LA = libnfft3_threads.la
else
  LIBNFFT3_THREADS_LA =
endif

SUBDIRS= 3rdparty include util kernel . examples applications $(EXTRA_DIRS) $(MATLAB_DIRS)

lib_LTLIBRARIES = libnfft3.la $(LIBNFFT3_THREADS_LA)

libnfft3_la_SOURCES =
libnfft3_la_LIBADD = 3rdparty/lib3rdparty.la util/libutil.la kernel/libkernel.la @fftw3_LIBS@
libnfft3_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@ @fftw3_LDFLAGS@

if HAVE_THREADS
libnfft3_threads_la_SOURCES =
libnfft3_threads_la_LIBADD = 3rdparty/lib3rdparty.la util/libutil.la kernel/libkernel.la @fftw3_threads_LIBS@ @fftw3_LIBS@
libnfft3_threads_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@ @fftw3_LDFLAGS@
endif

EXTRA_DIST = bootstrap.sh VERSION doxygen.c nfft3.pc.in VERSION.in doxygen/doxygen.Doxyfile.in

install-data-hook:
	mkdir -p $(datadir)/nfft
	cp -R doc $(datadir)/nfft

maintainer-clean-local:
	rm -f doc/api/html/*.html
	rm -f doc/api/html/*.png
	rm -f doc/api/html/*.gif
	rm -f doc/api/html/doxygen.css
	rm -f doc/api/html/formula.repository
	rm -f doc/api/latex/*.tex
	rm -f doc/api/latex/*.sty
	rm -f doc/api/latex/*.ttf
	rm -f doc/api/latex/*.pdf
	rm -f doc/api/latex/Makefile

dist-hook:
	support/toc.sed < doxygen.c
	doxygen doxygen/doxygen.Doxyfile
	rm -f doc/api/html/toc.txt
	support/rem.sed < doc/api/html/index.html > doc/api/html/index.html.new
	mv doc/api/html/index.html.new doc/api/html/index.html
	mkdir -p $(distdir)/doc
	mkdir -p $(distdir)/doc/api
	cp -R doc/api/html $(distdir)/doc/api
	cp -R doc/api/latex $(distdir)/doc/api
	cp -R doc/tutorial $(distdir)/doc
	rm -Rf $(distdir)/doc/.svn
	rm -Rf $(distdir)/doc/tutorial/.svn
	rm -Rf $(distdir)/doc/api/html/.svn
	rm -Rf $(distdir)/doc/api/html/images/.svn
	rm -Rf $(distdir)/doc/api/latex/.svn
	rm -Rf $(distdir)/doc/api/latex/images/.svn
	rm -Rf $(distdir)/doxygen
	rm -Rf $(distdir)/VERSION.in
	cp VERSION $(distdir)

# Flags to be passed to aclocal.
ACLOCAL_AMFLAGS = -I m4
