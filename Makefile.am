# $Id$

LIBTOOL_DEPS = @LIBTOOL_DEPS@

libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck

# Subdirectories
DIST_SUBDIRS=3rdparty include util kernel . tests examples applications matlab support doxygen

if HAVE_MATLAB
  MATLAB_DIRS=matlab
  LIBNFFT3_MATLAB_LA=libnfft3_matlab.la
else
  MATLAB_DIRS=
  LIBNFFT3_MATLAB_LA=
endif

if HAVE_THREADS
  LIBNFFT3_THREADS_LA = libnfft3_threads.la
else
  LIBNFFT3_THREADS_LA =
endif

SUBDIRS= 3rdparty include util kernel . tests examples applications $(MATLAB_DIRS)

lib_LTLIBRARIES = libnfft3.la $(LIBNFFT3_THREADS_LA)
noinst_LTLIBRARIES = $(LIBNFFT3_MATLAB_LA)

libnfft3_la_SOURCES = doxygen.h
libnfft3_la_LIBADD = 3rdparty/lib3rdparty.la util/libutil.la kernel/libkernel.la @fftw3_LIBS@ -lm
libnfft3_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@ @fftw3_LDFLAGS@

if HAVE_THREADS
libnfft3_threads_la_SOURCES =
libnfft3_threads_la_LIBADD = 3rdparty/lib3rdparty.la util/libutil.la kernel/libkernel.la @fftw3_threads_LIBS@ @fftw3_LIBS@ -lm
libnfft3_threads_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@ @fftw3_LDFLAGS@
endif

if HAVE_MATLAB
libnfft3_matlab_la_SOURCES =
libnfft3_matlab_la_LIBADD = 3rdparty/lib3rdparty.la util/libutil.la kernel/libkernel.la @fftw3_LIBS@ -lm
libnfft3_matlab_la_LDFLAGS = $(fftw3_LDFLAGS)
endif

EXTRA_DIST = bootstrap.sh doxygen.c nfft3.pc.in

pkgconfigdir=$(libdir)/pkgconfig
pkgconfig_DATA = nfft3.pc

#install-data-hook:
#	$(MKDIR_P) $(docdir)
#	cp -R $(srcdir)/doc/* $(docdir)/
#	chmod -R u+rwX,g+rX,g-w,o+rX,o-w $(docdir)

#uninstall-hook:
#	rm -Rf $(docdir)

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
	echo `pwd`
	support/toc.sed < doxygen.c
	doxygen doxygen/doxygen.Doxyfile
	rm -f doc/api/html/toc.txt
	support/rem.sed < doc/api/html/index.html > doc/api/html/index.html.new
	mv doc/api/html/index.html.new doc/api/html/index.html
	cp -R $(srcdir)/doc $(distdir)/
#	mkdir -p $(distdir)/doc
#	mkdir -p $(distdir)/doc/api
#	cp -R doc/api/html $(distdir)/doc/api
#	cp -R doc/api/latex $(distdir)/doc/api
#	cp -R doc/tutorial $(distdir)/doc
#	rm -Rf $(distdir)/doxygen
	rm -Rf `find $(distdir) -name .svn -type d`

# Flags to be passed to aclocal.
ACLOCAL_AMFLAGS = -I m4
