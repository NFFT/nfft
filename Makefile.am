# $Id$
LIBTOOL_DEPS = @LIBTOOL_DEPS@

libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck

# Subdirectories
DIST_SUBDIRS=3rdparty include kernel . tests examples applications matlab support doxygen

if HAVE_EXAMPLES
  EXAMPLE_DIRS=examples
else
  EXAMPLE_DIRS=
endif

if HAVE_APPLICATIONS
  APPLICATION_DIRS=applications
else
  APPLICATION_DIRS=
endif

if HAVE_MATLAB
  MATLAB_DIRS=matlab
  LIBNFFT3_MATLAB_LA=libnfft3@PREC_SUFFIX@_matlab.la
else
  MATLAB_DIRS=
  LIBNFFT3_MATLAB_LA=
endif

if HAVE_THREADS
  LIBNFFT3_THREADS_LA = libnfft3@PREC_SUFFIX@_threads.la
else
  LIBNFFT3_THREADS_LA =
endif

SUBDIRS= 3rdparty include kernel . tests $(EXAMPLE_DIRS) $(APPLICATION_DIRS) $(MATLAB_DIRS)

lib_LTLIBRARIES = libnfft3@PREC_SUFFIX@.la $(LIBNFFT3_THREADS_LA)
noinst_LTLIBRARIES = $(LIBNFFT3_MATLAB_LA)

libnfft3@PREC_SUFFIX@_la_SOURCES =
libnfft3@PREC_SUFFIX@_la_LIBADD = 3rdparty/lib3rdparty.la kernel/libkernel.la @fftw3_LIBS@ -lm
libnfft3@PREC_SUFFIX@_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@ @fftw3_LDFLAGS@

if HAVE_THREADS
libnfft3@PREC_SUFFIX@_threads_la_SOURCES =
libnfft3@PREC_SUFFIX@_threads_la_LIBADD = 3rdparty/lib3rdparty.la kernel/libkernel_threads.la @fftw3_threads_LIBS@ -lm
libnfft3@PREC_SUFFIX@_threads_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@ @fftw3_LDFLAGS@
if HAVE_OPENMP
  libnfft3@PREC_SUFFIX@_threads_la_CFLAGS = $(OPENMP_CFLAGS)
endif
endif

if HAVE_MATLAB
  libnfft3@PREC_SUFFIX@_matlab_la_SOURCES =
if HAVE_MATLAB_THREADS
    libnfft3@PREC_SUFFIX@_matlab_la_LIBADD = 3rdparty/lib3rdparty.la kernel/libkernel_threads.la @matlab_fftw3_LIBS@ -lm
else
    libnfft3@PREC_SUFFIX@_matlab_la_LIBADD = 3rdparty/lib3rdparty.la kernel/libkernel.la @matlab_fftw3_LIBS@ -lm
endif
  libnfft3@PREC_SUFFIX@_matlab_la_LDFLAGS = @matlab_fftw3_LDFLAGS@
endif

EXTRA_DIST = bootstrap.sh doxygen.dox nfft3.pc.in

nfft3@PREC_SUFFIX@.pc: nfft3.pc
	cp -f nfft3.pc nfft3@PREC_SUFFIX@.pc
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = nfft3@PREC_SUFFIX@.pc

install-data-hook:
	$(MKDIR_P) $(DESTDIR)$(docdir)
	cp -R $(srcdir)/doc/* $(DESTDIR)$(docdir)/

uninstall-hook:
	rm -Rf $(docdir)

maintainer-clean-local:
	rm -f $(srcdir)/doc/api/html/*.html
	rm -f $(srcdir)/doc/api/html/*.png
	rm -f $(srcdir)/doc/api/html/*.gif
	rm -f $(srcdir)/doc/api/html/doxygen.css
	rm -f $(srcdir)/doc/api/html/formula.repository
	rm -f $(srcdir)/doc/api/latex/*.tex
	rm -f $(srcdir)/doc/api/latex/*.sty
	rm -f $(srcdir)/doc/api/latex/*.ttf
	rm -f $(srcdir)/doc/api/latex/*.pdf
	rm -f $(srcdir)/doc/api/latex/Makefile

dist-hook:
	if test ! -f "$(srcdir)/doc/api/html/index.html"; then \
		chmod -Rf u+rwX $(srcdir)/doc; \
		sed 's/^ \* \\section \([a-zA-Z0-9_]*\) \(.*\)/<ul><li><a href="#\1">\2<\/a><\/li><\/ul>/' < $(srcdir)/doxygen.dox > $(srcdir)/doc/api/html/toc.txt; \
		sed 's/^ \* \\subsection \([a-zA-Z0-9_]*\) \(.*\)/<ul><ul><li><a href="#\1">\2<\/a><\/li><\/ul><\/ul>/' < $(srcdir)/doxygen.dox > $(srcdir)/doc/api/html/toc.txt; \
		sed 's/^ \* \\subsubsection \([a-zA-Z0-9_]*\) \(.*\)/<ul><ul><ul><li><a href="#\1">\2<\/a><\/li><\/ul><\/ul><\/ul>/' < $(srcdir)/doxygen.dox > $(srcdir)/doc/api/html/toc.txt; \
		doxygen $(builddir)/doxygen/doxygen.Doxyfile; \
		rm -f $(srcdir)/doc/api/html/toc.txt; \
		$(srcdir)/support/rem.sed < $(srcdir)/doc/api/html/index.html > $(srcdir)/doc/api/html/index.html.new; \
		mv $(srcdir)/doc/api/html/index.html.new $(srcdir)/doc/api/html/index.html; \
		cp -R $(srcdir)/doc $(distdir)/; \
		rm -Rf `find $(distdir) -name .svn -type d`; \
	else \
		cp -R $(srcdir)/doc $(distdir)/; \
	fi

doc:

# Flags to be passed to aclocal.
ACLOCAL_AMFLAGS = -I m4
