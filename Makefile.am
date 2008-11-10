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

SUBDIRS= 3rdparty include util kernel . examples applications $(EXTRA_DIRS) $(MATLAB_DIRS)

lib_LTLIBRARIES = libnfft3.la

libnfft3_la_SOURCES =
libnfft3_la_LIBADD = 3rdparty/lib3rdparty.la util/libutil.la kernel/libkernel.la -lm
libnfft3_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@

EXTRA_DIST = bootstrap.sh VERSION doxygen.c

install-data-hook:
	mkdir -p $(pkgdatadir)
	cp -R doc $(pkgdatadir)

clean-local:
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
	doxygen doxygen/doxygen.Doxyfile
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
