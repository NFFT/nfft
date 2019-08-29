LIBTOOL_DEPS = @LIBTOOL_DEPS@

libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck

# Subdirectories
DIST_SUBDIRS=include kernel . tests examples applications matlab julia support doxygen

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

if HAVE_JULIA
  JULIA_DIRS=julia
  LIBNFFT3_JULIA_LA=libnfft3@PREC_SUFFIX@_julia.la
else
  JULIA_DIRS=
  LIBNFFT3_JULIA_LA=
endif

if HAVE_THREADS
  LIBNFFT3_THREADS_LA = libnfft3@PREC_SUFFIX@_threads.la
else
  LIBNFFT3_THREADS_LA =
endif

SUBDIRS= include kernel . tests $(EXAMPLE_DIRS) $(APPLICATION_DIRS) $(MATLAB_DIRS) $(JULIA_DIRS)

lib_LTLIBRARIES = libnfft3@PREC_SUFFIX@.la $(LIBNFFT3_THREADS_LA)
noinst_LTLIBRARIES = $(LIBNFFT3_MATLAB_LA) $(LIBNFFT3_JULIA_LA)

libnfft3@PREC_SUFFIX@_la_SOURCES =
libnfft3@PREC_SUFFIX@_la_LIBADD = kernel/libkernel.la @fftw3_LIBS@ -lm
libnfft3@PREC_SUFFIX@_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@ @fftw3_LDFLAGS@
# Removed since it is only used by iterS2:  3rdparty/lib3rdparty.la

if HAVE_THREADS
libnfft3@PREC_SUFFIX@_threads_la_SOURCES =
libnfft3@PREC_SUFFIX@_threads_la_LIBADD = kernel/libkernel_threads.la @fftw3_threads_LIBS@ -lm
libnfft3@PREC_SUFFIX@_threads_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@ @fftw3_LDFLAGS@
if HAVE_OPENMP
  libnfft3@PREC_SUFFIX@_threads_la_CFLAGS = $(OPENMP_CFLAGS)
endif
endif

if HAVE_MATLAB
  libnfft3@PREC_SUFFIX@_matlab_la_SOURCES =
if HAVE_MATLAB_THREADS
    libnfft3@PREC_SUFFIX@_matlab_la_LIBADD = kernel/libkernel_threads.la @matlab_fftw3_LIBS@ -lm
else
    libnfft3@PREC_SUFFIX@_matlab_la_LIBADD = kernel/libkernel.la @matlab_fftw3_LIBS@ -lm
endif
  libnfft3@PREC_SUFFIX@_matlab_la_LDFLAGS = @matlab_fftw3_LDFLAGS@
endif

if HAVE_JULIA
  libnfft3@PREC_SUFFIX@_julia_la_SOURCES =
if HAVE_THREADS
    libnfft3@PREC_SUFFIX@_julia_la_LIBADD = kernel/libkernel_threads.la @fftw3_threads_LIBS@ -lm
else
    libnfft3@PREC_SUFFIX@_julia_la_LIBADD = kernel/libkernel.la @fftw3_LIBS@ -lm
endif
  libnfft3@PREC_SUFFIX@_julia_la_LDFLAGS = @fftw3_LDFLAGS@
endif

EXTRA_DIST = bootstrap.sh doxygen.dox nfft3.pc.in doc macos-build-mex.sh linux-build-mex.sh windows-build-dll.sh

if HAVE_NON_DOUBLE_PRECISION
nfft3@PREC_SUFFIX@.pc: nfft3.pc
	cp -f nfft3.pc nfft3@PREC_SUFFIX@.pc
endif

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = nfft3@PREC_SUFFIX@.pc

install-data-hook:
	if test -f "$(abs_top_builddir)/doc/html/index.html"; then \
		$(MKDIR_P) $(DESTDIR)$(docdir); \
		cp -Rf $(abs_top_builddir)/doc/html $(DESTDIR)$(docdir)/; \
	elif test -f "$(abs_top_srcdir)/doc/html/index.html"; then \
		$(MKDIR_P) $(DESTDIR)$(docdir); \
		cp -Rf $(abs_top_srcdir)/doc/html $(DESTDIR)$(docdir)/; \
	fi

uninstall-hook:
	chmod -Rf u+rwX $(DESTDIR)$(docdir)
	rm -Rf $(DESTDIR)$(docdir)

clean-local:
	rm -f $(abs_top_builddir)/doc/doxygen_sqlite3.db
	rm -f $(abs_top_builddir)/doc/nfft.pdf
	rm -f $(abs_top_builddir)/doc/nfft.tag
	rm -Rf $(abs_top_builddir)/doc/html
	rm -Rf $(abs_top_builddir)/doc/latex

dist-hook: doc
	rm -f @DX_DOCDIR@/@PACKAGE@.tag
	cp -R $(abs_top_builddir)/doc $(distdir)/
	rm -Rf `find $(distdir) -name .svn -type d`

include aminclude.am

if DX_COND_doc
doc: doxygen-doc
	rm -f $(abs_top_builddir)/doc/html/*.md5
	rm -f $(abs_top_builddir)/doc/html/*.map
endif

# Flags to be passed to aclocal.
ACLOCAL_AMFLAGS = -I m4

AM_DISTCHECK_CONFIGURE_FLAGS = --enable-all
