# Subdirectories
SUBDIRS = 3rdparty include util kernel . examples applications

# Development code extra directories.
if WITH_DEVCODE
  nodist_SUBDIRS = tests
endif  

lib_LTLIBRARIES = libnfft3.la      

libnfft3_la_SOURCES =

libnfft3_la_LIBADD = 3rdparty/lib3rdparty.la util/libutil.la kernel/libkernel.la

EXTRA_DIST = bootstrap.sh VERSION

install-data-hook:
	mkdir -p $(pkgdatadir)
	cp -R doc $(pkgdatadir)

dist-hook:
	doxygen doxygen/doxygen.Doxyfile
	mkdir $(distdir)/doc
	mkdir $(distdir)/doc/api
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
	