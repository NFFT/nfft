# Subdirectories
if WITH_TESTS
  SUBDIRS =  3rdparty include util kernel lowlevel highlevel . examples tests
#  DIST_SUBDIRS = include util kernel lowlevel highlevel . examples
else
  SUBDIRS = 3rdparty include util kernel lowlevel highlevel . examples
endif  

lib_LTLIBRARIES = libnfft3.la      

libnfft3_la_SOURCES =

libnfft3_la_LIBADD = 3rdparty/lib3rdparty.la util/libutil.la kernel/libkernel.la lowlevel/liblowlevel.la highlevel/libhighlevel.la

dist-hook:
	doxygen ./doxygen.Doxyfile
	mkdir $(distdir)/doc
	cp -R doc/html $(distdir)/doc/
