# Subdirectories
if WITH_TESTS
  SUBDIRS = include util kernel lowlevel highlevel . examples tests
else
  SUBDIRS = include util kernel lowlevel highlevel . examples
endif  

lib_LTLIBRARIES = libnfft3.la      

libnfft3_la_SOURCES =

libnfft3_la_LIBADD = util/libutil.la kernel/libkernel.la lowlevel/liblowlevel.la highlevel/libhighlevel.la

