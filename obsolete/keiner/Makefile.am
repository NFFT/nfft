lib_LTLIBRARIES = libnfsft.la
      
SUBDIRS = kernel . examples

EXTRA_DIST = nfsft.lsm.in nfsft.spec.in main.1 manpage.tpl

libnfsft_la_SOURCES =

libnfsft_la_LIBADD = kernel/libkernel.la

libnfsft_la_LDFLAGS = -version-info 3:0:0
