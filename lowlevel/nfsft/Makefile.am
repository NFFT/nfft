AM_CFLAGS = -I$(top_srcdir)/include

SUBDIRS = obsolete .

noinst_LTLIBRARIES = libnfsft.la

libnfsft_la_SOURCES = api.h nfsft.c legendre.h legendre.c
#infsft.h infsft.c util.h util.c flft.h flft.c c2f.h c2f.c  u.h u.c direct.h direct.c
#libnfsft_la_SOURCES = 

libnfsft_la_LIBADD = obsolete/libnfsftold.la
