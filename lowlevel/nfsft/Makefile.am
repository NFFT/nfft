AM_CFLAGS = -I$(top_srcdir)/include

SUBDIRS = obsolete .

noinst_LTLIBRARIES = libnfsft.la

libnfsft_la_SOURCES = nfsft.c legendre.h legendre.c api.h

libnfsft_la_LIBADD = obsolete/libnfsftold.la
