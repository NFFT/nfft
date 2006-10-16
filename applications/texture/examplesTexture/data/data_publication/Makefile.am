AM_CFLAGS = -Wall
AM_CPPFLAGS = -I../util -I$(top_srcdir)/include -I../../.. \
							-I../../texture_util

noinst_PROGRAMS = create

create_SOURCES = create.c

create_LDADD = ../../texture_util/libtexture_util.la ../util/libdata_util.la
