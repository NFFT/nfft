AM_CFLAGS = -I$(top_srcdir)/kernel
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = accuracy

## The source files.
accuracy_SOURCES = accuracy.c

#  add " -lm" to the following if sfft requires the math library
accuracy_LDADD = $(top_builddir)/libnfsft.la