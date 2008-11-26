AM_CPPFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

TESTS = test

## The program name
nodist_PROGRAMS = test
nodistdir = $(PWD) 

## The source files.
test_SOURCES = test.c

test_LDADD = $(top_builddir)/libnfft3.la
