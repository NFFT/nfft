AM_CFLAGS = -I$(top_srcdir)/include

noinst_LTLIBRARIES = libsolver.la

libsolver_la_SOURCES = solver.c solver2.c
if MAINTAINER_MODE
  nodist_libsolver_la_SOURCES = solver_adjoint.c
endif

