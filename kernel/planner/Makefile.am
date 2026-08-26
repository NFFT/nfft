AM_CPPFLAGS = -I$(top_srcdir)/include

if ENABLE_OPENMP
  LIBPLANNER_THREADS_LA = libplanner_threads.la
else
  LIBPLANNER_THREADS_LA =
endif

noinst_LTLIBRARIES = libplanner.la $(LIBPLANNER_THREADS_LA)

libplanner_la_SOURCES = md5.c print.c printers.c scan.c scanners.c \
  primes.c solver.c solvtab.c planner.c tensor.c problem.c plan.c \
  theplanner.c timer.c

if ENABLE_OPENMP
libplanner_threads_la_SOURCES = $(libplanner_la_SOURCES)
libplanner_threads_la_CFLAGS = $(OPENMP_CFLAGS)
endif
