SUBDIRS = sys err ieee-utils utils specfunc .

noinst_LTLIBRARIES = libgsl.la

libgsl_la_SOURCES =

libgsl_la_LIBADD = err/libgslerr.la ieee-utils/libgslieeeutils.la \
  specfunc/libgslspecfunc.la sys/libgslsys.la utils/libutils.la

noinst_HEADERS = gsl_machine.h gsl_math.h gsl_mode.h gsl_nan.h gsl_pow_int.h\
  gsl_precision.h gsl_types.h
