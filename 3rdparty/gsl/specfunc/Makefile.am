noinst_LTLIBRARIES = libgslspecfunc.la 

noinst_HEADERS = gsl_sf.h gsl_sf_bessel.h gsl_sf_debye.h gsl_sf_dilog.h gsl_sf_elementary.h gsl_sf_exp.h gsl_sf_gamma.h gsl_sf_log.h gsl_sf_pow_int.h gsl_sf_psi.h gsl_sf_result.h  gsl_sf_trig.h gsl_sf_zeta.h gsl_specfunc.h bessel_amp_phase.h eval.h chebyshev.h cheb_eval.c cheb_eval_mode.c check.h error.h

libgslspecfunc_la_SOURCES = bessel.c bessel.h bessel_In.c bessel_amp_phase.c bessel_amp_phase.h bessel_i.c bessel_sequence.c debye.c elementary.c exp.c gamma.c log.c poch.c sf_pow_int.c psi.c recurse.h result.c trig.c zeta.c

# TESTS = $(check_PROGRAMS)

# check_PROGRAMS = test

# test_LDADD = libgslspecfunc.la ../ieee-utils/libgslieeeutils.la  ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

# test_SOURCES = test_sf.c test_sf.h test_airy.c test_bessel.c test_coulomb.c test_dilog.c test_gamma.c test_hyperg.c test_legendre.c

AM_CFLAGS = -I$(top_srcdir)/include
