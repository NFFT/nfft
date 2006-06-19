/* specfunc/bessel.c
 * 
 * Copyright (C) 1996,1997,1998,1999,2000,2001,2002,2003 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */
/* Miscellaneous support functions for Bessel function evaluations.
 */
#include "config.h"
#include "../gsl_math.h"
#include "../err/gsl_errno.h"
#include "gsl_sf_elementary.h"
#include "gsl_sf_exp.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_trig.h"

#include "error.h"

#include "bessel_amp_phase.h"
#include "bessel.h"

#define CubeRoot2_  1.25992104989487316476721060728



/* Debye functions [Abramowitz+Stegun, 9.3.9-10] */

inline static double 
debye_u1(const double * tpow)
{
  return (3.0*tpow[1] - 5.0*tpow[3])/24.0;
}

inline static double 
debye_u2(const double * tpow)
{
  return (81.0*tpow[2] - 462.0*tpow[4] + 385.0*tpow[6])/1152.0;
}

inline
static double debye_u3(const double * tpow)
{
  return (30375.0*tpow[3] - 369603.0*tpow[5] + 765765.0*tpow[7] - 425425.0*tpow[9])/414720.0;
}

inline
static double debye_u4(const double * tpow)
{
  return (4465125.0*tpow[4] - 94121676.0*tpow[6] + 349922430.0*tpow[8] - 
          446185740.0*tpow[10] + 185910725.0*tpow[12])/39813120.0;
}

inline
static double debye_u5(const double * tpow)
{
  return (1519035525.0*tpow[5]     - 49286948607.0*tpow[7] + 
          284499769554.0*tpow[9]   - 614135872350.0*tpow[11] + 
          566098157625.0*tpow[13]  - 188699385875.0*tpow[15])/6688604160.0;
}

#if 0
inline
static double debye_u6(const double * tpow)
{
  return (2757049477875.0*tpow[6] - 127577298354750.0*tpow[8] + 
          1050760774457901.0*tpow[10] - 3369032068261860.0*tpow[12] + 
          5104696716244125.0*tpow[14] - 3685299006138750.0*tpow[16] + 
          1023694168371875.0*tpow[18])/4815794995200.0;
}
#endif


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_bessel_IJ_taylor_e(const double nu, const double x,
                             const int sign,
                             const int kmax,
                             const double threshold,
                             gsl_sf_result * result
                             )
{
  /* CHECK_POINTER(result) */

  if(nu < 0.0 || x < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x == 0.0) {
    if(nu == 0.0) {
      result->val = 1.0;
      result->err = 0.0;
    }
    else {
      result->val = 0.0;
      result->err = 0.0;
    }
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result prefactor;   /* (x/2)^nu / Gamma(nu+1) */
    gsl_sf_result sum;

    int stat_pre;
    int stat_sum;
    int stat_mul;

    if(nu == 0.0) {
      prefactor.val = 1.0;
      prefactor.err = 0.0;
      stat_pre = GSL_SUCCESS;
    }
    else if(nu < INT_MAX-1) {
      /* Separate the integer part and use
       * y^nu / Gamma(nu+1) = y^N /N! y^f / (N+1)_f,
       * to control the error.
       */
      const int    N = (int)floor(nu + 0.5);
      const double f = nu - N;
      gsl_sf_result poch_factor;
      gsl_sf_result tc_factor;
      const int stat_poch = gsl_sf_poch_e(N+1.0, f, &poch_factor);
      const int stat_tc   = gsl_sf_taylorcoeff_e(N, 0.5*x, &tc_factor);
      const double p = pow(0.5*x,f);
      prefactor.val  = tc_factor.val * p / poch_factor.val;
      prefactor.err  = tc_factor.err * p / poch_factor.val;
      prefactor.err += fabs(prefactor.val) / poch_factor.val * poch_factor.err;
      prefactor.err += 2.0 * GSL_DBL_EPSILON * fabs(prefactor.val);
      stat_pre = GSL_ERROR_SELECT_2(stat_tc, stat_poch);
    }
    else {
      gsl_sf_result lg;
      const int stat_lg = gsl_sf_lngamma_e(nu+1.0, &lg);
      const double term1  = nu*log(0.5*x);
      const double term2  = lg.val;
      const double ln_pre = term1 - term2;
      const double ln_pre_err = GSL_DBL_EPSILON * (fabs(term1)+fabs(term2)) + lg.err;
      const int stat_ex = gsl_sf_exp_err_e(ln_pre, ln_pre_err, &prefactor);
      stat_pre = GSL_ERROR_SELECT_2(stat_ex, stat_lg);
    }

    /* Evaluate the sum.
     * [Abramowitz+Stegun, 9.1.10]
     * [Abramowitz+Stegun, 9.6.7]
     */
    {
      const double y = sign * 0.25 * x*x;
      double sumk = 1.0;
      double term = 1.0;
      int k;

      for(k=1; k<=kmax; k++) {
        term *= y/((nu+k)*k);
        sumk += term;
        if(fabs(term/sumk) < threshold) break;
      }

      sum.val = sumk;
      sum.err = threshold * fabs(sumk);

      stat_sum = ( k >= kmax ? GSL_EMAXITER : GSL_SUCCESS );
    }

    stat_mul = gsl_sf_multiply_err_e(prefactor.val, prefactor.err,
                                        sum.val, sum.err,
                                        result);

    return GSL_ERROR_SELECT_3(stat_mul, stat_pre, stat_sum);
  }
}

/* nu -> Inf; uniform in x > 0  [Abramowitz+Stegun, 9.7.7]
 *
 * error:
 *   The error has the form u_N(t)/nu^N  where  0 <= t <= 1.
 *   It is not hard to show that |u_N(t)| is small for such t.
 *   We have N=6 here, and |u_6(t)| < 0.025, so the error is clearly
 *   bounded by 0.025/nu^6. This gives the asymptotic bound on nu
 *   seen below as nu ~ 100. For general MACH_EPS it will be 
 *                     nu > 0.5 / MACH_EPS^(1/6)
 *   When t is small, the bound is even better because |u_N(t)| vanishes
 *   as t->0. In fact u_N(t) ~ C t^N as t->0, with C ~= 0.1.
 *   We write
 *                     err_N <= min(0.025, C(1/(1+(x/nu)^2))^3) / nu^6
 *   therefore
 *                     min(0.29/nu^2, 0.5/(nu^2+x^2)) < MACH_EPS^{1/3}
 *   and this is the general form.
 *
 * empirical error analysis, assuming 14 digit requirement:
 *   choose   x > 50.000 nu   ==>  nu >   3
 *   choose   x > 10.000 nu   ==>  nu >  15
 *   choose   x >  2.000 nu   ==>  nu >  50
 *   choose   x >  1.000 nu   ==>  nu >  75
 *   choose   x >  0.500 nu   ==>  nu >  80
 *   choose   x >  0.100 nu   ==>  nu >  83
 *
 * This makes sense. For x << nu, the error will be of the form u_N(1)/nu^N,
 * since the polynomial term will be evaluated near t=1, so the bound
 * on nu will become constant for small x. Furthermore, increasing x with
 * nu fixed will decrease the error.
 */
int
gsl_sf_bessel_Inu_scaled_asymp_unif_e(const double nu, const double x, gsl_sf_result * result)
{
  int i;
  double z = x/nu;
  double root_term = hypot(1.0,z);
  double pre = 1.0/sqrt(2.0*M_PI*nu * root_term);
  double eta = root_term + log(z/(1.0+root_term));
  double ex_arg = ( z < 1.0/GSL_ROOT3_DBL_EPSILON ? nu*(-z + eta) : -0.5*nu/z*(1.0 - 1.0/(12.0*z*z)) );
  gsl_sf_result ex_result;
  int stat_ex = gsl_sf_exp_e(ex_arg, &ex_result);
  if(stat_ex == GSL_SUCCESS) {
    double t = 1.0/root_term;
    double sum;
    double tpow[16];
    tpow[0] = 1.0;
    for(i=1; i<16; i++) tpow[i] = t * tpow[i-1];
    sum = 1.0 + debye_u1(tpow)/nu + debye_u2(tpow)/(nu*nu) + debye_u3(tpow)/(nu*nu*nu)
          + debye_u4(tpow)/(nu*nu*nu*nu) + debye_u5(tpow)/(nu*nu*nu*nu*nu);
    result->val  = pre * ex_result.val * sum;
    result->err  = pre * ex_result.val / (nu*nu*nu*nu*nu*nu);
    result->err += pre * ex_result.err * fabs(sum);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = 0.0;
    result->err = 0.0;
    return stat_ex;
  }
}
