/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <CUnit/CUnit.h>

#include "nfft3.h"
#include "infft.h"
#include "window.h"

/* The window macros in infft.h fall through to Kaiser-Bessel. */
#if !defined(DIRAC_DELTA) && !defined(GAUSSIAN) && !defined(B_SPLINE) \
    && !defined(SINC_POWER)
#define WINDOW_IS_KAISER_BESSEL 1
#endif

#define ERR(x,y) IF(ABS(y) == K(0.0), ABS((x) - (y)), ABS((x) - (y)) / ABS(y))

/* The test build carries -ffast-math, which rules out isnan and isinf, so bound
 * the magnitude instead. Every peak-normalized window value is at most one. */
#define FINITE(v) (ABS(v) < K(1.0E30))

/* Slack for the shape checks, which compare neighbouring window values rather
 * than measuring accuracy. */
static R bound(void)
{
  return K(64.0) * Y(float_property)(NFFT_EPSILON);
}

/* Accuracy against a reference. Worst measured across the checks that use it is
 * 2.05 ulp, so this bites if anything moves. */
static R tight_bound(void)
{
  return K(16.0) * Y(float_property)(NFFT_EPSILON);
}

/* Arguments straddling NFFT_I0_ASYMP_SPLIT in every precision, up to a peak
 * far beyond what I0 itself can represent. */
static const R i0_args[] =
{
  K(0.0), K(0.5), K(1.0), K(3.0), K(7.0), K(14.5), K(15.0), K(15.5),
  K(24.0), K(25.5), K(40.0), K(80.0), K(200.0)
};

/* The two scaled forms describe the same function, so exp(logtail(x)) is
 * exp_scaled(x). Checking them against each other covers the arguments the
 * reference tables in bessel.c do not reach, in particular a peak far past
 * where I0 itself would overflow. The exponential is taken here rather than
 * inside either helper, so its own error is the test's, not the library's. */
void X(check_bessel_i0_scaling)(void)
{
  const R b = tight_bound();
  unsigned int j;

  printf("BESSEL I0 SCALING\n-----------------\n");

  for (j = 0; j < sizeof(i0_args) / sizeof(i0_args[0]); j++)
  {
    const R x = i0_args[j];
    const R e = Y(bessel_i0_exp_scaled)(x);
    const R lt = Y(bessel_i0_logtail)(x);
    /* EXP(lt) inherits an error that grows with |lt|, which grows with x */
    const R tol = b * (K(1.0) + ABS(lt));
    R err;
    int ok;

    /* exp(-x) I0(x) lies in (0, 1] whatever x, and never forms I0(x) */
    ok = IF(FINITE(e) && e > K(0.0) && e <= K(1.0), 1, 0);
    printf("i0e[" __FE__ "] = " __FE__ " in (0,1] -> %-4s\n", x, e,
      IF(ok == 0, "FAIL", "OK"));
    CU_ASSERT(ok == 1);

    /* log I0(x) - x is at most zero and agrees with the scaled form */
    err = ERR(EXP(lt), e);
    ok = IF(err < tol && FINITE(lt) && lt <= K(0.0), 1, 0);
    printf("i0_logtail[" __FE__ "] = " __FE__ " err_rel = " __FE__ " %-2s "
      __FE__ " -> %-4s\n", x, lt, err, IF(ok == 0, ">=", "<"), tol,
      IF(ok == 0, "FAIL", "OK"));
    CU_ASSERT(ok == 1);
  }
}

/* Oversampling factors and cut-offs. The large cut-offs put I0(m b) far above
 * the overflow threshold of every precision. */
static const R kb_sigma[] = {K(1.25), K(1.5), K(2.0)};
static const int kb_m[] = {2, 4, 6, 8, 12, 16, 20, 24};

void X(check_kaiser_bessel_peak)(void)
{
  const R b = bound();
  const INT N = 64;
  unsigned int s, i;

  printf("KAISER-BESSEL PEAK\n------------------\n");

  for (s = 0; s < sizeof(kb_sigma) / sizeof(kb_sigma[0]); s++)
  {
    for (i = 0; i < sizeof(kb_m) / sizeof(kb_m[0]); i++)
    {
      const int m = kb_m[i];
      const INT n = 2 * (INT)(kb_sigma[s] * (R)N / K(2.0));
      const R sh = KPI * (K(2.0) - (R)N / (R)n);
      const R lt = Y(bessel_i0_logtail)((R)m * sh);
      const R pk = K(1.0) / Y(bessel_i0_exp_scaled)((R)m * sh);
      const R pki = EXP(-(R)m * sh - lt);
      const R peak = Y(kb_phi_hut)(sh, pk, (R)m, (R)n, K(0.0));
      const R err = ERR(peak, K(1.0));
      const int ok = IF(err < b, 1, 0);
      R prev = peak;
      INT k;
      int l;

      printf("phi_hut[sigma=" __FE__ ", m=%2d, k=0] = " __FE__ " err_rel = "
        __FE__ " %-2s " __FE__ " -> %-4s\n", (R)n / (R)N, m, peak, err,
        IF(ok == 0, ">=", "<"), b, IF(ok == 0, "FAIL", "OK"));
      CU_ASSERT(ok == 1);

      /* phi_hut decays from the peak across the band and stays finite */
      for (k = 1; k <= N / 2; k++)
      {
        const R v = Y(kb_phi_hut)(sh, pk, (R)m, (R)n, (R)k);
        const int okk = IF(FINITE(v) && v > K(0.0)
            && v <= prev * (K(1.0) + b), 1, 0);
        if (okk == 0)
          printf("phi_hut[sigma=" __FE__ ", m=%2d, k=" __D__ "] = " __FE__
            " -> FAIL\n", (R)n / (R)N, m, k, v);
        CU_ASSERT(okk == 1);
        prev = v;
      }

      /* phi peaks at the node and decays to zero across its support; at the
       * far edge the peak-normalized value is below the smallest number the
       * precision can hold, so it may be zero */
      prev = Y(kb_phi)(sh, lt, pki, (R)m, K(0.0));
      CU_ASSERT(FINITE(prev) && prev > K(0.0));
      for (l = 1; l <= m; l++)
      {
        const R v = Y(kb_phi)(sh, lt, pki, (R)m, (R)l);
        const int okl = IF(FINITE(v) && v >= K(0.0)
            && v <= prev * (K(1.0) + b), 1, 0);
        if (okl == 0)
          printf("phi[sigma=" __FE__ ", m=%2d, l=%d] = " __FE__ " -> FAIL\n",
            (R)n / (R)N, m, l, v);
        CU_ASSERT(okl == 1);
        prev = v;
      }
    }
  }
}

/* phi_hut for n = 128, N = 64, m = 6, k = 0 .. 32, at 45 digits. */
static const R kb_phi_hut_ref_128_64_6[] =
{
  K(1.00000000000000000000000000000000000000000000),
  K(9.98494489089584995517772679023520262620547445e-1),
  K(9.93991063421376749393128110159622542611031748e-1),
  K(9.86528860853616151095968774036530557865710218e-1),
  K(9.76172503508485311703567660370461667798253571e-1),
  K(9.63011198343454949996506783403547278610540124e-1),
  K(9.47157501729044250478264403746813636845055708e-1),
  K(9.28745771143545606954894413431302977632503544e-1),
  K(9.07930332765614390471081612878464847372408219e-1),
  K(8.84883398751903021414118839751469598193588794e-1),
  K(8.59792772219546791899966821635962597227125582e-1),
  K(8.32859381323620690092249720637173789815349855e-1),
  K(8.04294686262121018881310156140543140233668599e-1),
  K(7.74318004514606635981658095666373156325724635e-1),
  K(7.43153800109706175321610582899518609301241107e-1),
  K(7.11028982230902751137988415595282930898507300e-1),
  K(6.78170257043523882963137139910475559686833000e-1),
  K(6.44801574315932798858343361214745551677727178e-1),
  K(6.11141707292725979593819178171097615746106469e-1),
  K(5.77402000453743929691421216729615246765526617e-1),
  K(5.43784315371465826826544276254058433839502290e-1),
  K(5.10479199984059489936203273658596322747366992e-1),
  K(4.77664301363033583070705281307420845696971326e-1),
  K(4.45503036608129462145906473582831271759099813e-1),
  K(4.14143530982990040955928383684286730919146642e-1),
  K(3.83717826944883261415883852730364716376773438e-1),
  K(3.54341362444930115926805283810669803207466001e-1),
  K(3.26112711896312364923809210421517855897623211e-1),
  K(2.99113578628357982932064370688739584153949427e-1),
  K(2.73409023550682881336433098289080836077217039e-1),
  K(2.49047911213415444713042653987731682982830827e-1),
  K(2.26063551518812896793627739933671092376900165e-1),
  K(2.04474513049807467159431457718724882098004999e-1),
};

void X(check_kaiser_bessel_reference)(void)
{
  const int m = 6;
  const INT N = 64, n = 128;
  const R sh = KPI * (K(2.0) - (R)N / (R)n);
  const R pk = K(1.0) / Y(bessel_i0_exp_scaled)((R)m * sh);
  /* against a table, not against the library's own I0, so this is a real
   * accuracy bound. Worst measured 2.01 ulp. */
  const R b = tight_bound();
  unsigned int k;

  printf("KAISER-BESSEL REFERENCE\n-----------------------\n");

  for (k = 0; k < sizeof(kb_phi_hut_ref_128_64_6) / sizeof(R); k++)
  {
    const R ref = kb_phi_hut_ref_128_64_6[k];
    const R v = Y(kb_phi_hut)(sh, pk, (R)m, (R)n, (R)k);
    const R err = ERR(v, ref);
    const int ok = IF(err < b, 1, 0);
    printf("phi_hut[k=%u] = " __FE__ " err_rel = " __FE__ " %-2s " __FE__
      " -> %-4s\n", k, v, err, IF(ok == 0, ">=", "<"), b,
      IF(ok == 0, "FAIL", "OK"));
    CU_ASSERT(ok == 1);
  }
}

/* phi_hut against reference I0(m*sqrt(b^2-t^2))/I0(m*b) at 60 digits, for the near-peak
 * bins where the two Bessel arguments nearly coincide, plus the band edge.
 *
 * Regenerate with b = pi*(2 - N/n), t = 2*pi*k/n and
 *   mpmath.besseli(0, m*mpmath.sqrt(b*b - t*t)) / mpmath.besseli(0, m*b)
 * at mpmath.mp.dps = 60. */
static const R kb_phi_hut_ref_512_256_8[] =
{
  K(1.0),
  K(9.99873882985617506643870290279057361423044e-1),
  K(9.99495624837643910978594435458575184590944e-1),
  K(9.98865504130588397987934385630536020655212e-1),
  K(9.97983984785438854687931348979505266186542e-1)
};

/* Value at the band edge, same parameters. */
#define KB_PHI_HUT_EDGE_512_256_8 K(1.19266391939960090608558187858420025398744e-1)

void X(check_kaiser_bessel_cancellation)(void)
{
  const INT N = 256, n = 512;
  const int m = 8;
  const R sh = KPI * (K(2.0) - (R)N / (R)n);
  const R lt = Y(bessel_i0_logtail)((R)m * sh);
  const R pk = K(1.0) / Y(bessel_i0_exp_scaled)((R)m * sh);
  /* leaves room for the rounding of b and t, which no formulation can undo */
  const R b = tight_bound();
  unsigned int k;

  printf("KAISER-BESSEL CANCELLATION\n--------------------------\n");

  for (k = 0; k < sizeof(kb_phi_hut_ref_512_256_8) / sizeof(R); k++)
  {
    const R ref = kb_phi_hut_ref_512_256_8[k];
    const R v = Y(kb_phi_hut)(sh, pk, (R)m, (R)n, (R)k);
    const R err = ERR(v, ref);
    const int ok = IF(err < b, 1, 0);
    printf("phi_hut[k=%u] = " __FE__ " err_rel = " __FE__ " %-2s " __FE__
      " -> %-4s\n", k, v, err, IF(ok == 0, ">=", "<"), b,
      IF(ok == 0, "FAIL", "OK"));
    CU_ASSERT(ok == 1);
  }

  {
    const R ref = KB_PHI_HUT_EDGE_512_256_8;
    const R v = Y(kb_phi_hut)(sh, pk, (R)m, (R)n, (R)(N / 2));
    const R err = ERR(v, ref);
    const int ok = IF(err < b, 1, 0);
    printf("phi_hut[k=" __D__ "] = " __FE__ " err_rel = " __FE__ " %-2s " __FE__
      " -> %-4s\n", N / 2, v, err, IF(ok == 0, ">=", "<"), b,
      IF(ok == 0, "FAIL", "OK"));
    CU_ASSERT(ok == 1);
  }
}

void X(check_kaiser_bessel_nfft)(void)
{
#if defined(WINDOW_IS_KAISER_BESSEL)
  /* Cut-offs for which the unnormalized I0(m b) raised to the power d leaves
   * the range of single precision. */
  static const int cut_off[] = {6, 8, 10, 12, 14};
  const int d = 3, M = 24;
  int N[3] = {16, 16, 16}, n[3] = {32, 32, 32};
  unsigned int i;

  printf("KAISER-BESSEL NFFT\n------------------\n");

  for (i = 0; i < sizeof(cut_off) / sizeof(cut_off[0]); i++)
  {
    const R sigma = (R)n[0] / (R)N[0];
    /* roundoff plus the truncation error of the window itself */
    const R b = K(1000.0) * Y(float_property)(NFFT_EPSILON)
        + K(4.0) * EXP(K(-2.0) * KPI * (R)cut_off[i]
            * SQRT(K(1.0) - K(1.0) / sigma));
    Y(plan) p;
    C *ref;
    R err;
    int j, ok;

    Y(init_guru)(&p, d, N, M, n, cut_off[i],
        PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT,
        FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

    Y(vrand_shifted_unit_double)(p.x, p.d * p.M_total);
    Y(vrand_unit_complex)((C*)p.f_hat, p.N_total);
    Y(precompute_one_psi)(&p);

    ref = (C*)Y(malloc)((size_t)p.M_total * sizeof(C));
    Y(trafo_direct)(&p);
    for (j = 0; j < p.M_total; j++)
      ref[j] = ((C*)p.f)[j];

    Y(trafo)(&p);
    err = Y(error_l_infty_1_complex)(ref, (C*)p.f, p.M_total, (C*)p.f_hat,
        p.N_total);
    ok = IF(FINITE(err) && err < b, 1, 0);
    printf("nfft_3d[m=%2d] err = " __FE__ " %-2s " __FE__ " -> %-4s\n",
      cut_off[i], err, IF(ok == 0, ">=", "<"), b, IF(ok == 0, "FAIL", "OK"));
    CU_ASSERT(ok == 1);

    Y(free)(ref);
    Y(finalize)(&p);
  }
#else
  printf("KAISER-BESSEL NFFT\n------------------\nskipped: %s window\n",
    Y(get_window_name)());
#endif
}

/* kb_phi_run splits a run into a guarded stretch and a branch-free interior.
 * The two must agree, at every offset of the run against the support and at
 * every cutoff, including the offsets that put ra at or next to zero. */
void X(check_kaiser_bessel_run)(void)
{
#if defined(WINDOW_IS_KAISER_BESSEL)
  static const int kb_run_m[] = {1, 2, 4, 8, 11, 16};
  /* fractions of a grid cell, the ends included: those are where ra vanishes */
  static const R kb_run_frac[] =
  {
    K(0.0), K(1.0E-12), K(1.0E-6), K(0.25), K(0.5), K(0.75),
    K(1.0) - K(1.0E-6), K(1.0)
  };
  const R b = tight_bound();
  const R sh = KPI * K(1.5); /* sigma = 2 */
  unsigned int i, f;

  printf("KAISER-BESSEL RUN\n-----------------\n");

  for (i = 0; i < sizeof(kb_run_m) / sizeof(kb_run_m[0]); i++)
  {
    const INT m = (INT)kb_run_m[i];
    const R lt = Y(bessel_i0_logtail)((R)m * sh);
    const R pki = EXP(-(R)m * sh - lt);
    const R peak = Y(kb_phi)(sh, lt, pki, (R)m, K(0.0));
    R *run = (R*)Y(malloc)((size_t)(2 * m + 2) * sizeof(R));

    for (f = 0; f < sizeof(kb_run_frac) / sizeof(kb_run_frac[0]); f++)
    {
      const R nx0 = (R)m + kb_run_frac[f];
      R err = K(0.0);
      INT l, worst = 0;
      int ok;

      Y(kb_phi_run)(run, sh, lt, pki, (R)m, m, nx0);

      /* Scaled by the peak of the run, not by each value: down in the tail the
       * exponent carries m b, so two roundings of it differ by more than an ulp
       * of a value that is itself 1E-15 of the peak. What the psi table owes
       * the transform is accuracy against the peak. */
      for (l = 0; l < 2 * m + 2; l++)
      {
        const R ref = Y(kb_phi)(sh, lt, pki, (R)m, nx0 - (R)l);
        R e;

        if (!(FINITE(run[l])))
        {
          err = K(1.0);
          worst = l;
          break;
        }
        e = ABS(run[l] - ref) / peak;
        if (e > err)
        {
          err = e;
          worst = l;
        }
      }

      ok = IF(err < b, 1, 0);
      printf("run[m=" __D__ ", frac=" __FE__ "] max err_peak = " __FE__ " %-2s "
        __FE__ " at l = " __D__ " -> %-4s\n", m, kb_run_frac[f], err,
        IF(ok == 0, ">=", "<"), b, worst, IF(ok == 0, "FAIL", "OK"));
      CU_ASSERT(ok == 1);
    }

    Y(free)(run);
  }
#else
  printf("KAISER-BESSEL RUN\n-----------------\nskipped: %s window\n",
    Y(get_window_name)());
#endif
}


/* kb_phi over one full run for n = 512, m = 8, sigma = 2, with the node three
 * tenths of a cell past a grid point, at 45 digits. The run straddles the
 * support, so it covers the interior kernel, both guarded ends and the
 * out-of-support tail in one sweep. Scaled by the window peak, which is what
 * the psi table owes the transform: relative error is not the yardstick out in
 * a tail that is 1E-19 of the peak. */
static const R kb_phi_ref_512_8_fr3[] =
{
  K(-7.85892543958984626234760839987669431974422115e-17),
  K(1.57695746334769672676281540979331619263915038e-10),
  K(2.58719341293312107085540589899031766706724005e-7),
  K(3.17336750299790857820829439712639005468306356e-5),
  K(9.82642935193592753169649014657611814547481724e-4),
  K(1.16729248907278442847979688855228428645887835e-2),
  K(6.48668395087801286575741926583042819130430870e-2),
  K(1.87379688129696004223212763593580831265537010e-1),
  K(2.97382623930255136946603580473539321010687990e-1),
  K(2.65093116557523689981137022210302751246001787e-1),
  K(1.32016745641728771175013488455801323624395690e-1),
  K(3.54925666531525603998211931133479680415846357e-2),
  K(4.79138596108146813158167803477516328795267173e-3),
  K(2.83776673136889391974806346382884468525264064e-4),
  K(5.67562308714283037084637665709343759678174769e-6),
  K(2.09422510977472083746492640250005379927681671e-8),
  K(1.31872120324319329691486991491412240554570553e-12),
  K(-2.38079531924638347327199506798677762538172144e-17),
};

#define KB_PHI_PEAK_512_8 K(3.05158861803586006004326791417863091093860567e-1)

void X(check_kaiser_bessel_phi)(void)
{
#if defined(WINDOW_IS_KAISER_BESSEL)
  const INT m = 8, n = 512;
  const R sh = KPI * K(1.5); /* sigma = 2 */
  const R lt = Y(bessel_i0_logtail)((R)m * sh);
  const R pki = EXP(-(R)m * sh - lt);
  const R peak = KB_PHI_PEAK_512_8;
  /* worst measured 2.05 ulp of the peak */
  const R b = tight_bound();
  const R nx0 = (R)m + K(0.3);
  unsigned int l;

  printf("KAISER-BESSEL PHI\n-----------------\n");

  for (l = 0; l < sizeof(kb_phi_ref_512_8_fr3) / sizeof(R); l++)
  {
    const R ref = kb_phi_ref_512_8_fr3[l];
    const R v = Y(kb_phi)(sh, lt, pki, (R)m, nx0 - (R)l);
    const R err = ABS(v - ref) / peak;
    const int ok = IF(FINITE(v) && err < b, 1, 0);
    printf("phi[l=%u] = " __FE__ " err_peak = " __FE__ " %-2s " __FE__
      " -> %-4s\n", l, v, err, IF(ok == 0, ">=", "<"), b,
      IF(ok == 0, "FAIL", "OK"));
    CU_ASSERT(ok == 1);
  }
  (void)n;
#else
  printf("KAISER-BESSEL PHI\n-----------------\nskipped: %s window\n",
    Y(get_window_name)());
#endif
}
