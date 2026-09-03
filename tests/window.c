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

/* log|sinc(x)| at 60 digits. The arguments are dyadic, so their literals are
 * exact in every precision. Regenerate with mpmath.log(abs(mpmath.sin(x)/x))
 * at mpmath.mp.dps = 60. */
static const R log_sinc_args[] =
{
  K(0.0), K(0.0009765625), K(0.0625), K(0.25), K(0.5), K(0.75), K(1.0),
  K(1.5), K(1.75), K(1.9375), K(2.0), K(2.125), K(2.5), K(3.0), K(3.5),
  K(4.0), K(5.0), K(6.0)
};

static const R log_sinc_ref[] =
{
  K(0.0),
  K(-1.589457244537903157973261989827636306810e-7),
  K(-6.511264587477414330056542128013714207687e-4),
  K(-1.043845457789935401184945573043609341100e-2),
  K(-4.201950582536896172579838403790203712454e-2),
  K(-9.557336640065549370377720803908670442410e-2),
  K(-1.726037462690916785134109758639090698401e-1),
  K(-4.079732642997607841582797016799435332387e-1),
  K(-5.757594515994677430825497749105942216683e-1),
  K(-7.301975001705864507030295078523610852806e-1),
  K(-7.882302166551059406604484922728009331910e-1),
  K(-9.159145793125751502284428482141424133358e-1),
  K(-1.429666029791948529048265147097077863789),
  K(-3.056756918278195617196798690265269719773),
  K(-2.300349799725624325702134320468842225703),
  K(-1.664947325187014385896374493824454668416),
  K(-1.651381082462686277424410260467381509305),
  K(-3.066814833351888510250775163544230639745)
};

void X(check_log_sinc)(void)
{
  const R b = bound();
  unsigned int j;

  printf("LOG SINC\n--------\n");

  for (j = 0; j < SIZE(log_sinc_args); j++)
  {
    const R x = log_sinc_args[j];
    const R v = Y(log_sinc)(x);
    const R err = ERR(v, log_sinc_ref[j]);
    /* Even in x, so the mirrored argument must give the identical bits. */
    const int ok = IF(FINITE(v) && err < b && v == Y(log_sinc)(-x), 1, 0);
    printf("log_sinc(" __FE__ ") = " __FE__ " err_rel = " __FE__ " %-2s " __FE__
      " -> %-4s\n", x, v, err, IF(ok == 0, ">=", "<"), b,
      IF(ok == 0, "FAIL", "OK"));
    CU_ASSERT(ok == 1);
  }
}

/* EXP(log_sinc(x)) must return |sin(x)/x|, across the split and past the first
 * zero of sinc. */
void X(check_log_sinc_exp)(void)
{
  const R b = bound();
  R worst = K(0.0);
  int j, ok;

  printf("LOG SINC ROUND TRIP\n-------------------\n");

  for (j = 1; j <= 60; j++)
  {
    const R x = (R)j / K(20.0);
    const R ref = FABS(SIN(x) / x);
    const R err = ERR(EXP(Y(log_sinc)(x)), ref);
    worst = FMAX(worst, err);
  }

  ok = IF(FINITE(worst) && worst < b, 1, 0);
  printf("exp(log_sinc) err_rel = " __FE__ " %-2s " __FE__ " -> %-4s\n",
    worst, IF(ok == 0, ">=", "<"), b, IF(ok == 0, "FAIL", "OK"));
  CU_ASSERT(ok == 1);
}

#define BSPLINE_RUN_MAX_M 16

/* Every independent Clenshaw evaluation in the run must reproduce what one
 * evaluation per point gives, at every offset of the run against the support
 * and at every cutoff, and the run must sum to one, because the cardinal
 * B-splines of one order partition unity. */
void X(check_bspline_run)(void)
{
#if defined(B_SPLINE)
  static const int bs_run_m[] = {1, 2, 4, 8, 11, BSPLINE_RUN_MAX_M};
  /* fractions of a grid cell, the ends included */
  static const R bs_run_frac[] =
  {
    K(0.0), K(1.0E-12), K(1.0E-6), K(0.25), K(0.5), K(0.75),
    K(1.0) - K(1.0E-6), K(1.0)
  };
  const R b = tight_bound();
  R inv[2 * BSPLINE_RUN_MAX_M];
  unsigned int i, f;
  INT j;

  printf("B-SPLINE RUN\n------------\n");

  for (j = 0; j < 2 * BSPLINE_RUN_MAX_M; j++)
    inv[j] = K(1.0) / (R)(j + 1);

  for (i = 0; i < SIZE(bs_run_m); i++)
  {
    const INT m = (INT)bs_run_m[i];
    R *run = (R*)Y(malloc)((size_t)(2 * m + 2) * sizeof(R));

    for (f = 0; f < SIZE(bs_run_frac); f++)
    {
      const R t = (R)(2 * m) + bs_run_frac[f];
      R peak = K(0.0), err = K(0.0), sum = K(0.0);
      INT l, worst = 0;
      int ok;

      Y(bspline_phi_run)(run, inv, m, t, K(1.0));

      for (l = 0; l < 2 * m + 2; l++)
        peak = FMAX(peak, Y(bsplines)(2 * m, t - (R)l));

      for (l = 0; l < 2 * m + 2; l++)
      {
        const R ref = Y(bsplines)(2 * m, t - (R)l);
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
        sum += run[l];
      }

      ok = IF(err < b && ERR(sum, K(1.0)) < b, 1, 0);
      printf("bspline_run[m=" __D__ ", frac=" __FE__ "] max err_peak = " __FE__
        " %-2s " __FE__ " at l = " __D__ " sum-1 = " __FE__ " -> %-4s\n", m,
        bs_run_frac[f], err, IF(ok == 0, ">=", "<"), b, worst,
        ERR(sum, K(1.0)), IF(ok == 0, "FAIL", "OK"));
      CU_ASSERT(ok == 1);
    }

    Y(free)(run);
  }
#else
  printf("B-SPLINE RUN\n------------\nskipped: %s window\n",
    Y(get_window_name)());
#endif
}

/* B-spline phi_hut for n = 512, N = 256, m = 11, k = 0, 8, .. 128, at
 * 45 digits. Regenerate with tests/windowref. */
static const R bspline_phi_hut_ref_512_256_11[] =
{
  K(1.95312500000000000000000000000000000000000000e-3),
  K(1.93594358630009186311578288650563870965255356e-3),
  K(1.88528482941776629943034204572880157565517174e-3),
  K(1.80373089487902124776434144254828551822251216e-3),
  K(1.69534629235609927237473028524735291254398690e-3),
  K(1.56534783216229833627136086623454819426842273e-3),
  K(1.41969556338760304110678980475386601477386333e-3),
  K(1.26464859743901316525070422565733692355519045e-3),
  K(1.10633101190372100065318355064685183166551524e-3),
  K(9.50348555876092574771234895494703173290163051e-4),
  K(8.01487674996058986126093358201508137622569129e-4),
  K(6.63516121348329183487788911569770406464546355e-4),
  K(5.39091100831158274053159970819464337440066295e-4),
  K(4.29768483105561443818441647473641523891266777e-4),
  K(3.36096640830027773092193674133496053347832340e-4),
  K(2.57772004112209819983699482678435346029178899e-4),
  K(1.93830759132412245295343433343789551401816276e-4),
};

/* The evaluation carries the rounding of KPI as well as its own. Over this
 * range |t| <= pi/4, where the exponent 2m amplifies a perturbation of t by at
 * most 2m |t d(log sinc)/dt| = 3.4, so that costs under two ulp. */
void X(check_bspline_phi_hut_reference)(void)
{
  const INT m = 11, n = 512;
  const R b = tight_bound();
  unsigned int j;

  printf("B-SPLINE PHI HUT REFERENCE\n--------------------------\n");

  for (j = 0; j < SIZE(bspline_phi_hut_ref_512_256_11); j++)
  {
    const INT k = (INT)(8 * j);
    const R ref = bspline_phi_hut_ref_512_256_11[j];
    const R v = Y(bspline_phi_hut)((R)m, (R)n, (R)k);
    const R err = ERR(v, ref);
    const int ok = IF(FINITE(v) && err < b, 1, 0);

    printf("phi_hut[k=" __D__ "] = " __FE__ " err_rel = " __FE__ " %-2s " __FE__
      " -> %-4s\n", k, v, err, IF(ok == 0, ">=", "<"), b,
      IF(ok == 0, "FAIL", "OK"));
    CU_ASSERT(ok == 1);
  }
}

/* B-spline phi over one full run for n = 512, m = 11, with the node 5/16 of a
 * cell past a grid point, exactly. Regenerate with tests/windowref. */
static const R bspline_phi_ref_512_11[] =
{
  K(0.0),
  K(1.46255726667649444039432977010989116915111891e-26),
  K(2.26206706815254378402002903440249673559203781e-18),
  K(3.96423328726049195101177642945694358070951166e-14),
  K(2.95893189828169424184368140050253227154580188e-11),
  K(4.03959793477833029716510367329659454224574999e-9),
  K(1.76294252635112805078724397469070352661777913e-7),
  K(3.22236783135750553361558570780220201720795257e-6),
  K(2.86646368560905105806986148831801990663641313e-5),
  K(1.35726504120086485471928294660028230483461825e-4),
  K(3.61154774091844053591416348689386924276002766e-4),
  K(5.56905475628424193528322711193526483252857425e-4),
  K(5.04089129801325900969047549242893230987443202e-4),
  K(2.67057365274301469578621297112843567486302795e-4),
  K(8.12375513875158214806960760367346961404177694e-5),
  K(1.36514280899762779180172930711925560347649505e-5),
  K(1.18675675043218049293885956310162680560032029e-6),
  K(4.78975993762742430677902934623596881193538917e-8),
  K(7.45921482721972249638727884062332599728723649e-10),
  K(3.16656253752052469510544008518504721286688634e-12),
  K(1.69058667978236500053514027018215199202892026e-15),
  K(1.15470668751465213119007510996947578407386178e-20),
  K(9.42402832342794418047808376805973030499025120e-34),
  K(0.0),
};

#define BSPLINE_PHI_PEAK_512_11 K(5.56905475628424193528322711193526483252857425e-4)

/* Scaled by the peak of the run, which is what the psi table owes the
 * transform: relative error is not the yardstick in a tail that is 1e-30 of
 * the peak. */
void X(check_bspline_phi_reference)(void)
{
  const INT m = 11, n = 512;
  const R t = (R)(2 * m) + K(0.3125);
  const R peak = BSPLINE_PHI_PEAK_512_11;
  const R b = tight_bound();
  R inv[2 * 11];
  R run[2 * 11 + 2];
  INT j;
  unsigned int l;

  printf("B-SPLINE PHI REFERENCE\n----------------------\n");

  for (j = 0; j < 2 * m; j++)
    inv[j] = K(1.0) / (R)(j + 1);

  Y(bspline_phi_run)(run, inv, m, t, K(1.0) / (R)n);

  for (l = 0; l < SIZE(bspline_phi_ref_512_11); l++)
  {
    const R ref = bspline_phi_ref_512_11[l];
    const R err = ABS(run[l] - ref) / peak;
    const int ok = IF(FINITE(run[l]) && err < b, 1, 0);

    printf("phi[l=%u] = " __FE__ " err_peak = " __FE__ " %-2s " __FE__
      " -> %-4s\n", l, run[l], err, IF(ok == 0, ">=", "<"), b,
      IF(ok == 0, "FAIL", "OK"));
    CU_ASSERT(ok == 1);
  }
}

/* Sinc-power phi over one full run for n = 512, m = 8, sigma = 2,
 * with the node 5/16 of a cell past a grid point, at 45 digits.
 * Regenerate with tests/windowref. */
static const R sincpow_phi_ref_512_8[] =
{
  K(4.36207586065430746387232746447349724977706968e-11),
  K(2.43680743518711969357518281964372277366506755e-8),
  K(2.34293077514470900494033946774452399401838925e-6),
  K(7.26304931890152642363910375424280653843661635e-5),
  K(9.80501338998901095072353779317351126248453582e-4),
  K(6.79832785934676595981396341059056602713432049e-3),
  K(2.66797992154586308167968186733703615422369400e-2),
  K(6.28119218563758742590312193724812453766520929e-2),
  K(9.16553687098782541881852710555349615598479185e-2),
  K(8.40277875338755162795744869185860849194731528e-2),
  K(4.82504217468975899417126337655798080331273772e-2),
  K(1.70056629677184414621021078479783954399443935e-2),
  K(3.53036381578621282853178780309425694356361411e-3),
  K(4.02254150508882857825322138174645139839879033e-4),
  K(2.23572346927523570911968586768458633293762271e-5),
  K(4.94132555907791169878213439563993391417467121e-7),
  K(2.94595404023347758253385857899695838024837600e-9),
  K(1.97379694955786135670254676056348595056548479e-12),
};

#define SINCPOW_PHI_PEAK_512_8 K(9.16553687098782541881852710555349615598479185e-2)

/* m is 8 rather than the default cutoff so that w = (2 sigma - 1)/(2 m sigma)
 * is 3/32 and exactly representable; at 11 it is 3/44, whose rounding would be
 * charged to the window. The scalar evaluation is checked against the run, so
 * the batching cannot drift from what PHI gives. */
void X(check_sincpow_phi_reference)(void)
{
  const INT m = 8;
  const R sigma = K(2.0);
  const R w = (K(2.0) * sigma - K(1.0)) / (K(2.0) * (R)m * sigma);
  const R nx0 = (R)m + K(0.3125);
  const R peak = SINCPOW_PHI_PEAK_512_8;
  const R b = tight_bound();
  R run[2 * 8 + 2];
  unsigned int l;

  printf("SINC-POWER PHI REFERENCE\n------------------------\n");

  Y(sincpow_phi_run)(run, w, (R)m, m, KPI * w, nx0);

  for (l = 0; l < SIZE(sincpow_phi_ref_512_8); l++)
  {
    const R ref = sincpow_phi_ref_512_8[l];
    const R s = Y(sincpow_phi)(w, (R)m, KPI * w * (nx0 - (R)l));
    const R err = ABS(run[l] - ref) / peak;
    const int ok = IF(FINITE(run[l]) && err < b
        && ABS(run[l] - s) / peak < b, 1, 0);

    printf("phi[l=%u] = " __FE__ " err_peak = " __FE__ " %-2s " __FE__
      " -> %-4s\n", l, run[l], err, IF(ok == 0, ">=", "<"), b,
      IF(ok == 0, "FAIL", "OK"));
    CU_ASSERT(ok == 1);
  }
}

/* Sinc-power phi_hut for n = 512, N = 256, m = 8, sigma = 2,
 * k = 0, 8, .. 128, exactly. Regenerate with tests/windowref. */
static const R sincpow_phi_hut_ref_512_256_8[] =
{
  K(3.42240261355340720420085499450578815658180738e-1),
  K(3.38824153447004483507437458309156379687532505e-1),
  K(3.28774151231265859059605533349723993187870495e-1),
  K(3.12666606251760809745782502788951201649614348e-1),
  K(2.91402434821646596718351229654690575010814302e-1),
  K(2.66125543472254185414009769273148172864898261e-1),
  K(2.38123194910707311500962294613088263881914676e-1),
  K(2.08720123267057691655396762694173722991948452e-1),
  K(1.79178216565491298033718038325929439740481808e-1),
  K(1.50611949803996986842046998813698549148284598e-1),
  K(1.23926850826279051318194677196824049501956566e-1),
  K(9.97846726607076617126524811846554674654411158e-2),
  K(7.85952538667485757432847379937327027274117221e-2),
  K(6.05318412360053820176892995127200283934300358e-2),
  K(4.55643445808147538625550334425132732457154837e-2),
  K(3.35038025716498355259280229311644655559999475e-2),
  K(2.40512578451359977196034432985131716577929348e-2),
};

/* The argument carries the rounding of 1/w, which the plan stores rather than
 * dividing per evaluation. */
void X(check_sincpow_phi_hut_reference)(void)
{
  const INT m = 8, n = 512;
  const R sigma = K(2.0);
  const R w = (K(2.0) * sigma - K(1.0)) / (K(2.0) * (R)m * sigma);
  const R winv = K(1.0) / w;
  const R b = tight_bound();
  unsigned int j;

  printf("SINC-POWER PHI HUT REFERENCE\n----------------------------\n");

  for (j = 0; j < SIZE(sincpow_phi_hut_ref_512_256_8); j++)
  {
    const INT k = (INT)(8 * j);
    const R ref = sincpow_phi_hut_ref_512_256_8[j];
    const R v = Y(bsplines)(2 * m, (R)k * winv / (R)n + (R)m);
    const R err = ERR(v, ref);
    const int ok = IF(FINITE(v) && err < b, 1, 0);

    printf("phi_hut[k=" __D__ "] = " __FE__ " err_rel = " __FE__ " %-2s " __FE__
      " -> %-4s\n", k, v, err, IF(ok == 0, ">=", "<"), b,
      IF(ok == 0, "FAIL", "OK"));
    CU_ASSERT(ok == 1);
  }
}

/* The Chebyshev table must reproduce what de Boor gives, at the accuracy each
 * of them owes its caller: against the peak everywhere, because the table is
 * built for phi and psi weights, and relatively over the pieces the guard
 * admits, because sinc-power phi_hut is divided by. The tail pieces the guard
 * rejects are where B_k runs orders below its peak; there the fallback returns
 * de Boor itself, so agreement is exact. */
void X(check_bspline_cheb)(void)
{
  static const int cheb_m[] = {1, 2, 4, 8, 11, 16};
  /* against the peak this is an ordinary evaluation; relatively it is only
   * promised what the guard admits */
  const R b = K(64.0) * Y(float_property)(NFFT_EPSILON);
  const R br = NFFT_CHEB_GUARD * Y(float_property)(NFFT_EPSILON);
  unsigned int i;
  INT t;

  printf("B-SPLINE CHEBYSHEV\n------------------\n");

  for (i = 0; i < SIZE(cheb_m); i++)
  {
    const INT m = (INT)cheb_m[i], k = 2 * m;
    R *tab = (R*)Y(malloc)((size_t)(k * k) * sizeof(R));
    R peak = K(0.0), wpeak = K(0.0), wrel = K(0.0);
    INT i0, guarded = 0;
    int ok;

    Y(bspline_cheb_init)(tab, k);
    i0 = Y(bspline_cheb_guard)(tab, k, NFFT_CHEB_GUARD);

    for (t = 0; t <= 64 * k; t++)
      peak = FMAX(peak, Y(bsplines)(k, (R)t / K(64.0)));

    for (t = 1; t < 64 * k; t++)
    {
      const R x = (R)t / K(64.0);
      const R ref = Y(bsplines)(k, x);
      const R v = Y(bspline_cheb)(tab, k, x);
      const R g = Y(bspline_cheb_rel)(tab, k, i0, x);

      if (!(FINITE(v)) || !(FINITE(g)))
      {
        wpeak = K(1.0);
        break;
      }
      wpeak = FMAX(wpeak, ABS(v - ref) / peak);

      if ((INT)FLOOR(x) >= i0 && (INT)FLOOR(x) < k - i0)
        wrel = FMAX(wrel, ERR(g, ref));
      else
      {
        guarded++;
        /* outside the admitted pieces the fallback is de Boor itself */
        if (g != ref)
          wrel = K(1.0);
      }
    }

    /* the ends are always rejected: B_k reaches zero there. Everything the
     * guard does admit is contiguous around the middle, so i0 bounds it. */
    ok = IF(wpeak < b && wrel < br && i0 > 0 && i0 <= k / 2 && guarded > 0,
        1, 0);
    printf("cheb[m=" __D__ "] i0 = " __D__ " err_peak = " __FE__ " err_rel = "
      __FE__ " %-2s " __FE__ " -> %-4s\n", m, i0, wpeak, wrel,
      IF(ok == 0, ">=", "<"), br, IF(ok == 0, "FAIL", "OK"));
    CU_ASSERT(ok == 1);

    Y(free)(tab);
  }
}
