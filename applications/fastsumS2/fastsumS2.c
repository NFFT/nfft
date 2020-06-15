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

/**
 * \defgroup applications_fastsumS2_test fastsumS2_matlab
 * \ingroup applications_fastsumS2
 * \{
 */
#include "config.h"

/* standard headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

/* NFFT3 header */
#include "nfft3.h"

/* NFFT3 utilities */
#include "infft.h"

/* Fourier-Legendre coefficients for Abel-Poisson kernel */
#define SYMBOL_ABEL_POISSON(k,h) (pow(h,k))

/* Fourier-Legendre coefficients for singularity kernel */
#define SYMBOL_SINGULARITY(k,h) ((2.0/(2*k+1))*pow(h,k))

/* Flags for the different kernel functions */

/** Abel-Poisson kernel */
#define KT_ABEL_POISSON (0)
/** Singularity kernel */
#define KT_SINGULARITY  (1)
/** Locally supported kernel */
#define KT_LOC_SUPP     (2)
/** Gaussian kernel */
#define KT_GAUSSIAN     (3)

/** Enumeration type for yes/no/both-type parameters */
enum pvalue {NO = 0, YES = 1, BOTH = 2};

static inline int scaled_modified_bessel_i_series(const R x, const R alpha,
  const int nb, const int ize, R *b)
{
  const R enmten = K(4.0)*nfft_float_property(NFFT_R_MIN);
  R tempa = K(1.0), empal = K(1.0) + alpha, halfx = K(0.0), tempb = K(0.0);
  int n, ncalc = nb;

  if (enmten < x)
    halfx = x/K(2.0);

  if (alpha != K(0.0))
    tempa = POW(halfx, alpha)/TGAMMA(empal);

  if (ize == 2)
    tempa *= EXP(-x);

  if (K(1.0) < x + K(1.0))
    tempb = halfx*halfx;

  b[0] = tempa + tempa*tempb/empal;

  if (x != K(0.0) && b[0] == K(0.0))
    ncalc = 0;

  if (nb == 1)
    return ncalc;

  if (K(0.0) < x)
  {
    R tempc = halfx, tover = (enmten + enmten)/x;

    if (tempb != K(0.0))
      tover = enmten/tempb;

    for (n = 1; n < nb; n++)
    {
      tempa /= empal;
      empal += K(1.0);
      tempa *= tempc;

      if (tempa <= tover*empal)
        tempa = K(0.0);

      b[n] = tempa + tempa*tempb/empal;

      if (b[n] == K(0.0) && n < ncalc)
        ncalc = n;
    }
  }
  else
    for (n = 1; n < nb; n++)
      b[n] = K(0.0);

  return ncalc;
}

static inline void scaled_modified_bessel_i_normalize(const R x,
  const R alpha, const int nb, const int ize, R *b, const R sum_)
{
  const R enmten = K(4.0)*nfft_float_property(NFFT_R_MIN);
  R sum = sum_, tempa;
  int n;

  /* Normalize, i.e., divide all b[n] by sum */
  if (alpha != K(0.0))
    sum = sum * TGAMMA(K(1.0) + alpha) * POW(x/K(2.0), -alpha);

  if (ize == 1)
    sum *= EXP(-x);

  tempa = enmten;

  if (K(1.0) < sum)
    tempa *= sum;

  for (n = 1; n <= nb; n++)
  {
    if (b[n-1] < tempa)
      b[n-1] = K(0.0);

    b[n-1] /= sum;
  }
}

/**
 * Calculates the modified bessel function \f$I_{n+\alpha}(x)\f$, possibly
 * scaled by \f$\mathrm{e}^{-x}\f$, for real non-negative \f$x,alpha\f$ with
 * \f$0 \le \alpha < 1\f$, and \f$n=0,1,\ldots,nb-1\f$.
 *
 * \arg[in] \c x non-negative real number in \f$I_{n+\alpha}(x)\f$
 * \arg[in] \c alpha non-negative real number with \f$0 \le \alpha < 1\f$ in
 *   \f$I_{n+\alpha}(x)\f$
 * \arg[in] \c nb number of functions to be calculated
 * \arg[in] \c ize switch between no scaling (\c ize = 1) and exponential
 *   scaling (\c ize = 2)
 * \arg[out] \c b real output vector to contain \f$I_{n+\alpha}(x)\f$,
 *   \f$n=0,1,\ldots,nb-1\f$
 * \return error indicator. Only if this value is identical to \c nb, then all
 *   values in \c b have been calculated to full accuracy. If not, errors are
 *   indicated using the following scheme:
 *   - ncalc < 0: At least one of the arguments was out of range (e.g. nb <= 0,
 *     ize neither equals 1 nor 2, \f$|x| \ge exparg\f$). In this case, the
 *     output vector b is not calculated and \c ncalc is set to
 *     \f$\min(nb,0)-1\f$.
 *   - 0 < ncalc < nb: Not all requested functions could be calculated to full
 *     accuracy. This can occur when nb is much larger than |x|. in this case,
 *     the values \f$I_{n+\alpha}(x)\f$ are calculated to full accuracy for
 *     \f$n=0,1,\ldots,ncalc\f$. The rest of the values up to
 *     \f$n=0,1,\ldots,nb-1\f$ is calculated to a lower accuracy.
 *
 * \acknowledgement
 *
 * This program is based on a program written by David J. Sookne [2] that
 * computes values of the Bessel functions \f$J_{\nu}(x)\f$ or \f$I_{\nu}(x)\f$
 * for real argument \f$x\f$ and integer order \f$\nu\f$. modifications include
 * the restriction of the computation to the Bessel function \f$I_{\nu}(x)\f$
 * for non-negative real argument, the extension of the computation to arbitrary
 * non-negative orders \f$\nu\f$, and the elimination of most underflow.
 *
 * References:
 * [1] F. W. J. Olver and D. J. Sookne, A note on backward recurrence
 *   algorithms", Math. Comput. (26), 1972, pp 125 -- 132.
 * [2] D. J. Sookne, "Bessel functions of real argument and int order", NBS
 *   Jour. of Res. B. (77B), 1973, pp. 125 -- 132.
 *
 * Modified by W. J. Cody, Applied Mathematics Division, Argonne National
 *   Laboratory, Argonne, IL, 60439, USA
 *
 * Modified by Jens Keiner, Institute of Mathematics, University of Lübeck,
 *   23560 Lübeck, Germany
 */
static int smbi(const R x, const R alpha, const int nb, const int ize, R *b)
{
  /* machine dependent parameters */
  /* NSIG   - DECIMAL SIGNIFICANCE DESIRED.  SHOULD BE SET TO */
  /*          IFIX(ALOG10(2)*NBIT+1), WHERE NBIT IS THE NUMBER OF */
  /*          BITS IN THE MANTISSA OF A WORKING PRECISION VARIABLE. */
  /*          SETTING NSIG LOWER WILL RESULT IN DECREASED ACCURACY */
  /*          WHILE SETTING NSIG HIGHER WILL INCREASE CPU TIME */
  /*          WITHOUT INCREASING ACCURACY.  THE TRUNCATION ERROR */
  /*          IS LIMITED TO A RELATIVE ERROR OF T=.5*10**(-NSIG). */
  /* ENTEN  - 10.0 ** K, WHERE K IS THE LARGEST int SUCH THAT */
  /*          ENTEN IS MACHINE-REPRESENTABLE IN WORKING PRECISION. */
  /* ENSIG  - 10.0 ** NSIG. */
  /* RTNSIG - 10.0 ** (-K) FOR THE SMALLEST int K SUCH THAT */
  /*          K .GE. NSIG/4. */
  /* ENMTEN - THE SMALLEST ABS(X) SUCH THAT X/4 DOES NOT UNDERFLOW. */
  /* XLARGE - UPPER LIMIT ON THE MAGNITUDE OF X WHEN IZE=2.  BEAR */
  /*          IN MIND THAT IF ABS(X)=N, THEN AT LEAST N ITERATIONS */
  /*          OF THE BACKWARD RECURSION WILL BE EXECUTED. */
  /* EXPARG - LARGEST WORKING PRECISION ARGUMENT THAT THE LIBRARY */
  /*          EXP ROUTINE CAN HANDLE AND UPPER LIMIT ON THE */
  /*          MAGNITUDE OF X WHEN IZE=1. */
  const int nsig = MANT_DIG + 2;
  const R enten = nfft_float_property(NFFT_R_MAX);
  const R ensig = POW(K(10.0),(R)nsig);
  const R rtnsig = POW(K(10.0),-CEIL((R)nsig/K(4.0)));
  const R xlarge = K(1E4);
  const R exparg = FLOOR(LOG(POW(K(R_RADIX),K(DBL_MAX_EXP-1))));

  /* System generated locals */
  int l, n, nend, magx, nbmx, ncalc, nstart;
  R p, em, en, sum, pold, test, empal, tempa, tempb, tempc, psave, plast, tover,
    emp2al, psavel;

  magx = LRINT(FLOOR(x));

  /* return if x, nb, or ize out of range */
  if (   nb <= 0 || x < K(0.0) || alpha < K(0.0) || K(1.0) <= alpha
      || ((ize != 1 || exparg < x) && (ize != 2 || xlarge < x)))
    return (MIN(nb,0) - 1);

  /* 2-term ascending series for small x */
  if (x < rtnsig)
    return scaled_modified_bessel_i_series(x,alpha,nb,ize,b);

  ncalc = nb;
  /* forward sweep, Olver's p-sequence */

  nbmx = nb - magx;
  n = magx + 1;

  en = (R) (n+n) + (alpha+alpha);
  plast = K(1.0);
  p = en/x;

  /* significance test */
  test = ensig + ensig;

  if ((5*nsig) < (magx << 1))
    test = SQRT(test*p);
  else
    test /= POW(K(1.585),(R)magx);

  if (3 <= nbmx)
  {
    /* calculate p-sequence until n = nb-1 */
    tover = enten/ensig;
    nstart = magx+2;
    nend = nb - 1;

    for (n = nstart; n <= nend; n++)
    {
      en += K(2.0);
      pold = plast;
      plast = p;
      p = en*plast/x + pold;
      if (p > tover)
      {
        /* divide p-sequence by tover to avoid overflow. Calculate p-sequence
         * until 1 <= |p| */
        tover = enten;
        p /= tover;
        plast /= tover;
        psave = p;
        psavel = plast;
        nstart = n + 1;

        do
        {
          n++;
          en += K(2.0);
          pold = plast;
          plast = p;
          p = en*plast/x + pold;
        } while (p <= K(1.0));

        tempb = en/x;

        /* Backward test. Find ncalc as the largest n such that test is passed. */
        test = pold*plast*(K(0.5) - K(0.5)/(tempb * tempb))/ensig;
        p = plast*tover;
        n--;
        en -= K(2.0);
        nend = MIN(nb,n);

        for (ncalc = nstart; ncalc <= nend; ncalc++)
        {
          pold = psavel;
          psavel = psave;
          psave = en*psavel/x + pold;
          if (test < psave * psavel)
            break;
        }

        ncalc--;
        goto L80;
      }
    }

    n = nend;
    en = (R) (n+n) + (alpha+alpha);

    /* special significance test for 2 <= nbmx */
    test = FMAX(test,SQRT(plast*ensig)*SQRT(p+p));
  }

  /* calculate p-sequence until significance test is passed */
  do
  {
    n++;
    en += K(2.0);
    pold = plast;
    plast = p;
    p = en*plast/x + pold;
  } while (p < test);

  /* Initialize backward recursion and normalization sum. */
L80:
  n++;
  en += K(2.0);
  tempb = K(0.0);
  tempa = K(1.0)/p;
  em = (R)(n-1);
  empal = em + alpha;
  emp2al = em - K(1.0) + (alpha+alpha);
  sum = tempa*empal*emp2al/em;
  nend = n-nb;

  if (nend < 0)
  {
    /* We have n <= nb. So store b[n] and set higher orders to zero */
    b[n-1] = tempa;
    nend = -nend;
    for (l = 1; l <= nend; ++l)
      b[n-1 + l] = K(0.0);
  }
  else
  {
    if (nend != 0)
    {
      /* recur backward via difference equation, calculating b[n] until n = nb */
      for (l = 1; l <= nend; ++l)
      {
        n--;
        en -= K(2.0);
        tempc = tempb;
        tempb = tempa;
        tempa = en*tempb/x + tempc;
        em -= K(1.0);
        emp2al -= K(1.0);

        if (n == 1)
          break;

        if (n == 2)
          emp2al = K(1.0);

        empal -= K(1.0);
        sum = (sum + tempa*empal)*emp2al/em;
      }
    }

    /* store b[nb] */
    b[n-1] = tempa;

    if (nb <= 1)
    {
      sum = sum + sum + tempa;
      scaled_modified_bessel_i_normalize(x,alpha,nb,ize,b,sum);
      return ncalc;
    }

    /* calculate and store b[nb-1] */
    n--;
    en -= 2.0;
    b[n-1] = en*tempa/x + tempb;

    if (n == 1)
    {
      sum = sum + sum + b[0];
      scaled_modified_bessel_i_normalize(x,alpha,nb,ize,b,sum);
      return ncalc;
    }

    em -= K(1.0);
    emp2al -= K(1.0);

    if (n == 2)
      emp2al = K(1.0);

    empal -= K(1.0);
    sum = (sum + b[n-1]*empal)*emp2al/em;
  }

  nend = n - 2;

  if (nend != 0)
  {
    /* Calculate and store b[n] until n = 2. */
    for (l = 1; l <= nend; ++l)
    {
      n--;
      en -= K(2.0);
      b[n-1] = en*b[n]/x + b[n+1];
      em -= K(1.0);
      emp2al -= K(1.0);

      if (n == 2)
        emp2al = K(1.0);

      empal -= K(1.0);
      sum = (sum + b[n-1]*empal)*emp2al/em;
    }
  }

  /* calculate b[1] */
  b[0] = K(2.0)*empal*b[1]/x + b[2];
  sum = sum + sum + b[0];

  scaled_modified_bessel_i_normalize(x,alpha,nb,ize,b,sum);
  return ncalc;
}

/**
 * Computes the \f$\mathbb{R}^3\f$ standard inner product between two vectors
 * on the unit sphere \f$\mathbb{S}^2\f$ given in spherical coordinates.
 *
 * \arg phi1   The angle \f$\varphi_1 \in [-\pi,\pi)\f$ of the first vector
 * \arg theta1 The angle \f$\vartheta_1 \in [0,\pi]\f$ of the first vector
 * \arg phi2   The angle \f$\varphi_2 \in [-\pi,\pi)\f$ of the second vector
 * \arg theta2 The angle \f$\vartheta_2 \in [0,\pi]\f$ of the second vector
 *
 * \return The inner product \f$\cos \vartheta_1 \cos \vartheta_2 +
 *         \sin \vartheta_1 \sin(\vartheta_2 \cos(\varphi_1 - \varphi_2)\f$
 *
 * \author Jens Keiner
 */
static inline double innerProduct(const double phi1, const double theta1,
  const double phi2, const double theta2)
{
  double pi2theta1 = K2PI*theta1, pi2theta2 = K2PI*theta2;
  return (cos(pi2theta1)*cos(pi2theta2)
    + sin(pi2theta1)*sin(pi2theta2)*cos(K2PI*(phi1-phi2)));
}

/**
 * Evaluates the Poisson kernel \f$Q_h: [-1,1] \rightarrow \mathbb{R}\f$ at a
 * node \f$x \in [-1,1]\f$.
 *
 * \arg x The node \f$x \in [-1,1]\f$
 * \arg h The parameter \f$h \in (0,1)\f$
 *
 * \return The value of the Poisson kernel \f$Q_h(x)\f$ at the node \f$x\f$
 *
 * \author Jens Keiner
 */
static inline double poissonKernel(const double x, const double h)
{
  return (1.0/(K4PI))*((1.0-h)*(1.0+h))/pow(sqrt(1.0-2.0*h*x+h*h),3.0);
}

/**
 * Evaluates the singularity kernel \f$S_h: [-1,1] \rightarrow \mathbb{R}\f$ at
 * a node \f$x \in [-1,1]\f$.
 *
 * \arg x The node \f$x \in [-1,1]\f$
 * \arg h The parameter \f$h \in (0,1)\f$
 *
 * \return The value of the Poisson kernel \f$S_h(x)\f$ at the node \f$x\f$
 *
 * \author Jens Keiner
 */
static inline double singularityKernel(const double x, const double h)
{
  return (1.0/(K2PI))/sqrt(1.0-2.0*h*x+h*h);
}

/**
 * Evaluates the locally supported kernel \f$L_{h,\lambda}: [-1,1] \rightarrow
 * \mathbb{R}\f$ at a node \f$x \in [-1,1]\f$.
 *
 * \arg x The node \f$x \in [-1,1]\f$
 * \arg h The parameter \f$h \in (0,1)\f$
 * \arg lambda The parameter \f$\lambda \in \mathbb{N}_0\f$
 *
 * \return The value of the locally supported kernel \f$L_{h,\lambda}(x)\f$ at
 *   the node \f$x\f$
 *
 * \author Jens Keiner
 */
static inline double locallySupportedKernel(const double x, const double h,
  const double lambda)
{
  return (x<=h)?(0.0):(pow((x-h),lambda));
}

/**
 * Evaluates the spherical Gaussian kernel \f$G_\sigma: [-1,1] \rightarrow
 * \mathbb{R}\f$ at a node \f$x \in [-1,1]\f$.
 *
 * \arg x The node \f$x \in [-1,1]\f$
 * \arg sigma The parameter \f$\sigma \in \mathbb{R}_+\f$
 *
 * \return The value of the pherical Gaussian kernel \f$G_\sigma(x)\f$ at the
 *   node \f$x\f$
 *
 * \author Jens Keiner
 */
static inline double gaussianKernel(const double x, const double sigma)
{
   return exp(2.0*sigma*(x-1.0));
}

/**
 * The main program.
 *
 * \param argc The number of arguments
 * \param argv An array containing the arguments as C-strings
 *
 * \return Exit code
 *
 * \author Jens Keiner
 */
int main (int argc, char **argv)
{
  double **p;                  /* The array containing the parameter sets     *
                                * for the kernel functions                    */
  int *m;                      /* The array containing the cut-off degrees M  */
  int **ld;                    /* The array containing the numbers of source  *
                                * and target nodes, L and D                   */
  int ip;                      /* Index variable for p                        */
  int im;                      /* Index variable for m                        */
  int ild;                     /* Index variable for l                        */
  int ipp;                     /* Index for kernel parameters                 */
  int ip_max;                  /* The maximum index for p                     */
  int im_max;                  /* The maximum index for m                     */
  int ild_max;                 /* The maximum index for l                     */
  int ipp_max;                 /* The maximum index for ip                    */
  int tc_max;                  /* The number of testcases                     */
  int m_max;                   /* The maximum cut-off degree M for the        *
                                * current dataset                             */
  int l_max;                   /* The maximum number of source nodes L for    *
                                * the current dataset                         */
  int d_max;                   /* The maximum number of target nodes D for    *
                                * the current dataset                         */
  long ld_max_prec;            /* The maximum number of source and target     *
                                * nodes for precomputation multiplied         */
  long l_max_prec;             /* The maximum number of source nodes for      *
                                * precomputation                              */
  int tc;                      /* Index variable for testcases                */
  int kt;                      /* The kernel function                         */
  int cutoff;                  /* The current NFFT cut-off parameter          */
  double threshold;            /* The current NFSFT threshold parameter       */
  double t_d;                  /* Time for direct algorithm in seconds        */
  double t_dp;                 /* Time for direct algorithm with              *
                                  precomputation in seconds                   */
  double t_fd;                 /* Time for fast direct algorithm in seconds   */
  double t_f;                  /* Time for fast algorithm in seconds          */
  double temp;                 /*                                             */
  double err_f;                /* Error E_infty for fast algorithm            */
  double err_fd;               /* Error E_\infty for fast direct algorithm    */
  ticks t0, t1;                /*                                             */
  int precompute = NO;         /*                                             */
  fftw_complex *ptr;         /*                                             */
  double* steed;               /*                                             */
  fftw_complex *b;           /* The weights (b_l)_{l=0}^{L-1}               */
  fftw_complex *f_hat;       /* The spherical Fourier coefficients          */
  fftw_complex *a;           /* The Fourier-Legendre coefficients           */
  double *xi;                  /* Target nodes                                */
  double *eta;                 /* Source nodes                                */
  fftw_complex *f_m;         /* Approximate function values                 */
  fftw_complex *f;           /* Exact function values                       */
  fftw_complex *prec = NULL; /*                                             */
  nfsft_plan plan;             /* NFSFT plan                                  */
  nfsft_plan plan_adjoint;     /* adjoint NFSFT plan                          */
  int i;                       /*                                             */
  int k;                       /*                                             */
  int n;                       /*                                             */
  int d;                       /*                                             */
  int l;                       /*                                             */
  int use_nfsft;               /*                                             */
  int use_nfft;                /*                                             */
  int use_fpt;                 /*                                             */
  int rinc;                    /*                                             */
  double constant;             /*                                             */

  /* Read the number of testcases. */
  fscanf(stdin,"testcases=%d\n",&tc_max);
  fprintf(stdout,"%d\n",tc_max);

  /* Process each testcase. */
  for (tc = 0; tc < tc_max; tc++)
  {
    /* Check if the fast transform shall be used. */
    fscanf(stdin,"nfsft=%d\n",&use_nfsft);
    fprintf(stdout,"%d\n",use_nfsft);
    if (use_nfsft != NO)
    {
      /* Check if the NFFT shall be used. */
      fscanf(stdin,"nfft=%d\n",&use_nfft);
      fprintf(stdout,"%d\n",use_nfft);
      if (use_nfft != NO)
      {
        /* Read the cut-off parameter. */
        fscanf(stdin,"cutoff=%d\n",&cutoff);
        fprintf(stdout,"%d\n",cutoff);
      }
      else
      {
        /* TODO remove this */
        /* Initialize unused variable with dummy value. */
        cutoff = 1;
      }
      /* Check if the fast polynomial transform shall be used. */
      fscanf(stdin,"fpt=%d\n",&use_fpt);
      fprintf(stdout,"%d\n",use_fpt);
      /* Read the NFSFT threshold parameter. */
      fscanf(stdin,"threshold=%lf\n",&threshold);
      fprintf(stdout,"%lf\n",threshold);
    }
    else
    {
      /* TODO remove this */
      /* Set dummy values. */
      cutoff = 3;
      threshold = 1000000000000.0;
    }

    /* Initialize bandwidth bound. */
    m_max = 0;
    /* Initialize source nodes bound. */
    l_max = 0;
    /* Initialize target nodes bound. */
    d_max = 0;
    /* Initialize source nodes bound for precomputation. */
    l_max_prec = 0;
    /* Initialize source and target nodes bound for precomputation. */
    ld_max_prec = 0;

    /* Read the kernel type. This is one of KT_ABEL_POISSON, KT_SINGULARITY,
     * KT_LOC_SUPP and KT_GAUSSIAN. */
    fscanf(stdin,"kernel=%d\n",&kt);
    fprintf(stdout,"%d\n",kt);

    /* Read the number of parameter sets. */
    fscanf(stdin,"parameter_sets=%d\n",&ip_max);
    fprintf(stdout,"%d\n",ip_max);

    /* Allocate memory for pointers to parameter sets. */
    p = (double**) nfft_malloc(ip_max*sizeof(double*));

    /* We now read in the parameter sets. */

    /* Read number of parameters. */
    fscanf(stdin,"parameters=%d\n",&ipp_max);
    fprintf(stdout,"%d\n",ipp_max);

    for (ip = 0; ip < ip_max; ip++)
    {
      /* Allocate memory for the parameters. */
      p[ip] = (double*) nfft_malloc(ipp_max*sizeof(double));

      /* Read the parameters. */
      for (ipp = 0; ipp < ipp_max; ipp++)
      {
        /* Read the next parameter. */
        fscanf(stdin,"%lf\n",&p[ip][ipp]);
        fprintf(stdout,"%lf\n",p[ip][ipp]);
      }
    }

    /* Read the number of cut-off degrees. */
    fscanf(stdin,"bandwidths=%d\n",&im_max);
    fprintf(stdout,"%d\n",im_max);
    m = (int*) nfft_malloc(im_max*sizeof(int));

    /* Read the cut-off degrees. */
    for (im = 0; im < im_max; im++)
    {
      /* Read cut-off degree. */
      fscanf(stdin,"%d\n",&m[im]);
      fprintf(stdout,"%d\n",m[im]);
      m_max = MAX(m_max,m[im]);
    }

    /* Read number of node specifications. */
    fscanf(stdin,"node_sets=%d\n",&ild_max);
    fprintf(stdout,"%d\n",ild_max);
    ld = (int**) nfft_malloc(ild_max*sizeof(int*));

    /* Read the run specification. */
    for (ild = 0; ild < ild_max; ild++)
    {
      /* Allocate memory for the run parameters. */
      ld[ild] = (int*) nfft_malloc(5*sizeof(int));

      /* Read number of source nodes. */
      fscanf(stdin,"L=%d ",&ld[ild][0]);
      fprintf(stdout,"%d\n",ld[ild][0]);
      l_max = MAX(l_max,ld[ild][0]);

      /* Read number of target nodes. */
      fscanf(stdin,"D=%d ",&ld[ild][1]);
      fprintf(stdout,"%d\n",ld[ild][1]);
      d_max = MAX(d_max,ld[ild][1]);

      /* Determine whether direct and fast algorithm shall be compared. */
      fscanf(stdin,"compare=%d ",&ld[ild][2]);
      fprintf(stdout,"%d\n",ld[ild][2]);

      /* Check if precomputation for the direct algorithm is used. */
      if (ld[ild][2] == YES)
      {
        /* Read whether the precomputed version shall also be used. */
        fscanf(stdin,"precomputed=%d\n",&ld[ild][3]);
        fprintf(stdout,"%d\n",ld[ild][3]);

        /* Read the number of repetitions over which measurements are
         * averaged. */
        fscanf(stdin,"repetitions=%d\n",&ld[ild][4]);
        fprintf(stdout,"%d\n",ld[ild][4]);

        /* Update ld_max_prec and l_max_prec. */
        if (ld[ild][3] == YES)
        {
          /* Update ld_max_prec. */
          ld_max_prec = MAX(ld_max_prec,ld[ild][0]*ld[ild][1]);
          /* Update l_max_prec. */
          l_max_prec = MAX(l_max_prec,ld[ild][0]);
          /* Turn on the precomputation for the direct algorithm. */
          precompute = YES;
        }
      }
      else
      {
        /* Set default value for the number of repetitions. */
        ld[ild][4] = 1;
      }
    }

    /* Allocate memory for data structures. */
    b = (fftw_complex*) nfft_malloc(l_max*sizeof(fftw_complex));
    eta = (double*) nfft_malloc(2*l_max*sizeof(double));
    f_hat = (fftw_complex*) nfft_malloc(NFSFT_F_HAT_SIZE(m_max)*sizeof(fftw_complex));
    a = (fftw_complex*) nfft_malloc((m_max+1)*sizeof(fftw_complex));
    xi = (double*) nfft_malloc(2*d_max*sizeof(double));
    f_m = (fftw_complex*) nfft_malloc(d_max*sizeof(fftw_complex));
    f = (fftw_complex*) nfft_malloc(d_max*sizeof(fftw_complex));

    /* Allocate memory for precomputed data. */
    if (precompute == YES)
    {
      prec = (fftw_complex*) nfft_malloc(ld_max_prec*sizeof(fftw_complex));
    }

    /* Generate random source nodes and weights. */
    for (l = 0; l < l_max; l++)
    {
      b[l] = (((double)rand())/RAND_MAX) - 0.5;
      eta[2*l] = (((double)rand())/RAND_MAX) - 0.5;
      eta[2*l+1] = acos(2.0*(((double)rand())/RAND_MAX) - 1.0)/(K2PI);
    }

    /* Generate random target nodes. */
    for (d = 0; d < d_max; d++)
    {
      xi[2*d] = (((double)rand())/RAND_MAX) - 0.5;
      xi[2*d+1] = acos(2.0*(((double)rand())/RAND_MAX) - 1.0)/(K2PI);
    }

    /* Do precomputation. */
    nfsft_precompute(m_max,threshold,
      ((use_nfsft==NO)?(NFSFT_NO_FAST_ALGORITHM):(0U/*NFSFT_NO_DIRECT_ALGORITHM*/)), 0U);

    /* Process all parameter sets. */
    for (ip = 0; ip < ip_max; ip++)
    {
      /* Compute kernel coeffcients up to the maximum cut-off degree m_max. */
      switch (kt)
      {
        case KT_ABEL_POISSON:
          /* Compute Fourier-Legendre coefficients for the Poisson kernel. */
          for (k = 0; k <= m_max; k++)
            a[k] = SYMBOL_ABEL_POISSON(k,p[ip][0]);
          break;

        case KT_SINGULARITY:
          /* Compute Fourier-Legendre coefficients for the singularity
           * kernel. */
          for (k = 0; k <= m_max; k++)
            a[k] = SYMBOL_SINGULARITY(k,p[ip][0]);
          break;

        case KT_LOC_SUPP:
          /* Compute Fourier-Legendre coefficients for the locally supported
           * kernel. */
          a[0] = 1.0;
          if (1 <= m_max)
            a[1] = ((p[ip][1]+1+p[ip][0])/(p[ip][1]+2.0))*a[0];
          for (k = 2; k <= m_max; k++)
            a[k] = (1.0/(k+p[ip][1]+1))*((2*k-1)*p[ip][0]*a[k-1] -
              (k-p[ip][1]-2)*a[k-2]);
          break;

        case KT_GAUSSIAN:
          /* Fourier-Legendre coefficients */
          steed = (double*) nfft_malloc((m_max+1)*sizeof(double));
          smbi(2.0*p[ip][0],0.5,m_max+1,2,steed);
          for (k = 0; k <= m_max; k++)
            a[k] = K2PI*(sqrt(KPI/p[ip][0]))*steed[k];

          nfft_free(steed);
          break;
      }

      /* Normalize Fourier-Legendre coefficients. */
      for (k = 0; k <= m_max; k++)
        a[k] *= (2*k+1)/(K4PI);

      /* Process all node sets. */
      for (ild = 0; ild < ild_max; ild++)
      {
        /* Check if the fast algorithm shall be used. */
        if (ld[ild][2] != NO)
        {
          /* Check if the direct algorithm with precomputation should be
           * tested. */
          if (ld[ild][3] != NO)
          {
            /* Get pointer to start of data. */
            ptr = prec;
            /* Calculate increment from one row to the next. */
            rinc = l_max_prec-ld[ild][0];

            /* Process al target nodes. */
            for (d = 0; d < ld[ild][1]; d++)
            {
              /* Process all source nodes. */
              for (l = 0; l < ld[ild][0]; l++)
              {
                /* Compute inner product between current source and target
                 * node. */
                temp = innerProduct(eta[2*l],eta[2*l+1],xi[2*d],xi[2*d+1]);

                /* Switch by the kernel type. */
                switch (kt)
                {
                  case KT_ABEL_POISSON:
                    /* Evaluate the Poisson kernel for the current value. */
                    *ptr++ = poissonKernel(temp,p[ip][0]);
                   break;

                  case KT_SINGULARITY:
                    /* Evaluate the singularity kernel for the current
                     * value. */
                    *ptr++ = singularityKernel(temp,p[ip][0]);
                    break;

                  case KT_LOC_SUPP:
                     /* Evaluate the localized kernel for the current
                      * value. */
                    *ptr++ = locallySupportedKernel(temp,p[ip][0],p[ip][1]);
                    break;

                    case KT_GAUSSIAN:
                       /* Evaluate the spherical Gaussian kernel for the current
                        * value. */
                      *ptr++ = gaussianKernel(temp,p[ip][0]);
                       break;
                }
              }
              /* Increment pointer for next row. */
              ptr += rinc;
            }

            /* Initialize cumulative time variable. */
            t_dp = 0.0;

            /* Initialize time measurement. */
            t0 = getticks();

            /* Cycle through all runs. */
            for (i = 0; i < ld[ild][4]; i++)
            {

              /* Reset pointer to start of precomputed data. */
              ptr = prec;
              /* Calculate increment from one row to the next. */
              rinc = l_max_prec-ld[ild][0];

              /* Check if the localized kernel is used. */
              if (kt == KT_LOC_SUPP)
              {
                /* Perform final summation */

                /* Calculate the multiplicative constant. */
                constant = ((p[ip][1]+1)/(K2PI*pow(1-p[ip][0],p[ip][1]+1)));

                /* Process all target nodes. */
                for (d = 0; d < ld[ild][1]; d++)
                {
                  /* Initialize function value. */
                  f[d] = 0.0;

                  /* Process all source nodes. */
                  for (l = 0; l < ld[ild][0]; l++)
                    f[d] += b[l]*(*ptr++);

                  /* Multiply with the constant. */
                  f[d] *= constant;

                  /* Proceed to next row. */
                  ptr += rinc;
                }
              }
              else
              {
                /* Process all target nodes. */
                for (d = 0; d < ld[ild][1]; d++)
                {
                  /* Initialize function value. */
                  f[d] = 0.0;

                  /* Process all source nodes. */
                  for (l = 0; l < ld[ild][0]; l++)
                    f[d] += b[l]*(*ptr++);

                  /* Proceed to next row. */
                  ptr += rinc;
                }
              }
            }

            /* Calculate the time needed. */
            t1 = getticks();
            t_dp = nfft_elapsed_seconds(t1,t0);

            /* Calculate average time needed. */
            t_dp = t_dp/((double)ld[ild][4]);
          }
          else
          {
            /* Initialize cumulative time variable with dummy value. */
            t_dp = -1.0;
          }

          /* Initialize cumulative time variable. */
          t_d = 0.0;

          /* Initialize time measurement. */
          t0 = getticks();

          /* Cycle through all runs. */
          for (i = 0; i < ld[ild][4]; i++)
          {
            /* Switch by the kernel type. */
            switch (kt)
            {
              case KT_ABEL_POISSON:

                /* Process all target nodes. */
                for (d = 0; d < ld[ild][1]; d++)
                {
                  /* Initialize function value. */
                  f[d] = 0.0;

                  /* Process all source nodes. */
                  for (l = 0; l < ld[ild][0]; l++)
                  {
                    /* Compute the inner product for the current source and
                     * target nodes. */
                    temp = innerProduct(eta[2*l],eta[2*l+1],xi[2*d],xi[2*d+1]);

                    /* Evaluate the Poisson kernel for the current value and add
                     * to the result. */
                    f[d] += b[l]*poissonKernel(temp,p[ip][0]);
                  }
                }
                break;

              case KT_SINGULARITY:
                /* Process all target nodes. */
                for (d = 0; d < ld[ild][1]; d++)
                {
                  /* Initialize function value. */
                  f[d] = 0.0;

                  /* Process all source nodes. */
                  for (l = 0; l < ld[ild][0]; l++)
                  {
                    /* Compute the inner product for the current source and
                     * target nodes. */
                    temp = innerProduct(eta[2*l],eta[2*l+1],xi[2*d],xi[2*d+1]);

                    /* Evaluate the Poisson kernel for the current value and add
                     * to the result. */
                    f[d] += b[l]*singularityKernel(temp,p[ip][0]);
                  }
                }
                break;

              case KT_LOC_SUPP:
                /* Calculate the multiplicative constant. */
                constant = ((p[ip][1]+1)/(K2PI*pow(1-p[ip][0],p[ip][1]+1)));

                /* Process all target nodes. */
                for (d = 0; d < ld[ild][1]; d++)
                {
                  /* Initialize function value. */
                  f[d] = 0.0;

                  /* Process all source nodes. */
                  for (l = 0; l < ld[ild][0]; l++)
                  {
                    /* Compute the inner product for the current source and
                     * target nodes. */
                    temp = innerProduct(eta[2*l],eta[2*l+1],xi[2*d],xi[2*d+1]);

                    /* Evaluate the Poisson kernel for the current value and add
                     * to the result. */
                    f[d] += b[l]*locallySupportedKernel(temp,p[ip][0],p[ip][1]);
                  }

                  /* Multiply result with constant. */
                  f[d] *= constant;
                }
                break;

                case KT_GAUSSIAN:
                  /* Process all target nodes. */
                  for (d = 0; d < ld[ild][1]; d++)
                  {
                    /* Initialize function value. */
                    f[d] = 0.0;

                    /* Process all source nodes. */
                    for (l = 0; l < ld[ild][0]; l++)
                    {
                      /* Compute the inner product for the current source and
                       * target nodes. */
                      temp = innerProduct(eta[2*l],eta[2*l+1],xi[2*d],xi[2*d+1]);
                      /* Evaluate the Poisson kernel for the current value and add
                       * to the result. */
                      f[d] += b[l]*gaussianKernel(temp,p[ip][0]);
                    }
                  }
                  break;
            }
          }

          /* Calculate and add the time needed. */
          t1 = getticks();
          t_d = nfft_elapsed_seconds(t1,t0);
          /* Calculate average time needed. */
          t_d = t_d/((double)ld[ild][4]);
        }
        else
        {
          /* Initialize cumulative time variable with dummy value. */
          t_d = -1.0;
          t_dp = -1.0;
        }

        /* Initialize error and cumulative time variables for the fast
         * algorithm. */
        err_fd = -1.0;
        err_f = -1.0;
        t_fd = -1.0;
        t_f = -1.0;

        /* Process all cut-off bandwidths. */
        for (im = 0; im < im_max; im++)
        {
          /* Init transform plans. */
          nfsft_init_guru(&plan_adjoint, m[im],ld[ild][0],
            ((use_nfft!=0)?(0U):(NFSFT_USE_NDFT)) |
            ((use_fpt!=0)?(0U):(NFSFT_USE_DPT)),
            PRE_PHI_HUT | PRE_PSI | FFTW_INIT,
            cutoff);
          nfsft_init_guru(&plan,m[im],ld[ild][1],
            ((use_nfft!=0)?(0U):(NFSFT_USE_NDFT)) |
            ((use_fpt!=0)?(0U):(NFSFT_USE_DPT)),
            PRE_PHI_HUT | PRE_PSI | FFTW_INIT,
             cutoff);
          plan_adjoint.f_hat = f_hat;
          plan_adjoint.x = eta;
          plan_adjoint.f = b;
          plan.f_hat = f_hat;
          plan.x = xi;
          plan.f = f_m;
          nfsft_precompute_x(&plan_adjoint);
          nfsft_precompute_x(&plan);

          /* Check if direct algorithm shall also be tested. */
          if (use_nfsft == BOTH)
          {
            /* Initialize cumulative time variable. */
            t_fd = 0.0;

            /* Initialize time measurement. */
            t0 = getticks();

            /* Cycle through all runs. */
            for (i = 0; i < ld[ild][4]; i++)
            {

              /* Execute adjoint direct NDSFT transformation. */
              nfsft_adjoint_direct(&plan_adjoint);

              /* Multiplication with the Fourier-Legendre coefficients. */
              for (k = 0; k <= m[im]; k++)
                for (n = -k; n <= k; n++)
                  f_hat[NFSFT_INDEX(k,n,&plan_adjoint)] *= a[k];

              /* Execute direct NDSFT transformation. */
              nfsft_trafo_direct(&plan);

            }

            /* Calculate and add the time needed. */
            t1 = getticks();
            t_fd = nfft_elapsed_seconds(t1,t0);

            /* Calculate average time needed. */
            t_fd = t_fd/((double)ld[ild][4]);

            /* Check if error E_infty should be computed. */
            if (ld[ild][2] != NO)
            {
              /* Compute the error E_infinity. */
              err_fd = X(error_l_infty_1_complex)(f, f_m, ld[ild][1], b,
                ld[ild][0]);
            }
          }

          /* Check if the fast NFSFT algorithm shall also be tested. */
          if (use_nfsft != NO)
          {
            /* Initialize cumulative time variable for the NFSFT algorithm. */
            t_f = 0.0;
          }
          else
          {
            /* Initialize cumulative time variable for the direct NDSFT
             * algorithm. */
            t_fd = 0.0;
          }

          /* Initialize time measurement. */
          t0 = getticks();

          /* Cycle through all runs. */
          for (i = 0; i < ld[ild][4]; i++)
          {
            /* Check if the fast NFSFT algorithm shall also be tested. */
            if (use_nfsft != NO)
            {
              /* Execute the adjoint NFSFT transformation. */
              nfsft_adjoint(&plan_adjoint);
            }
            else
            {
              /* Execute the adjoint direct NDSFT transformation. */
              nfsft_adjoint_direct(&plan_adjoint);
            }

            /* Multiplication with the Fourier-Legendre coefficients. */
            for (k = 0; k <= m[im]; k++)
              for (n = -k; n <= k; n++)
                f_hat[NFSFT_INDEX(k,n,&plan_adjoint)] *= a[k];

            /* Check if the fast NFSFT algorithm shall also be tested. */
            if (use_nfsft != NO)
            {
              /* Execute the NFSFT transformation. */
              nfsft_trafo(&plan);
            }
            else
            {
              /* Execute the NDSFT transformation. */
              nfsft_trafo_direct(&plan);
            }
          }

          /* Check if the fast NFSFT algorithm has been used. */
          t1 = getticks();

          if (use_nfsft != NO)
            t_f = nfft_elapsed_seconds(t1,t0);
          else
            t_fd = nfft_elapsed_seconds(t1,t0);

          /* Check if the fast NFSFT algorithm has been used. */
          if (use_nfsft != NO)
          {
            /* Calculate average time needed. */
            t_f = t_f/((double)ld[ild][4]);
          }
          else
          {
            /* Calculate average time needed. */
            t_fd = t_fd/((double)ld[ild][4]);
          }

          /* Check if error E_infty should be computed. */
          if (ld[ild][2] != NO)
          {
            /* Check if the fast NFSFT algorithm has been used. */
            if (use_nfsft != NO)
            {
              /* Compute the error E_infinity. */
              err_f = X(error_l_infty_1_complex)(f, f_m, ld[ild][1], b,
                ld[ild][0]);
            }
            else
            {
              /* Compute the error E_infinity. */
              err_fd = X(error_l_infty_1_complex)(f, f_m, ld[ild][1], b,
                ld[ild][0]);
            }
          }

          /* Print out the error measurements. */
          fprintf(stdout,"%e\n%e\n%e\n%e\n%e\n%e\n\n",t_d,t_dp,t_fd,t_f,err_fd,
            err_f);

          /* Finalize the NFSFT plans */
          nfsft_finalize(&plan_adjoint);
          nfsft_finalize(&plan);
        } /* for (im = 0; im < im_max; im++) - Process all cut-off
           * bandwidths.*/
      } /* for (ild = 0; ild < ild_max; ild++) - Process all node sets. */
    } /* for (ip = 0; ip < ip_max; ip++) - Process all parameter sets. */

    /* Delete precomputed data. */
    nfsft_forget();

    /* Check if memory for precomputed data of the matrix K has been
     * allocated. */
    if (precompute == YES)
    {
      /* Free memory for precomputed matrix K. */
      nfft_free(prec);
    }
    /* Free data arrays. */
    nfft_free(f);
    nfft_free(f_m);
    nfft_free(xi);
    nfft_free(eta);
    nfft_free(a);
    nfft_free(f_hat);
    nfft_free(b);

    /* Free memory for node sets. */
    for (ild = 0; ild < ild_max; ild++)
      nfft_free(ld[ild]);
    nfft_free(ld);

    /* Free memory for cut-off bandwidths. */
    nfft_free(m);

    /* Free memory for parameter sets. */
    for (ip = 0; ip < ip_max; ip++)
      nfft_free(p[ip]);
    nfft_free(p);
  } /* for (tc = 0; tc < tc_max; tc++) - Process each testcase. */

  /* Return exit code for successful run. */
  return EXIT_SUCCESS;
}
/* \} */
