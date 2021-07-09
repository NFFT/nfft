/*
 * Copyright (c) 2002, 2019 Jens Keiner, Stefan Kunis, Daniel Potts
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
 * \file nfsft.c
 * \brief Implementation file for the NFSFT module
 * \author Jens Keiner
 */

#include "config.h"

/* Include standard C headers. */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/* Include NFFT3 utilities header. */

/* Include NFFT3 library header. */
#include "nfft3.h"

#include "infft.h"

/* Include private associated Legendre functions header. */
#include "legendre.h"

/* Include private API header. */
#include "api.h"

#ifdef _OPENMP
#include "../fpt/fpt.h"
#endif

/** \addtogroup nfsft
 * \{
 */

/**
 * The default NFFT cutoff parameter
 *
 * \author Jens Keiner
 */
#define NFSFT_DEFAULT_NFFT_CUTOFF 6

/**
 * The default threshold for the FPT
 *
 * \author Jens Keiner
 */
#define NFSFT_DEFAULT_THRESHOLD 1000

/**
 * The break-even bandwidth \f$N \in \mathbb{N}_0\f$
 *
 * \author Jens Keiner
 */
#define NFSFT_BREAK_EVEN 5

/**
 * The global wisdom structure for precomputed data. \c wisdom.initialized
 * is set to \c false and \c wisdom.flags is set to \c 0U.
 *
 * \author Jens Keiner
 */
#ifdef _OPENMP
  static struct nfsft_wisdom wisdom = {false,0U,-1,-1,0,0,0,0,0,0};
#else
  static struct nfsft_wisdom wisdom = {false,0U,-1,-1,0,0,0,0,0};
#endif

/**
 * Converts coefficients \f$\left(b_k^n\right)_{k=0}^M\f$ with
 * \f$M \in \mathbb{N}_0\f$, \f$-M \le n \le M\f$ from a linear combination
 * of Chebyshev polynomials
 * \f[
 *   f(\cos\vartheta) = \sum_{k=0}^{2\lfloor\frac{M}{2}\rfloor}
 *   a_k (\sin\vartheta)^{n\;\mathrm{mod}\;2} T_k(\cos\vartheta)
 * \f]
 * to coefficients \f$\left(c_k^n\right)_{k=0}^M\f$ matching the representation
 * by complex exponentials
 * \f[
 *   f(\cos\vartheta) = \sum_{k=-M}^{M} c_k \mathrm{e}^{\mathrm{i}k\vartheta}
 * \f]
 * for each order \f$n=-M,\ldots,M\f$.
 *
 * \arg plan The \c nfsft_plan containing the coefficients
 *      \f$\left(b_k^n\right)_{k=0,\ldots,M;n=-M,\ldots,M}\f$
 *
 * \remark The transformation is computed in place.
 *
 * \author Jens Keiner
 */
static inline void c2e(nfsft_plan *plan)
{
  int k;               /**< The degree k                                     */
  int n;               /**< The order k                                      */
  double _Complex last; /**< Stores temporary values                          */
  double _Complex act;  /**< Stores temporary values                          */
  double _Complex *xp;  /**< Auxilliary pointer                               */
  double _Complex *xm;  /**< Auxilliary pointer                               */
  int low;             /**< Lower loop bound                                 */
  int up;              /**< Upper loop bound                                 */
  int lowe;            /**< Lower loop bound for even terms                  */
  int upe;             /**< Upper loop bound for even terms                  */

  /* Set the first row to order to zero since it is unused. */
  memset(plan->f_hat_intern,0U,(2*plan->N+2)*sizeof(double _Complex));

  /* Determine lower and upper bounds for loop processing even terms. */
  lowe = -plan->N + (plan->N%2);
  upe = -lowe;

  /* Process even terms. */
  for (n = lowe; n <= upe; n += 2)
  {
    /* Compute new coefficients \f$\left(c_k^n\right)_{k=-M,\ldots,M}\f$ from
     * old coefficients $\left(b_k^n\right)_{k=0,\ldots,M}$. */
    xm = &(plan->f_hat_intern[NFSFT_INDEX(-1,n,plan)]);
    xp = &(plan->f_hat_intern[NFSFT_INDEX(+1,n,plan)]);
    for(k = 1; k <= plan->N; k++)
    {
      *xp *= 0.5;
      *xm-- = *xp++;
    }
    /* Set the first coefficient in the array corresponding to this order to
     * zero since it is unused. */
    *xm = 0.0;
  }

  /* Determine lower and upper bounds for loop processing odd terms. */
  low = -plan->N + (1-plan->N%2);
  up = -low;

  /* Process odd terms incorporating the additional sine term
   * \f$\sin \vartheta\f$. */
  for (n = low; n <= up; n += 2)
  {
    /* Compute new coefficients \f$\left(c_k^n\right)_{k=-M,\ldots,M}\f$ from
     * old coefficients $\left(b_k^n\right)_{k=0,\ldots,M-1}$ incorporating
     * the additional term \f$\sin \vartheta\f$. */
    plan->f_hat_intern[NFSFT_INDEX(0,n,plan)] *= 2.0;
    xp = &(plan->f_hat_intern[NFSFT_INDEX(-plan->N-1,n,plan)]);
    /* Set the first coefficient in the array corresponding to this order to zero
     * since it is unused. */
    *xp++ = 0.0;
    xm = &(plan->f_hat_intern[NFSFT_INDEX(plan->N,n,plan)]);
    last = *xm;
    *xm = 0.5 * _Complex_I * (0.5*xm[-1]);
    *xp++ = -(*xm--);
    for (k = plan->N-1; k > 0; k--)
    {
      act = *xm;
      *xm = 0.5 * _Complex_I * (0.5*(xm[-1] - last));
      *xp++ = -(*xm--);
      last = act;
    }
    *xm = 0.0;
  }
}

/**
 * Transposed version of the function \ref c2e
 *
 * \arg plan The \c nfsft_plan containing the coefficients
 *      \f$\left(c_k^n\right)_{k=-M,\ldots,M;n=-M,\ldots,M}\f$
 *
 * \remark The transformation is computed in place.
 *
 * \author Jens Keiner
 */
static inline void c2e_transposed(nfsft_plan *plan)
{
  int k;               /**< The degree k                                     */
  int n;               /**< The order k                                      */
  double _Complex last; /**< Stores temporary values                          */
  double _Complex act;  /**< Stores temporary values                          */
  double _Complex *xp;  /**< Auxilliary pointer                               */
  double _Complex *xm;  /**< Auxilliary pointer                               */
  int low;             /**< Lower loop bound                                 */
  int up;              /**< Upper loop bound                                 */
  int lowe;            /**< Lower loop bound for even terms                  */
  int upe;             /**< Upper loop bound for even terms                  */

  /* Determine lower and upper bounds for loop processing even terms. */
  lowe = -plan->N + (plan->N%2);
  upe = -lowe;

  /* Process even terms. */
  for (n = lowe; n <= upe; n += 2)
  {
    /* Compute new coefficients \f$\left(b_k^n\right)_{k=0,\ldots,M}\f$ from
     * old coefficients $\left(c_k^n\right)_{k=-M,\ldots,M}$. */
    xm = &(plan->f_hat[NFSFT_INDEX(-1,n,plan)]);
    xp = &(plan->f_hat[NFSFT_INDEX(+1,n,plan)]);
    for(k = 1; k <= plan->N; k++)
    {
      *xp += *xm--;
      *xp++ *= 0.5;;
    }
  }

  /* Determine lower and upper bounds for loop processing odd terms. */
  low = -plan->N + (1-plan->N%2);
  up = -low;

  /* Process odd terms. */
  for (n = low; n <= up; n += 2)
  {
    /* Compute new coefficients \f$\left(b_k^n\right)_{k=0,\ldots,M-1}\f$ from
     * old coefficients $\left(c_k^n\right)_{k=0,\ldots,M-1}$. */
    xm = &(plan->f_hat[NFSFT_INDEX(-1,n,plan)]);
    xp = &(plan->f_hat[NFSFT_INDEX(+1,n,plan)]);
    for(k = 1; k <= plan->N; k++)
    {
      *xp++ -= *xm--;
    }

    plan->f_hat[NFSFT_INDEX(0,n,plan)] =
      -0.25*_Complex_I*plan->f_hat[NFSFT_INDEX(1,n,plan)];
    last = plan->f_hat[NFSFT_INDEX(1,n,plan)];
    plan->f_hat[NFSFT_INDEX(1,n,plan)] =
      -0.25*_Complex_I*plan->f_hat[NFSFT_INDEX(2,n,plan)];

    xp = &(plan->f_hat[NFSFT_INDEX(2,n,plan)]);
    for (k = 2; k < plan->N; k++)
    {
      act = *xp;
      *xp = -0.25 * _Complex_I * (xp[1] - last);
      xp++;
      last = act;
    }
    *xp = 0.25 * _Complex_I * last;

    plan->f_hat[NFSFT_INDEX(0,n,plan)] *= 2.0;
  }
}

void nfsft_init(nfsft_plan *plan, int N, int M)
{
  /* Call nfsft_init_advanced with default flags. */
  nfsft_init_advanced(plan, N, M, NFSFT_MALLOC_X | NFSFT_MALLOC_F |
    NFSFT_MALLOC_F_HAT);
}

void nfsft_init_advanced(nfsft_plan* plan, int N, int M,
                         unsigned int flags)
{
  /* Call nfsft_init_guru with the flags and default NFFT cut-off. */
  nfsft_init_guru(plan, N, M, flags, PRE_PHI_HUT | PRE_PSI | FFTW_INIT | NFFT_OMP_BLOCKWISE_ADJOINT,
                         NFSFT_DEFAULT_NFFT_CUTOFF);
}

void nfsft_init_guru(nfsft_plan *plan, int N, int M, unsigned int flags,
  unsigned int nfft_flags, int nfft_cutoff)
{
  int *nfft_size; /*< NFFT size                                              */
  int *fftw_size; /*< FFTW size                                              */

  /* Save the flags in the plan. */
  plan->flags = flags;

  /* Save the bandwidth N and the number of samples M in the plan. */
  plan->N = N;
  plan->M_total = M;

  /* M is fixed for FSFT algorithm */
  if (plan->flags & NFSFT_EQUISPACED)
    plan->M_total = (2*plan->N+2)*(plan->N+2);

  /* Calculate the next greater power of two with respect to the bandwidth N
   * and the corresponding exponent. */
  //next_power_of_2_exp_int(plan->N,&plan->NPT,&plan->t);

  /* Save length of array of Fourier coefficients. Owing to the data layout the
   * length is (2N+2)(2N+2) */
  plan->N_total = (2*plan->N+2)*(2*plan->N+2);

  /* Allocate memory for auxiliary array of spherical Fourier coefficients,
   * if necessary. */
  if (plan->flags & NFSFT_PRESERVE_F_HAT)
  {
    plan->f_hat_intern = (double _Complex*) nfft_malloc(plan->N_total*
                                                  sizeof(double _Complex));
  }

  /* Allocate memory for spherical Fourier coefficients, if necessary. */
  if (plan->flags & NFSFT_MALLOC_F_HAT)
  {
    plan->f_hat = (double _Complex*) nfft_malloc(plan->N_total*
                                           sizeof(double _Complex));
  }

  /* Allocate memory for samples, if necessary. */
  if (plan->flags & NFSFT_MALLOC_F)
  {
    plan->f = (double _Complex*) nfft_malloc(plan->M_total*sizeof(double _Complex));
  }

  /* Allocate memory for nodes, if necessary. */
  if (plan->flags & NFSFT_MALLOC_X)
  {
    plan->x = (double*) nfft_malloc(plan->M_total*2*sizeof(double));
    if (plan->flags & NFSFT_EQUISPACED)
      /* Set equispaced nodes. This way also trafo_direct works correctly. */
      for (int i=0; i<2*plan->N+2; i++)
        for (int j=0; j<plan->N+2; j++)
        {
          plan->x[2*(i*(plan->N+1) + j)] = ((double)i-plan->N-1)/(2.0*plan->N+2);
          plan->x[2*(i*(plan->N+1) + j) + 1] = ((double)j)/(2.0*plan->N+2);
        }
  }

  /* Check if fast algorithm is activated. */
  if ((plan->flags & NFSFT_NO_FAST_ALGORITHM) || (plan->flags & NFSFT_EQUISPACED))
  {
  }
  else
  {
      nfft_size = (int*)nfft_malloc(2*sizeof(int));
      fftw_size = (int*)nfft_malloc(2*sizeof(int));

      /** \todo Replace 4*plan->N by next_power_of_2(2*this->n). */
      nfft_size[0] = 2*plan->N+2;
      nfft_size[1] = 2*plan->N+2;
      fftw_size[0] = 4*plan->N;
      fftw_size[1] = 4*plan->N;

      /** \todo NFSFT: Check NFFT flags. */
      nfft_init_guru(&plan->plan_nfft, 2, nfft_size, plan->M_total, fftw_size,
                     nfft_cutoff, nfft_flags,
                     FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

      /* Assign angle array. */
      plan->plan_nfft.x = plan->x;
      /* Assign result array. */
      plan->plan_nfft.f = plan->f;
      /* Assign Fourier coefficients array. */
      plan->plan_nfft.f_hat = plan->f_hat;

      /** \todo Add precomputation if neccessary and possible. */

      /* Precompute. */
      //nfft_precompute_one_psi(&plan->plan_nfft);

      /* Free auxilliary arrays. */
      nfft_free(nfft_size);
      nfft_free(fftw_size);
  }

  plan->mv_trafo = (void (*) (void* ))nfsft_trafo;
  plan->mv_adjoint = (void (*) (void* ))nfsft_adjoint;
}

void nfsft_precompute(int N, double kappa, unsigned int nfsft_flags,
  unsigned int fpt_flags)
{
  int n; /*< The order n                                                     */

  /*  Check if already initialized. */
  if (wisdom.initialized == true)
  {
    return;
  }

#ifdef _OPENMP
  #pragma omp parallel default(shared)
  {
    int nthreads = omp_get_num_threads();
    #pragma omp single
    {
      wisdom.nthreads = nthreads;
    }
  }
#endif

  /* Save the precomputation flags. */
  wisdom.flags = nfsft_flags;

  /* Compute and save N_max = 2^{\ceil{log_2 N}} as next greater
   * power of two with respect to N. */
  X(next_power_of_2_exp_int)(N,&wisdom.N_MAX,&wisdom.T_MAX);

  /* Check, if precomputation for direct algorithms needs to be performed. */
  if (wisdom.flags & NFSFT_NO_DIRECT_ALGORITHM)
  {
    wisdom.alpha = NULL;
    wisdom.beta = NULL;
    wisdom.gamma = NULL;
  }
  else
  {
    /* Allocate memory for three-term recursion coefficients. */
    wisdom.alpha = (double*) nfft_malloc((wisdom.N_MAX+1)*(wisdom.N_MAX+2)*
      sizeof(double));
    wisdom.beta = (double*) nfft_malloc((wisdom.N_MAX+1)*(wisdom.N_MAX+2)*
      sizeof(double));
    wisdom.gamma = (double*) nfft_malloc((wisdom.N_MAX+1)*(wisdom.N_MAX+2)*
      sizeof(double));
    /** \todo Change to functions which compute only for fixed order n. */
    /* Compute three-term recurrence coefficients alpha_k^n, beta_k^n, and
     * gamma_k^n. */
    alpha_al_all(wisdom.alpha,wisdom.N_MAX);
    beta_al_all(wisdom.beta,wisdom.N_MAX);
    gamma_al_all(wisdom.gamma,wisdom.N_MAX);
  }

  /* Check, if precomputation for fast algorithms needs to be performed. */
  if (wisdom.flags & NFSFT_NO_FAST_ALGORITHM)
  {
  }
  else if (wisdom.N_MAX >= NFSFT_BREAK_EVEN)
  {
    /* Precompute data for DPT/FPT. */

    /* Check, if recursion coefficients have already been calculated. */
    if (wisdom.alpha != NULL)
    {
#ifdef _OPENMP
      #pragma omp parallel default(shared) private(n)
      {
        int nthreads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        #pragma omp single
        {
          wisdom.nthreads = nthreads;
          wisdom.set_threads = (fpt_set*) nfft_malloc(nthreads*sizeof(fpt_set));
        }

        if (threadid == 0)
          wisdom.set_threads[threadid] = fpt_init(wisdom.N_MAX+1, wisdom.T_MAX,
            fpt_flags | FPT_AL_SYMMETRY | FPT_PERSISTENT_DATA);
        else
          wisdom.set_threads[threadid] = fpt_init(wisdom.N_MAX+1, wisdom.T_MAX,
            fpt_flags | FPT_AL_SYMMETRY | FPT_PERSISTENT_DATA | FPT_NO_INIT_FPT_DATA);

        #pragma omp barrier

        if (threadid != 0)
          wisdom.set_threads[threadid]->dpt = wisdom.set_threads[0]->dpt;

        /* Perform the first part of precomputation which contains most allocations in one thread */
        #pragma omp master
        {
          for (int n = 0; n <= wisdom.N_MAX; n++)
            fpt_precompute_1(wisdom.set_threads[0],n,n);
        }
        #pragma omp barrier

        #pragma omp for private(n) schedule(dynamic)
        for (n = 0; n <= wisdom.N_MAX; n++)
          fpt_precompute_2(wisdom.set_threads[threadid],n,&wisdom.alpha[ROW(n)],
            &wisdom.beta[ROW(n)], &wisdom.gamma[ROW(n)],n,kappa);
      }
#else
      /* Use the recursion coefficients to precompute FPT data using persistent
       * arrays. */
      wisdom.set = fpt_init(wisdom.N_MAX+1, wisdom.T_MAX,
        fpt_flags | FPT_AL_SYMMETRY | FPT_PERSISTENT_DATA);
      for (n = 0; n <= wisdom.N_MAX; n++)
      {
        /*fprintf(stderr,"%d\n",n);
        fflush(stderr);*/
        /* Precompute data for FPT transformation for order n. */
        fpt_precompute(wisdom.set,n,&wisdom.alpha[ROW(n)],&wisdom.beta[ROW(n)],
          &wisdom.gamma[ROW(n)],n,kappa);
      }
#endif
    }
    else
    {
#ifdef _OPENMP
      #pragma omp parallel default(shared) private(n)
      {
        double *alpha, *beta, *gamma;
        int nthreads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        #pragma omp single
        {
          wisdom.nthreads = nthreads;
          wisdom.set_threads = (fpt_set*) nfft_malloc(nthreads*sizeof(fpt_set));
        }

        alpha = (double*) nfft_malloc((wisdom.N_MAX+2)*sizeof(double));
        beta = (double*) nfft_malloc((wisdom.N_MAX+2)*sizeof(double));
        gamma = (double*) nfft_malloc((wisdom.N_MAX+2)*sizeof(double));

        if (threadid == 0)
          wisdom.set_threads[threadid] = fpt_init(wisdom.N_MAX+1, wisdom.T_MAX,
            fpt_flags | FPT_AL_SYMMETRY);
        else
          wisdom.set_threads[threadid] = fpt_init(wisdom.N_MAX+1, wisdom.T_MAX,
            fpt_flags | FPT_AL_SYMMETRY | FPT_NO_INIT_FPT_DATA);

        #pragma omp barrier

        if (threadid != 0)
          wisdom.set_threads[threadid]->dpt = wisdom.set_threads[0]->dpt;

        /* Perform the first part of precomputation which contains most allocations in one thread */
        #pragma omp master
        {
          for (int n = 0; n <= wisdom.N_MAX; n++)
            fpt_precompute_1(wisdom.set_threads[0],n,n);
        }
        #pragma omp barrier

        #pragma omp for private(n) schedule(dynamic)
        for (n = 0; n <= wisdom.N_MAX; n++)
        {
          alpha_al_row(alpha,wisdom.N_MAX,n);
          beta_al_row(beta,wisdom.N_MAX,n);
          gamma_al_row(gamma,wisdom.N_MAX,n);

          /* Precompute data for FPT transformation for order n. */
          fpt_precompute_2(wisdom.set_threads[threadid],n,alpha,beta,gamma,n,
                         kappa);
        }
        /* Free auxilliary arrays. */
        nfft_free(alpha);
        nfft_free(beta);
        nfft_free(gamma);
      }
#else
    /* Allocate memory for three-term recursion coefficients. */
      double *alpha, *beta, *gamma;
      alpha = (double*) nfft_malloc((wisdom.N_MAX+2)*sizeof(double));
      beta = (double*) nfft_malloc((wisdom.N_MAX+2)*sizeof(double));
      gamma = (double*) nfft_malloc((wisdom.N_MAX+2)*sizeof(double));
      wisdom.set = fpt_init(wisdom.N_MAX+1, wisdom.T_MAX,
        fpt_flags | FPT_AL_SYMMETRY);
      for (n = 0; n <= wisdom.N_MAX; n++)
      {
        /*fprintf(stderr,"%d NO_DIRECT\n",n);
        fflush(stderr);*/
        /* Compute three-term recurrence coefficients alpha_k^n, beta_k^n, and
         * gamma_k^n. */
        alpha_al_row(alpha,wisdom.N_MAX,n);
        beta_al_row(beta,wisdom.N_MAX,n);
        gamma_al_row(gamma,wisdom.N_MAX,n);

        /* Precompute data for FPT transformation for order n. */
        fpt_precompute(wisdom.set,n,alpha,beta,gamma,n,
                       kappa);
      }
      /* Free auxilliary arrays. */
      nfft_free(alpha);
      nfft_free(beta);
      nfft_free(gamma);
#endif
    }
  }

  /* Wisdom has been initialised. */
  wisdom.initialized = true;
}

void nfsft_forget(void)
{
  /* Check if wisdom has been initialised. */
  if (wisdom.initialized == false)
  {
    /* Nothing to do. */
    return;
  }

  /* Check, if precomputation for direct algorithms has been performed. */
  if (wisdom.flags & NFSFT_NO_DIRECT_ALGORITHM)
  {
  }
  else
  {
    /* Free arrays holding three-term recurrence coefficients. */
    nfft_free(wisdom.alpha);
    nfft_free(wisdom.beta);
    nfft_free(wisdom.gamma);
    wisdom.alpha = NULL;
    wisdom.beta = NULL;
    wisdom.gamma = NULL;
  }

  /* Check, if precomputation for fast algorithms has been performed. */
  if (wisdom.flags & NFSFT_NO_FAST_ALGORITHM)
  {
  }
  else if (wisdom.N_MAX >= NFSFT_BREAK_EVEN)
  {
#ifdef _OPENMP
    int k;
    for (k = 0; k < wisdom.nthreads; k++)
      fpt_finalize(wisdom.set_threads[k]);
    nfft_free(wisdom.set_threads);
#else
    /* Free precomputed data for FPT transformation. */
    fpt_finalize(wisdom.set);
#endif
  }

  /* Wisdom is now uninitialised. */
  wisdom.initialized = false;
}


void nfsft_finalize(nfsft_plan *plan)
{
  if (!plan)
    return;

  if (!(plan->flags & NFSFT_NO_FAST_ALGORITHM) && !(plan->flags & NFSFT_EQUISPACED))
  {
    /* Finalise the nfft plan. */
    nfft_finalize(&plan->plan_nfft);
  }

  /* De-allocate memory for auxilliary array of spherical Fourier coefficients,
   * if neccesary. */
  if (plan->flags & NFSFT_PRESERVE_F_HAT)
  {
    nfft_free(plan->f_hat_intern);
  }

  /* De-allocate memory for spherical Fourier coefficients, if necessary. */
  if (plan->flags & NFSFT_MALLOC_F_HAT)
  {
    //fprintf(stderr,"deallocating f_hat\n");
    nfft_free(plan->f_hat);
  }

  /* De-allocate memory for samples, if necessary. */
  if (plan->flags & NFSFT_MALLOC_F)
  {
    //fprintf(stderr,"deallocating f\n");
    nfft_free(plan->f);
  }

  /* De-allocate memory for nodes, if necessary. */
  if (plan->flags & NFSFT_MALLOC_X)
  {
    //fprintf(stderr,"deallocating x\n");
    nfft_free(plan->x);
  }
}

static void nfsft_set_f_nan(nfsft_plan *plan)
{
  int m;
  double nan_value = nan("");
  for (m = 0; m < plan->M_total; m++)
    plan->f[m] = nan_value;
}

void nfsft_trafo_direct(nfsft_plan *plan)
{
  int m;               /*< The node index                                    */
  int k;               /*< The degree k                                      */
  int n;               /*< The order n                                       */
  int n_abs;           /*< The absolute value of the order n, ie n_abs = |n| */
  double *alpha;       /*< Pointer to current three-term recurrence
                           coefficient alpha_k^n for associated Legendre
                           functions P_k^n                                   */
  double *gamma;       /*< Pointer to current three-term recurrence
                           coefficient beta_k^n for associated Legendre
                           functions P_k^n                                   */
  double _Complex *a;   /*< Pointer to auxiliary array for Clenshaw algor.   */
  double stheta;       /*< Current angle theta for Clenshaw algorithm        */
  double sphi;         /*< Current angle phi for Clenshaw algorithm          */

#ifdef MEASURE_TIME
  plan->MEASURE_TIME_t[0] = 0.0;
  plan->MEASURE_TIME_t[1] = 0.0;
  plan->MEASURE_TIME_t[2] = 0.0;
#endif

  if (wisdom.flags & NFSFT_NO_DIRECT_ALGORITHM)
  {
    nfsft_set_f_nan(plan);
    return;
  }

  /* Copy spherical Fourier coefficients, if necessary. */
  if (plan->flags & NFSFT_PRESERVE_F_HAT)
  {
    memcpy(plan->f_hat_intern,plan->f_hat,plan->N_total*
           sizeof(double _Complex));
  }
  else
  {
    plan->f_hat_intern = plan->f_hat;
  }

  /* Check, if we compute with L^2-normalized spherical harmonics. If so,
   * multiply spherical Fourier coefficients with corresponding normalization
   * weight. */
  if (plan->flags & NFSFT_NORMALIZED)
  {
    /* Traverse Fourier coefficients array. */
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k,n) schedule(dynamic)
#endif
    for (k = 0; k <= plan->N; k++)
    {
      for (n = -k; n <= k; n++)
      {
        /* Multiply with normalization weight. */
        plan->f_hat_intern[NFSFT_INDEX(k,n,plan)] *=
          sqrt((2*k+1)/(4.0*KPI));
      }
    }
  }

  /* Distinguish by bandwidth M. */
  if (plan->N == 0)
  {
    /* N = 0 */

    /* Constant function */
    for (m = 0; m < plan->M_total; m++)
    {
      plan->f[m] = plan->f_hat_intern[NFSFT_INDEX(0,0,plan)];
    }
  }
  else
  {
    /* N > 0 */

    /* Evaluate
     *   \sum_{k=0}^N \sum_{n=-k}^k a_k^n P_k^{|n|}(cos theta_m) e^{i n phi_m}
     *   = \sum_{n=-N}^N \sum_{k=|n|}^N a_k^n P_k^{|n|}(cos theta_m)
     *     e^{i n phi_m}.
     */
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(m,stheta,sphi,n,a,n_abs,alpha,gamma,k)
#endif
    for (m = 0; m < plan->M_total; m++)
    {
      /* Scale angle theta from [0,1/2] to [0,pi] and apply cosine. */
      stheta = cos(2.0*KPI*plan->x[2*m+1]);
      /* Scale angle phi from [-1/2,1/2] to [-pi,pi]. */
      sphi = 2.0*KPI*plan->x[2*m];

      /* Initialize result for current node. */
      plan->f[m] = 0.0;

      /* For n = -N,...,N, evaluate
       *   b_n := \sum_{k=|n|}^N a_k^n P_k^{|n|}(cos theta_m)
       * using Clenshaw's algorithm.
       */
      for (n = -plan->N; n <= plan->N; n++)
      {
        /* Get pointer to Fourier coefficients vector for current order n. */
        a = &(plan->f_hat_intern[NFSFT_INDEX(0,n,plan)]);

        /* Take absolute value of n. */
        n_abs = abs(n);

        /* Get pointers to three-term recurrence coefficients arrays. */
        alpha = &(wisdom.alpha[ROW(n_abs)]);
        gamma = &(wisdom.gamma[ROW(n_abs)]);

        if (plan->N > 1024)
        {
          /* Clenshaw's algorithm */
          long double _Complex it2 = a[plan->N];
          long double _Complex it1 = a[plan->N-1];
          for (k = plan->N; k > n_abs + 1; k--)
          {
            long double _Complex temp = a[k-2] + it2 * gamma[k];
            it2 = it1 + it2 * alpha[k] * stheta;
            it1 = temp;
          }

          /* Compute final step if neccesary. */
          if (n_abs < plan->N)
          {
            it2 = it1 + it2 * wisdom.alpha[ROWK(n_abs)+1] * stheta;
          }

          /* Compute final result by multiplying the fixed part
           *   gamma_|n| (1-cos^2(theta))^{|n|/2}
           * for order n and the exponential part
           *   e^{i n phi}
           * and add the result to f_m.
           */
          long double _Complex result = it2 * wisdom.gamma[ROW(n_abs)] *
            powl(1- stheta * stheta, 0.5*n_abs) * cexp(_Complex_I*n*sphi);

          plan->f[m] += result;
        }
        else
        {
          /* Clenshaw's algorithm */
          double _Complex it2 = a[plan->N];
          double _Complex it1 = a[plan->N-1];
          for (k = plan->N; k > n_abs + 1; k--)
          {
            double _Complex temp = a[k-2] + it2 * gamma[k];
            it2 = it1 + it2 * alpha[k] * stheta;
            it1 = temp;
          }

          /* Compute final step if neccesary. */
          if (n_abs < plan->N)
          {
            it2 = it1 + it2 * wisdom.alpha[ROWK(n_abs)+1] * stheta;
          }

          /* Compute final result by multiplying the fixed part
           *   gamma_|n| (1-cos^2(theta))^{|n|/2}
           * for order n and the exponential part
           *   e^{i n phi}
           * and add the result to f_m.
           */
          plan->f[m] += it2 * wisdom.gamma[ROW(n_abs)] *
            pow(1- stheta * stheta, 0.5*n_abs) * cexp(_Complex_I*n*sphi);
        }
      }

      /* Write result f_m for current node to array f. */
//      plan->f[m] = f_m;
    }
  }
}

static void nfsft_set_f_hat_nan(nfsft_plan *plan)
{
  int k, n;
  double nan_value = nan("");
  for (k = 0; k <= plan->N; k++)
    for (n = -k; n <= k; n++)
      plan->f_hat[NFSFT_INDEX(k,n,plan)] = nan_value;
}

void nfsft_adjoint_direct(nfsft_plan *plan)
{
  int m;               /*< The node index                                    */
  int k;               /*< The degree k                                      */
  int n;               /*< The order n                                       */
  int n_abs;           /*< The absolute value of the order n, ie n_abs = |n| */
  double *alpha;       /*< Pointer to current three-term recurrence
                           coefficient alpha_k^n for associated Legendre
                           functions P_k^n                                   */
  double *gamma;       /*< Pointer to current three-term recurrence
                           coefficient beta_k^n for associated Legendre
                           functions P_k^n                                   */

#ifdef MEASURE_TIME
  plan->MEASURE_TIME_t[0] = 0.0;
  plan->MEASURE_TIME_t[1] = 0.0;
  plan->MEASURE_TIME_t[2] = 0.0;
#endif

  if (wisdom.flags & NFSFT_NO_DIRECT_ALGORITHM)
  {
    nfsft_set_f_hat_nan(plan);
    return;
  }

  /* Initialise spherical Fourier coefficients array with zeros. */
  memset(plan->f_hat,0U,plan->N_total*sizeof(double _Complex));

  /* Distinguish by bandwidth N. */
  if (plan->N == 0)
  {
    /* N == 0 */

    /* Constant function */
    for (m = 0; m < plan->M_total; m++)
    {
      plan->f_hat[NFSFT_INDEX(0,0,plan)] += plan->f[m];
    }
  }
  else
  {
    /* N > 0 */

#ifdef _OPENMP
      /* Traverse all orders n. */
      #pragma omp parallel for default(shared) private(n,n_abs,alpha,gamma,m,k) schedule(dynamic)
      for (n = -plan->N; n <= plan->N; n++)
      {
        /* Take absolute value of n. */
        n_abs = abs(n);

        /* Get pointers to three-term recurrence coefficients arrays. */
        alpha = &(wisdom.alpha[ROW(n_abs)]);
        gamma = &(wisdom.gamma[ROW(n_abs)]);

        /* Traverse all nodes. */
        for (m = 0; m < plan->M_total; m++)
        {
          /* Scale angle theta from [0,1/2] to [0,pi] and apply cosine. */
          double stheta = cos(2.0*KPI*plan->x[2*m+1]);
          /* Scale angle phi from [-1/2,1/2] to [-pi,pi]. */
          double sphi = 2.0*KPI*plan->x[2*m];

          /* Transposed Clenshaw algorithm */
          if (plan->N > 1024)
          {
            /* Initial step */
            long double _Complex it1 = plan->f[m] * wisdom.gamma[ROW(n_abs)] *
              powl(1 - stheta * stheta, 0.5*n_abs) * cexpl(-_Complex_I*n*sphi);
            plan->f_hat[NFSFT_INDEX(n_abs,n,plan)] += it1;
            long double _Complex it2 = 0.0;

            if (n_abs < plan->N)
            {
              it2 = it1 * wisdom.alpha[ROWK(n_abs)+1] * stheta;
              plan->f_hat[NFSFT_INDEX(n_abs+1,n,plan)] += it2;
            }

            /* Loop for transposed Clenshaw algorithm */
            for (k = n_abs+2; k <= plan->N; k++)
            {
              long double _Complex temp = it2;
              it2 = alpha[k] * stheta * it2 + gamma[k] * it1;
              it1 = temp;
              plan->f_hat[NFSFT_INDEX(k,n,plan)] += it2;
            }
          }
          else
          {
            /* Initial step */
            double _Complex it1 = plan->f[m] * wisdom.gamma[ROW(n_abs)] *
              pow(1 - stheta * stheta, 0.5*n_abs) * cexp(-_Complex_I*n*sphi);
            plan->f_hat[NFSFT_INDEX(n_abs,n,plan)] += it1;
            double _Complex it2 = 0.0;

            if (n_abs < plan->N)
            {
              it2 = it1 * wisdom.alpha[ROWK(n_abs)+1] * stheta;
              plan->f_hat[NFSFT_INDEX(n_abs+1,n,plan)] += it2;
            }

            /* Loop for transposed Clenshaw algorithm */
            for (k = n_abs+2; k <= plan->N; k++)
            {
              double _Complex temp = it2;
              it2 = alpha[k] * stheta * it2 + gamma[k] * it1;
              it1 = temp;
              plan->f_hat[NFSFT_INDEX(k,n,plan)] += it2;
            }
          }
        }
      }
#else
    /* Traverse all nodes. */
    for (m = 0; m < plan->M_total; m++)
    {
      /* Scale angle theta from [0,1/2] to [0,pi] and apply cosine. */
      double stheta = cos(2.0*KPI*plan->x[2*m+1]);
      /* Scale angle phi from [-1/2,1/2] to [-pi,pi]. */
      double sphi = 2.0*KPI*plan->x[2*m];

      /* Traverse all orders n. */
      for (n = -plan->N; n <= plan->N; n++)
      {
        /* Take absolute value of n. */
        n_abs = abs(n);

        /* Get pointers to three-term recurrence coefficients arrays. */
        alpha = &(wisdom.alpha[ROW(n_abs)]);
        gamma = &(wisdom.gamma[ROW(n_abs)]);

        /* Transposed Clenshaw algorithm */
        if (plan->N > 1024)
        {
          /* Initial step */
          long double _Complex it1 = plan->f[m] * wisdom.gamma[ROW(n_abs)] *
            powl(1 - stheta * stheta, 0.5*n_abs) * cexpl(-_Complex_I*n*sphi);
          plan->f_hat[NFSFT_INDEX(n_abs,n,plan)] += it1;
          long double _Complex it2 = 0.0;

          if (n_abs < plan->N)
          {
            it2 = it1 * wisdom.alpha[ROWK(n_abs)+1] * stheta;
            plan->f_hat[NFSFT_INDEX(n_abs+1,n,plan)] += it2;
          }

          /* Loop for transposed Clenshaw algorithm */
          for (k = n_abs+2; k <= plan->N; k++)
          {
            long double _Complex temp = it2;
            it2 = alpha[k] * stheta * it2 + gamma[k] * it1;
            it1 = temp;
            plan->f_hat[NFSFT_INDEX(k,n,plan)] += it2;
          }
        }
        else
        {
          /* Initial step */
          double _Complex it1 = plan->f[m] * wisdom.gamma[ROW(n_abs)] *
            pow(1 - stheta * stheta, 0.5*n_abs) * cexp(-_Complex_I*n*sphi);
          plan->f_hat[NFSFT_INDEX(n_abs,n,plan)] += it1;
          double _Complex it2 = 0.0;

          if (n_abs < plan->N)
          {
            it2 = it1 * wisdom.alpha[ROWK(n_abs)+1] * stheta;
            plan->f_hat[NFSFT_INDEX(n_abs+1,n,plan)] += it2;
          }

          /* Loop for transposed Clenshaw algorithm */
          for (k = n_abs+2; k <= plan->N; k++)
          {
            double _Complex temp = it2;
            it2 = alpha[k] * stheta * it2 + gamma[k] * it1;
            it1 = temp;
            plan->f_hat[NFSFT_INDEX(k,n,plan)] += it2;
          }
        }
      }
    }
#endif
  }

  /* Check, if we compute with L^2-normalized spherical harmonics. If so,
   * multiply spherical Fourier coefficients with corresponding normalization
   * weight. */
  if (plan->flags & NFSFT_NORMALIZED)
  {
    /* Traverse Fourier coefficients array. */
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k,n) schedule(dynamic)
#endif
    for (k = 0; k <= plan->N; k++)
    {
      for (n = -k; n <= k; n++)
      {
        /* Multiply with normalization weight. */
        plan->f_hat[NFSFT_INDEX(k,n,plan)] *=
          sqrt((2*k+1)/(4.0*KPI));
      }
    }
  }

  /* Set unused coefficients to zero. */
  if (plan->flags & NFSFT_ZERO_F_HAT)
  {
    for (n = -plan->N; n <= plan->N+1; n++)
    {
      memset(&plan->f_hat[NFSFT_INDEX(-plan->N-1,n,plan)],0U,
        (plan->N+1+abs(n))*sizeof(double _Complex));
    }
  }
}

void nfsft_trafo(nfsft_plan *plan)
{
  int k; /*< The degree k                                                    */
  int n; /*< The order n                                                     */
#ifdef MEASURE_TIME
  ticks t0, t1;
#endif
  #ifdef DEBUG
    double t, t_pre, t_nfft, t_fpt, t_c2e, t_norm;
    t_pre = 0.0;
    t_norm = 0.0;
    t_fpt = 0.0;
    t_c2e = 0.0;
    t_nfft = 0.0;
  #endif

#ifdef MEASURE_TIME
  plan->MEASURE_TIME_t[0] = 0.0;
  plan->MEASURE_TIME_t[1] = 0.0;
  plan->MEASURE_TIME_t[2] = 0.0;
#endif

  if ((wisdom.flags & NFSFT_NO_FAST_ALGORITHM) || (plan->flags & NFSFT_NO_FAST_ALGORITHM))
  {
    nfsft_set_f_nan(plan);
    return;
  }

  /* Check, if precomputation was done and that the bandwidth N is not too
   * big.
   */
  if (wisdom.initialized == 0 || plan->N > wisdom.N_MAX)
  {
    nfsft_set_f_nan(plan);
    return;
  }

  /* Check, if slow transformation should be used due to small bandwidth. */
  if (plan->N < NFSFT_BREAK_EVEN)
  {
    /* Use NDSFT. */
    nfsft_trafo_direct(plan);
  }

  /* Check for correct value of the bandwidth N. */
  else if (plan->N <= wisdom.N_MAX)
  {
    /* Copy spherical Fourier coefficients, if necessary. */
    if (plan->flags & NFSFT_PRESERVE_F_HAT)
    {
      memcpy(plan->f_hat_intern,plan->f_hat,plan->N_total*
             sizeof(double _Complex));
    }
    else
    {
      plan->f_hat_intern = plan->f_hat;
    }

    /* Propagate pointer values to the internal NFFT plan to assure
     * consistency. Pointers may have been modified externally.
     */
    if (!(plan->flags & NFSFT_EQUISPACED))
    {
      plan->plan_nfft.x = plan->x;
      plan->plan_nfft.f = plan->f;
      plan->plan_nfft.f_hat = plan->f_hat_intern;
    }

    /* Check, if we compute with L^2-normalized spherical harmonics. If so,
     * multiply spherical Fourier coefficients with corresponding normalization
     * weight. */
    if (plan->flags & NFSFT_NORMALIZED)
    {
      /* Traverse Fourier coefficients array. */
#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k,n) schedule(dynamic)
#endif
      for (k = 0; k <= plan->N; k++)
      {
        for (n = -k; n <= k; n++)
        {
          /* Multiply with normalization weight. */
          plan->f_hat_intern[NFSFT_INDEX(k,n,plan)] *=
            sqrt((2*k+1)/(4.0*KPI));
        }
      }
    }

#ifdef MEASURE_TIME
    t0 = getticks();
#endif
    /* Check, which polynomial transform algorithm should be used. */
    if (plan->flags & NFSFT_USE_DPT)
    {
#ifdef _OPENMP
      n = 0;
      fpt_trafo_direct(wisdom.set_threads[0],abs(n),
        &plan->f_hat_intern[NFSFT_INDEX(abs(n),n,plan)],
        &plan->f_hat_intern[NFSFT_INDEX(0,n,plan)],
        plan->N,0U);

      int n_abs;
      #pragma omp parallel for default(shared) private(n_abs,n) num_threads(wisdom.nthreads) schedule(dynamic)
      for (n_abs = 1; n_abs <= plan->N; n_abs++)
      {
        int threadid = omp_get_thread_num();
        n = -n_abs;
        fpt_trafo_direct(wisdom.set_threads[threadid],abs(n),
          &plan->f_hat_intern[NFSFT_INDEX(abs(n),n,plan)],
          &plan->f_hat_intern[NFSFT_INDEX(0,n,plan)],
          plan->N,0U);
        n = n_abs;
        fpt_trafo_direct(wisdom.set_threads[threadid],abs(n),
          &plan->f_hat_intern[NFSFT_INDEX(abs(n),n,plan)],
          &plan->f_hat_intern[NFSFT_INDEX(0,n,plan)],
          plan->N,0U);
      }
#else
      /* Use direct discrete polynomial transform DPT. */
      for (n = -plan->N; n <= plan->N; n++)
      {
        //fprintf(stderr,"nfsft_trafo: n = %d\n",n);
        //fflush(stderr);
        fpt_trafo_direct(wisdom.set,abs(n),
          &plan->f_hat_intern[NFSFT_INDEX(abs(n),n,plan)],
          &plan->f_hat_intern[NFSFT_INDEX(0,n,plan)],
          plan->N,0U);
      }
#endif
    }
    else
    {
#ifdef _OPENMP
      n = 0;
      fpt_trafo(wisdom.set_threads[0],abs(n),
        &plan->f_hat_intern[NFSFT_INDEX(abs(n),n,plan)],
        &plan->f_hat_intern[NFSFT_INDEX(0,n,plan)],
        plan->N,0U);

      int n_abs;
      #pragma omp parallel for default(shared) private(n_abs,n) num_threads(wisdom.nthreads) schedule(dynamic)
      for (n_abs = 1; n_abs <= plan->N; n_abs++)
      {
        int threadid = omp_get_thread_num();
        n = -n_abs;
        fpt_trafo(wisdom.set_threads[threadid],abs(n),
          &plan->f_hat_intern[NFSFT_INDEX(abs(n),n,plan)],
          &plan->f_hat_intern[NFSFT_INDEX(0,n,plan)],
          plan->N,0U);
        n = n_abs;
        fpt_trafo(wisdom.set_threads[threadid],abs(n),
          &plan->f_hat_intern[NFSFT_INDEX(abs(n),n,plan)],
          &plan->f_hat_intern[NFSFT_INDEX(0,n,plan)],
          plan->N,0U);
      }
#else
      /* Use fast polynomial transform FPT. */
      for (n = -plan->N; n <= plan->N; n++)
      {
        //fprintf(stderr,"nfsft_trafo: n = %d\n",n);
        //fflush(stderr);
        fpt_trafo(wisdom.set,abs(n),
          &plan->f_hat_intern[NFSFT_INDEX(abs(n),n,plan)],
          &plan->f_hat_intern[NFSFT_INDEX(0,n,plan)],
          plan->N,0U);
      }
#endif
    }
#ifdef MEASURE_TIME
    t1 = getticks();
    plan->MEASURE_TIME_t[0] = nfft_elapsed_seconds(t1,t0);
#endif

#ifdef MEASURE_TIME
    t0 = getticks();
#endif
    /* Convert Chebyshev coefficients to Fourier coefficients. */
    c2e(plan);
#ifdef MEASURE_TIME
    t1 = getticks();
    plan->MEASURE_TIME_t[1] = nfft_elapsed_seconds(t1,t0);
#endif

#ifdef MEASURE_TIME
    t0 = getticks();
#endif
    if (plan->flags & NFSFT_EQUISPACED)
    {
      /* Algorithm for equispaced nodes.
       * plan_nfft is not initialized if NFSFT_EQUISPACED is set. */

#ifdef _OPENMP
      INT nthreads = Y(get_num_threads)();
#endif

      int N[2];
      N[0] = 2*plan->N+2;
      N[1] = 2*plan->N+2;
      fftw_plan plan_fftw;

      for (int j=0; j<N[0]; j++)
        for (int k=0; k<N[1]; k++)
          if ((j+k)%2)
            plan->f_hat_intern[j*N[1]+k] *= -1;
//	  f_hat[j*N[1]+k] = plan->f_hat_intern[j*N[1]+k] * CEXP(II*KPI*(j+k));

#ifdef _OPENMP
#pragma omp critical (nfft_omp_critical_fftw_plan)
      {
        FFTW(plan_with_nthreads)(nthreads);
#endif
        plan_fftw = fftw_plan_dft(2, N, plan->f_hat_intern, plan->f_hat_intern, FFTW_FORWARD, FFTW_ESTIMATE);
#ifdef _OPENMP
      }
#endif
      fftw_execute(plan_fftw);
      for (int j=0; j<N[0]; j++)
//	for (int k=j%2-1; k<N[1]; k+=2)
//	    plan->f[j*N[1]+k] *= -1;
        for (int k=N[1]/2; k<N[1]+1; k++)
          plan->f[j*(N[1]/2+1)+(k-N[1]/2)] = plan->f_hat_intern[j*N[1]+k%N[1]] * ((j+k)%2 ? -1 : 1);
//          plan->f[j*N[1]+k] *= CEXP(II*KPI*(j-N[0]/2 + k-N[1]/2));
#ifdef _OPENMP
#pragma omp critical (nfft_omp_critical_fftw_plan)
#endif
      fftw_destroy_plan(plan_fftw);
    }
    /* Check, which nonequispaced discrete Fourier transform algorithm should
     * be used.
     */
    else if (plan->flags & NFSFT_USE_NDFT)
    {
      /* Use NDFT. */
      nfft_trafo_direct(&plan->plan_nfft);
    }
    else
    {
      /* Use NFFT. */
      //fprintf(stderr,"nfsft_adjoint: nfft_trafo\n");
      nfft_trafo_2d(&plan->plan_nfft);
    }
#ifdef MEASURE_TIME
    t1 = getticks();
    plan->MEASURE_TIME_t[2] = nfft_elapsed_seconds(t1,t0);
#endif
  }
}

void nfsft_adjoint(nfsft_plan *plan)
{
  int k; /*< The degree k                                                    */
  int n; /*< The order n                                                     */
#ifdef MEASURE_TIME
  ticks t0, t1;
#endif

#ifdef MEASURE_TIME
  plan->MEASURE_TIME_t[0] = 0.0;
  plan->MEASURE_TIME_t[1] = 0.0;
  plan->MEASURE_TIME_t[2] = 0.0;
#endif

  if ((wisdom.flags & NFSFT_NO_FAST_ALGORITHM) || (plan->flags & NFSFT_NO_FAST_ALGORITHM))
  {
    nfsft_set_f_hat_nan(plan);
    return;
  }

  /* Check, if precomputation was done and that the bandwidth N is not too
   * big.
   */
  if (wisdom.initialized == 0 || plan->N > wisdom.N_MAX)
  {
    nfsft_set_f_hat_nan(plan);
    return;
  }

  /* Check, if slow transformation should be used due to small bandwidth. */
  if (plan->N < NFSFT_BREAK_EVEN)
  {
    /* Use adjoint NDSFT. */
    nfsft_adjoint_direct(plan);
  }
  /* Check for correct value of the bandwidth N. */
  else if (plan->N <= wisdom.N_MAX)
  {
    //fprintf(stderr,"nfsft_adjoint: Starting\n");
    //fflush(stderr);
    /* Propagate pointer values to the internal NFFT plan to assure
     * consistency. Pointers may have been modified externally.
     */
    if (!(plan->flags & NFSFT_EQUISPACED))
    {
      plan->plan_nfft.x = plan->x;
      plan->plan_nfft.f = plan->f;
      plan->plan_nfft.f_hat = plan->f_hat;
    }

#ifdef MEASURE_TIME
    t0 = getticks();
#endif
    if (plan->flags & NFSFT_EQUISPACED)
    {
      /* Algorithm for equispaced nodes.
       * plan_nfft is not initialized if NFSFT_EQUISPACED is set. */
      int N[2];
      N[0] = 2*plan->N+2;
      N[1] = 2*plan->N+2;

      for (int j=0; j<N[0]; j++)
      {
        for (int k=0; k<N[1]/2+1; k++)
          plan->f_hat[j*N[1]+k] = 0;
        for (int k=N[1]/2; k<N[1]+1; k++)
          plan->f_hat[j*N[1]+k%N[1]] = plan->f[j*(N[1]/2+1)+k-N[1]/2] * ((j+k)%2 ? -1 : 1);
      }
      fftw_plan plan_fftw = FFTW(plan_dft)(2, N, plan->f_hat, plan->f_hat, FFTW_BACKWARD, FFTW_ESTIMATE);
      fftw_execute(plan_fftw);
      for (int j=0; j<N[0]; j++)
        for (int k=0; k<N[1]; k++)
          if ((j+k)%2)
            plan->f_hat[j*N[1]+k] *= -1;
      fftw_destroy_plan(plan_fftw);
    }
    /* Check, which adjoint nonequispaced discrete Fourier transform algorithm
     * should be used.
     */
    else if (plan->flags & NFSFT_USE_NDFT)
    {
      //fprintf(stderr,"nfsft_adjoint: Executing nfft_adjoint_direct\n");
      //fflush(stderr);
      /* Use adjoint NDFT. */
      nfft_adjoint_direct(&plan->plan_nfft);
    }
    else
    {
      //fprintf(stderr,"nfsft_adjoint: Executing nfft_adjoint\n");
      //fflush(stderr);
      //fprintf(stderr,"nfsft_adjoint: nfft_adjoint\n");
      /* Use adjoint NFFT. */
      nfft_adjoint_2d(&plan->plan_nfft);
    }
#ifdef MEASURE_TIME
    t1 = getticks();
    plan->MEASURE_TIME_t[2] = nfft_elapsed_seconds(t1,t0);
#endif

    //fprintf(stderr,"nfsft_adjoint: Executing c2e_transposed\n");
    //fflush(stderr);
#ifdef MEASURE_TIME
    t0 = getticks();
#endif
    /* Convert Fourier coefficients to Chebyshev coefficients. */
    c2e_transposed(plan);
#ifdef MEASURE_TIME
    t1 = getticks();
    plan->MEASURE_TIME_t[1] = nfft_elapsed_seconds(t1,t0);
#endif

#ifdef MEASURE_TIME
    t0 = getticks();
#endif
    /* Check, which transposed polynomial transform algorithm should be used */
    if (plan->flags & NFSFT_USE_DPT)
    {
#ifdef _OPENMP
      n = 0;
      fpt_transposed_direct(wisdom.set_threads[0],abs(n),
        &plan->f_hat[NFSFT_INDEX(abs(n),n,plan)],
        &plan->f_hat[NFSFT_INDEX(0,n,plan)],
        plan->N,0U);

      int n_abs;
      #pragma omp parallel for default(shared) private(n_abs,n) num_threads(wisdom.nthreads) schedule(dynamic)
      for (n_abs = 1; n_abs <= plan->N; n_abs++)
      {
        int threadid = omp_get_thread_num();
        n = -n_abs;
        fpt_transposed_direct(wisdom.set_threads[threadid],abs(n),
          &plan->f_hat[NFSFT_INDEX(abs(n),n,plan)],
          &plan->f_hat[NFSFT_INDEX(0,n,plan)],
          plan->N,0U);
        n = n_abs;
        fpt_transposed_direct(wisdom.set_threads[threadid],abs(n),
          &plan->f_hat[NFSFT_INDEX(abs(n),n,plan)],
          &plan->f_hat[NFSFT_INDEX(0,n,plan)],
          plan->N,0U);
      }
#else
      /* Use transposed DPT. */
      for (n = -plan->N; n <= plan->N; n++)
      {
        //fprintf(stderr,"nfsft_adjoint: Executing dpt_transposed\n");
        //fflush(stderr);
        fpt_transposed_direct(wisdom.set,abs(n),
          &plan->f_hat[NFSFT_INDEX(abs(n),n,plan)],
          &plan->f_hat[NFSFT_INDEX(0,n,plan)],
          plan->N,0U);
      }
#endif
    }
    else
    {
#ifdef _OPENMP
      n = 0;
      fpt_transposed(wisdom.set_threads[0],abs(n),
        &plan->f_hat[NFSFT_INDEX(abs(n),n,plan)],
        &plan->f_hat[NFSFT_INDEX(0,n,plan)],
        plan->N,0U);

      int n_abs;
       #pragma omp parallel for default(shared) private(n_abs,n) num_threads(wisdom.nthreads) schedule(dynamic)
       for (n_abs = 1; n_abs <= plan->N; n_abs++)
       {
         int threadid = omp_get_thread_num();
         n = -n_abs;
         fpt_transposed(wisdom.set_threads[threadid],abs(n),
           &plan->f_hat[NFSFT_INDEX(abs(n),n,plan)],
           &plan->f_hat[NFSFT_INDEX(0,n,plan)],
           plan->N,0U);
         n = n_abs;
         fpt_transposed(wisdom.set_threads[threadid],abs(n),
           &plan->f_hat[NFSFT_INDEX(abs(n),n,plan)],
           &plan->f_hat[NFSFT_INDEX(0,n,plan)],
           plan->N,0U);
       }
#else
      //fprintf(stderr,"nfsft_adjoint: fpt_transposed\n");
      /* Use transposed FPT. */
      for (n = -plan->N; n <= plan->N; n++)
      {
        //fprintf(stderr,"nfsft_adjoint: Executing fpt_transposed\n");
        //fflush(stderr);
        fpt_transposed(wisdom.set,abs(n),
          &plan->f_hat[NFSFT_INDEX(abs(n),n,plan)],
          &plan->f_hat[NFSFT_INDEX(0,n,plan)],
          plan->N,0U);
      }
#endif
    }
#ifdef MEASURE_TIME
    t1 = getticks();
    plan->MEASURE_TIME_t[0] = nfft_elapsed_seconds(t1,t0);
#endif

    /* Check, if we compute with L^2-normalized spherical harmonics. If so,
     * multiply spherical Fourier coefficients with corresponding normalization
     * weight. */
    if (plan->flags & NFSFT_NORMALIZED)
    {
      //fprintf(stderr,"nfsft_adjoint: Normalizing\n");
      //fflush(stderr);
      /* Traverse Fourier coefficients array. */
#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k,n) schedule(dynamic)
#endif
      for (k = 0; k <= plan->N; k++)
      {
        for (n = -k; n <= k; n++)
        {
          /* Multiply with normalization weight. */
          plan->f_hat[NFSFT_INDEX(k,n,plan)] *=
            sqrt((2*k+1)/(4.0*KPI));
        }
      }
    }

    /* Set unused coefficients to zero. */
    if (plan->flags & NFSFT_ZERO_F_HAT)
    {
      //fprintf(stderr,"nfsft_adjoint: Setting to zero\n");
      //fflush(stderr);
      for (n = -plan->N; n <= plan->N+1; n++)
      {
        memset(&plan->f_hat[NFSFT_INDEX(-plan->N-1,n,plan)],0U,
          (plan->N+1+abs(n))*sizeof(double _Complex));
      }
    }
    //fprintf(stderr,"nfsft_adjoint: Finished\n");
    //fflush(stderr);
  }
}

void nfsft_precompute_x(nfsft_plan *plan)
{
  if ((plan->flags & NFSFT_NO_FAST_ALGORITHM) || (plan->flags & NFSFT_EQUISPACED))
    return;

  /* Pass angle array to NFFT plan. */
  plan->plan_nfft.x = plan->x;

  /* Precompute. */
  if (plan->plan_nfft.flags & PRE_ONE_PSI)
    nfft_precompute_one_psi(&plan->plan_nfft);
}
/* \} */
