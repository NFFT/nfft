/* $Id$
 *
 * fastsumS2 - Fast summation of spherical radial functions
 *
 * Copyright (C) 2005 Jens Keiner
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/* Include standard C headers. */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

/* Include NFFT3 library header. */
#include "nfft3.h"

/* Include NFFT 3 utilities headers. */
#include "util.h"

#include "../3rdparty/gsl/specfunc/gsl_sf_bessel.h"


/** The Fourier-Legendre coefficients of the Abel-Poisson kernel */
#define SYMBOL_ABEL_POISSON(k,h) (pow(h,k))

/** The Fourier-Legendre coefficients of the singularity kernel */
#define SYMBOL_SINGULARITY(k,h) ((2.0/(2*k+1))*pow(h,k))


/* Flags for the different kernel types */

/** Abel-Poisson kernel */
#define KT_ABEL_POISSON (0)
/** Singularity kernel */
#define KT_SINGULARITY  (1)
/** Locally supported kernel */
#define KT_LOC_SUPP     (2)
/** Gaussian kernel */
#define KT_GAUSSIAN     (3)

/** Enumerations for parameter values */
enum pvalue {NO = 0, YES = 1, BOTH = 2};

/**
 * Computes the \f$\mathbb{R}^3\f$ standard inner product between two vectors
 * given in spherical coordinates.
 *
 * \arg phi1   The angle \f$\varphi_1 \in [-\pi,\pi)\f$ of the first vector
 * \arg theta1 The angle \f$\vartheta_1 \in [0,\pi]\f$ of the first vector
 * \arg phi2   The angle \f$\varphi_" \in [-\pi,\pi)\f$ of the second vector
 * \arg theta2 The angle \f$\vartheta_" \in [0,\pi]\f$ of the second vector
 *
 * \return The inner product \f$\cos\vartheta_1\cos\vartheta_2 +
 *   \sin\vartheta_1\sin(\vartheta_2\cos(\varphi_1-\varphi_2)\f$
 */
double innerProduct(const double phi1, const double theta1, const double phi2,
  const double theta2)
{
  return cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2);
}

/**
 * Evaluates the Poisson kernel \f$Q_h: [-1,1] \rightarrow \mathbb{R}\f$ at a
 * node \f$x \in [-1,1]\f$.
 *
 * \arg x The node \f$x \in [-1,1]\f$
 * \arg h The parameter \f$h \in (0,1)\f$
 *
 * \return The value of the Poisson kernel \f$Q_h(x)\f$ at the node \f$x\f$
 */
double poissonKernel(const double x, const double h)
{
  return (1.0/(4*PI))*(1-h*h)/pow(sqrt(1-2*h*x+h*h),3);
}

/**
 * Evaluates the singularity kernel \f$S_h: [-1,1] \rightarrow \mathbb{R}\f$ at
 * a node \f$x \in [-1,1]\f$.
 *
 * \arg x The node \f$x \in [-1,1]\f$
 * \arg h The parameter \f$h \in (0,1)\f$
 *
 * \return The value of the Poisson kernel \f$S_h(x)\f$ at the node \f$x\f$
 */
double singularityKernel(const double x, const double h)
{
  return (1.0/(2*PI))/sqrt(1-2*h*x+h*h);
}

/**
 * Evaluates the locally supported kernel \f$L_{h,\lambda}: [-1,1] \rightarrow
 * \mathbb{R}\f$ at a node \f$x \in [-1,1]\f$.
 *
 * \arg x The node \f$x \in [-1,1]\f$
 * \arg h The parameter \f$h \in (0,1)\f$
 * \arg lambda The parameter \f$lambda \in \mathbb{N}_0\f$
 *
 * \return The value of the locally supported kernel \f$L_{h,\lambda}(x)\f$ at
 *   the node \f$x\f$
 */
double locallySupportedKernel(const double x, const double h,
  const double lambda)
{
  return (x<=h)?(0.0):(pow((x-h),lambda));
}

  /**
   * Evaluates the spherical Gaussian kernel \f$G_\sigma: [-1,1] \rightarrow
   * \mathbb{R}\f$ at a node \f$x \in [-1,1]\f$.
   *
   * \arg x The node \f$x \in [-1,1]\f$
   * \arg h The parameter \f$\sigma \in \mathbb{R}_+\f$
   *
   * \return The value of the pherical Gaussian kernel \f$G_\sigma(x)\f$ at the
   *   node \f$x\f$
   */
  double gaussianKernel(const double x, const double sigma)
  {
     return exp(2.0*sigma*(x-1));
  }

/**
 * The main program.
 *
 * \param argc The number of arguments
 * \param argv An array containing the arguments as C-strings
 *
 * \return Exit code
 */
int main (int argc, char **argv)
{
  double **p;                  /**< The array containing the parameter sets        *
                                    for the kernels.                               */
  int *m;                      /**< The array containing the cut-off degrees       *
                                    \f$M\f$.                                       */
  int **ld;                    /**< The array containing the numbers of            *
                                    source and target nodes, \f$L \in              *
                                    \mathbb{N}\f$ and \f$D \in \mathbb{N}\f$,      *
                                    respectively.                                  */
  int ip;                      /**< Index variable for \code p                     */
  int im;                      /**< Index variable for \code m                     */
  int ild;                     /**< Index variable for \code l                     */
  int ipp;                     /**< Index for kernel parameters                    */
  int ip_max;                  /**< The maximum index for \code p                  */
  int im_max;                  /**< The maximum index for \code m                  */
  int ild_max;                 /**< The maximum index for \code l                  */
  int ipp_max;                 /**< The maximum index for \code ipp                */
  int tc_max;                  /**< The number of testcases                        */
  int m_max;                   /**< The maximum cut-off degree \f$M\f$ for the     *
                                    current dataset                                */
  int l_max;                   /**< The maximum number of source nodes             *
                                    \f$L\f$ for the current dataset                */
  int d_max;                   /**< The maximum number of target nodes \f$D\f$     *
                                    for the current dataset                        */
  long ld_max_prec;            /**< The maximum number of source and target        *
                                    nodes for precomputation multiplied            */
  long l_max_prec;             /**< The maximum number of source nodes for         *
                                    precomputation                                 */
  int tc;                      /**< Index variable for testcases                   */
  int kt;                      /**< The kernel type                                */
  int cutoff;                  /**< The current NFFT cut-off parameter             */
  double threshold;            /**< The current NFSFT threshold parameter          */
  int n_max;                   /**< Next greater power of two with respect to      *
                                    m_max                                          */
  double t_d;                  /**< Time for direct algorithm in seconds           */
  double t_dp;                 /**< Time for direct algorithm with                 *
                                    precomputation in seconds                      */
  double t_fd;                 /**< Time for fast direct algorithm in seconds      */
  double t_f;                  /**< Time for fast algorithm in seconds             */
  double nfactor;              /**<                                                */
  double temp;                 /**<                                                */
  double err_f;                /**< Error \f$E_\infty\f$ for fast algorithm        */
  double err_fd;               /**< Error \f$E_\infty\f$ for fast direct           *
                                    algorithm                                      */
  double t;                    /**<                                                */
  int precompute = NO;         /**<                                                */
  complex *ptr;                /**<                                                */
  double* steed;               /**<                                                */
  double* steed2;               /**<                                                */
  complex *b;                  /**< The weights \f$\left(b_l\right)_{l=0}          *
                                    ^{L-1}\f$                                      */
  complex *f_hat;              /**< The spherical Fourier coefficients             */
  complex *a;                  /**< The Fourier-Legendre coefficients              */
  double *xi;                  /**< Target nodes                                   */
  double *eta;                 /**< Source nodes                                   */
  complex *f_m;                /**< Approximate function values                    */
  complex *f;                  /**< Exact function values                          */
  complex *prec;               /**<                                                */
  nfsft_plan plan;             /**< NFSFT plan                                     */
  nfsft_plan plan_adjoint;     /**< adjoint NFSFT plan                             */
  int i;                       /**                                                 */
  int j;                       /**                                                 */
  int k;                       /**                                                 */
  int n;                       /**                                                 */
  int d;                       /**                                                 */
  int l;                       /**                                                 */
  int use_nfsft;               /**                                                 */
  int use_nfft;                /**                                                 */
  int use_fpt;                 /**                                                 */
  int nsymbols;                /**                                                 */
  long index;                  /**                                                 */
  int rinc;                    /**                                                 */
  FILE *file_gaussian;         /**                                                 */
  char filename_tex[100];      /**                                                 */
  char filename_dat[100];      /**                                                 */
  char filename_gaussian[100]; /**                                                 */
  double constant;             /**                                                 */

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
    p = (double**) malloc(ip_max*sizeof(double*));

    /* We now read in the parameter sets. */

    /* Read number of parameters. */
    fscanf(stdin,"parameters=%d\n",&ipp_max);
    fprintf(stdout,"%d\n",ipp_max);

    for (ip = 0; ip < ip_max; ip++)
    {
      /* Allocate memory for the parameters. */
      p[ip] = (double*) malloc(ipp_max*sizeof(double));

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
    m = (int*) malloc(im_max*sizeof(int));

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
    ld = (int**) malloc(ild_max*sizeof(int*));

    /* Read the run specification. */
    for (ild = 0; ild < ild_max; ild++)
    {
      /* Allocate memory for the run parameters. */
      ld[ild] = (int*) malloc(5*sizeof(int));

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
    b = (complex*) malloc(l_max*sizeof(complex));
    eta = (double*) malloc(2*l_max*sizeof(double));
    f_hat = (complex*) malloc(NFSFT_F_HAT_SIZE(m_max)*sizeof(complex));
    a = (complex*) malloc((m_max+1)*sizeof(complex));
    xi = (double*) malloc(2*d_max*sizeof(double));
    f_m = (complex*) malloc(d_max*sizeof(complex));
    f = (complex*) malloc(d_max*sizeof(complex));

    /* Allocate memory for precomputed data. */
    if (precompute == YES)
    {
      prec = (complex*) malloc(ld_max_prec*sizeof(complex));
    }

    /* Generate random source nodes and weights. */
    for (l = 0; l < l_max; l++)
    {
      b[l] = drand48() - 0.5;
      eta[2*l] = drand48() - 0.5;
      eta[2*l+1] = 0.5*drand48();
    }

    /* Generate random target nodes. */
    for (d = 0; d < d_max; d++)
    {
      xi[2*d] = drand48() - 0.5;
      xi[2*d+1] = 0.5*drand48();
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
          {
            a[k] = SYMBOL_ABEL_POISSON(k,p[ip][0]);
          }
          break;

        case KT_SINGULARITY:
          /* Compute Fourier-Legendre coefficients for the singularity
           * kernel. */
          for (k = 0; k <= m_max; k++)
          {
            a[k] = SYMBOL_SINGULARITY(k,p[ip][0]);
          }
          break;

        case KT_LOC_SUPP:
          /* Compute Fourier-Legendre coefficients for the locally supported
           * kernel. */
          for (k = 0; k <= m_max; k++)
          {
            /* First case k = 0 for initialization of three-term recurrence. */
            if (k == 0)
            {
              a[k] = 1.0;
            }
            /* Second case k = 1 for initialization of three-term recurrence. */
            else if (k == 1)
            {
              a[k] = ((p[ip][1]+1+p[ip][0])/(p[ip][1]+2.0))*a[k-1];
            }
            /* Apply three-term recurrence. */
            else
            {
              a[k] = (1.0/(k+p[ip][1]+1))*((2*k-1)*p[ip][0]*a[k-1] -
                (k-p[ip][1]-2)*a[k-2]);
            }
          }
          break;

          case KT_GAUSSIAN:
            /* Compute Fourier-Legendre coefficients for the locally supported
             * kernel. */
            steed = (double*) malloc((m_max+1)*sizeof(double));
            steed2 = (double*) malloc((m_max+1)*sizeof(double));
            gsl_sf_bessel_il_scaled_array(m_max,2.0*p[ip][0],steed);
            for (k = 0; k <= m_max; k++)
            {
              steed[k] = 4.0*PI;
              a[k] = steed[k];
            }
            for (k = 0; k <= m_max; k++)
            {
              steed[k] *= 4.0*PI;
              a[k] = steed[k];
            }

            free(steed);
            break;
      }

      /* Normalize Fourier-Legendre coefficients. */
      for (k = 0; k <= m_max; k++)
      {
        a[k] *= (2*k+1)/(4*PI);
      }

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
                temp = innerProduct(2*PI*eta[2*l],2*PI*eta[2*l+1],
                  2*PI*xi[2*d],2*PI*xi[2*d+1]);

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
            t = second();

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
                constant = ((p[ip][1]+1)/(2*PI*pow(1-p[ip][0],p[ip][1]+1)));

                /* Process all target nodes. */
                for (d = 0; d < ld[ild][1]; d++)
                {
                  /* Initialize function value. */
                  f[d] = 0.0;

                  /* Process all source nodes. */
                  for (l = 0; l < ld[ild][0]; l++)
                  {
                    f[d] += b[l]*(*ptr++);
                  }

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
                  {
                    f[d] += b[l]*(*ptr++);
                  }

                  /* Proceed to next row. */
                  ptr += rinc;
                }
              }
            }

            /* Calculate the time needed. */
            t_dp = second() - t;

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
          t = second();

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
                    temp = innerProduct(2*PI*eta[2*l],2*PI*eta[2*l+1],
                      2*PI*xi[2*d],2*PI*xi[2*d+1]);

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
                    temp = innerProduct(2*PI*eta[2*l],2*PI*eta[2*l+1],
                      2*PI*xi[2*d],2*PI*xi[2*d+1]);

                    /* Evaluate the Poisson kernel for the current value and add
                     * to the result. */
                    f[d] += b[l]*singularityKernel(temp,p[ip][0]);
                  }
                }
                break;

              case KT_LOC_SUPP:
                /* Calculate the multiplicative constant. */
                constant = ((p[ip][1]+1)/(2*PI*pow(1-p[ip][0],p[ip][1]+1)));

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
                    temp = innerProduct(2*PI*eta[2*l],2*PI*eta[2*l+1],
                      2*PI*xi[2*d],2*PI*xi[2*d+1]);

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
                      temp = innerProduct(2*PI*eta[2*l],2*PI*eta[2*l+1],
                        2*PI*xi[2*d],2*PI*xi[2*d+1]);
                      /* Evaluate the Poisson kernel for the current value and add
                       * to the result. */
                      f[d] += b[l]*gaussianKernel(temp,p[ip][0]);
                    }
                  }
                  break;
            }
          }

          /* Calculate and add the time needed. */
          t_d = second() - t;
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
            PRE_PHI_HUT | PRE_PSI | FFTW_INIT |
            FFT_OUT_OF_PLACE, cutoff);
          nfsft_init_guru(&plan,m[im],ld[ild][1],
            ((use_nfft!=0)?(0U):(NFSFT_USE_NDFT)) |
            ((use_fpt!=0)?(0U):(NFSFT_USE_DPT)),
            PRE_PHI_HUT | PRE_PSI | FFTW_INIT |
            FFT_OUT_OF_PLACE,
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
            t = second();

            /* Cycle through all runs. */
            for (i = 0; i < ld[ild][4]; i++)
            {

              /* Execute adjoint direct NDSFT transformation. */
              ndsft_adjoint(&plan_adjoint);

              /* Multiplication with the Fourier-Legendre coefficients. */
              for (k = 0; k <= m[im]; k++)
              {
                for (n = -k; n <= k; n++)
                {
                  f_hat[NFSFT_INDEX(k,n,&plan_adjoint)] *= a[k];
                }
              }

              /* Execute direct NDSFT transformation. */
              ndsft_trafo(&plan);

            }

            /* Calculate and add the time needed. */
            t_fd = second() - t;

            /* Calculate average time needed. */
            t_fd = t_fd/((double)ld[ild][4]);

            /* Check if error E_infty should be computed. */
            if (ld[ild][2] != NO)
            {
              /* Compute the error E_infinity. */
              err_fd = error_l_infty_1_complex(f, f_m, ld[ild][1], b,
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
          t = second();

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
              ndsft_adjoint(&plan_adjoint);
            }

            /* Multiplication with the Fourier-Legendre coefficients. */
            for (k = 0; k <= m[im]; k++)
            {
              for (n = -k; n <= k; n++)
              {
                f_hat[NFSFT_INDEX(k,n,&plan_adjoint)] *= a[k];
               }
             }

            /* Check if the fast NFSFT algorithm shall also be tested. */
            if (use_nfsft != NO)
            {
              /* Execute the NFSFT transformation. */
              nfsft_trafo(&plan);
            }
            else
            {
              /* Execute the NDSFT transformation. */
              ndsft_trafo(&plan);
            }

            /*for (d = 0; d < ld[ild][1]; d++)
            {
              fprintf(stderr,"f_ref[%d] = %le + I*%le, f[%d] = %le + I*%le\n",
                d,creal(f[d]),cimag(f[d]),d,creal(f_m[d]),cimag(f_m[d]));
            }*/
          }

          /* Check if the fast NFSFT algorithm has been used. */
          if (use_nfsft != NO)
          {
            /* Calculate and add the time needed. */
            t_f = second() - t;
          }
          else
          {
            /* Calculate and add the time needed. */
            t_fd = second() - t;
          }

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
              err_f = error_l_infty_1_complex(f, f_m, ld[ild][1], b,
                ld[ild][0]);
            }
            else
            {
              /* Compute the error E_infinity. */
              err_fd = error_l_infty_1_complex(f, f_m, ld[ild][1], b,
                ld[ild][0]);
            }
          }

          /* Print out the error measurements. */
          fprintf(stdout,"%e\n%e\n%e\n%e\n%e\n%e\n\n",t_d,t_dp,t_fd,t_f,err_fd,
            err_f);
          fprintf(stderr,"%d: %e\t%e\t%e\t%e\t%e\t%e\n",m[im],t_d,t_dp,t_fd,t_f,
            err_fd,err_f);

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
      free(prec);
    }
    /* Free data arrays. */
    free(f);
    free(f_m);
    free(xi);
    free(eta);
    free(a);
    free(f_hat);
    free(b);

    /* Free memory for node sets. */
    for (ild = 0; ild < ild_max; ild++)
    {
      free(ld[ild]);
    }
    free(ld);

    /* Free memory for cut-off bandwidths. */
    free(m);

    /* Free memory for parameter sets. */
    for (ip = 0; ip < ip_max; ip++)
    {
      free(p[ip]);
    }
    free(p);
  } /* for (tc = 0; tc < tc_max; tc++) - Process each testcase. */

  /* Return exit code for successful run. */
  return EXIT_SUCCESS;
}
