/*
   FastsumS2 - Fast summation of spherical radial functions

   Copyright (C) 2005 Jens Keiner

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

#include <termios.h>
#include <grp.h>
#include <pwd.h>
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

/* Include GSL header. */
//#include <gsl/gsl_sf_bessel.h>

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

#define NO   (0)
#define YES  (1)
#define BOTH (2)

/**
 * Computes the \f$\mathbb{S}^2\f$ inner product between two vectors given in
 * spherical coordinates.
 *
 * \arg phi1   The angle \f$\varphi_1 \in [-\pi,\pi)\f$ of the first vector
 * \arg theta1 The angle \f$\vartheta_1 \in [0,\pi]\f$ of the first vector
 * \arg phi2   The angle \f$\varphi_" \in [-\pi,\pi)\f$ of the second vector
 * \arg theta2 The angle \f$\vartheta_" \in [0,\pi]\f$ of the second vector
 *
 * \return The inner product \f$\cos\vartheta_1\cos\vartheta_2 +
 *         \sin\vartheta_1\sin(\vartheta_2\cos(\varphi_1-\varphi_2)\f$
 */
inline double innerProduct(const double phi1, const double theta1,
                           const double phi2, const double theta2)
{
  return cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2);
}

/**
 * Evaluates the Poisson kernel \f$Q_h\f$ at a node \f$x \in [-1,1]\f$.
 *
 * \arg x The node \f$x \in [-1,1]\f$
 * \arg h The parameter \f$h \in (0,1)\f$
 *
 * \return The value of the Poisson kernel \f$Q_h\f$ at the node \f$x\f$
 */
inline double poissonKernel(const double x, const double h)
{
   return (1.0/(4*PI))*(1-h*h)/pow(sqrt(1-2*h*x+h*h),3);
}

/**
 * Evaluates the singularity kernel \f$S_h\f$ at a node \f$x \in [-1,1]\f$.
 *
 * \arg x The node \f$x \in [-1,1]\f$
 * \arg h The parameter \f$h \in (0,1)\f$
 *
 * \return The value of the Poisson kernel \f$S_h\f$ at the node \f$x\f$
 */
inline double singularityKernel(const double x, const double h)
{
  return (1.0/(2*PI))/sqrt(1-2*h*x+h*h);
}

/**
 * Evaluates the locally supported kernel \f$L_{h,\lambda}\f$ at a node \f$x \in
 * [-1,1]\f$.
 *
 * \arg x The node \f$x \in [-1,1]\f$
 * \arg h The parameter \f$h \in (0,1)\f$
 * \arg lambda The parameter \f$lambda \in \mathbb{N}_0\f$
 *
 * \return The value of the locally supported kernel \f$L_{h,\lambda}\f$ at the node
 * \f$x\f$
 */
inline double locSuppKernel(const double x, const double h, const double lambda)
{
   return (x<=h)?(0.0):(pow((x-h),lambda));
}

/**
 * Evaluates the spherical Gaussian kernel \f$G_\sigma\f$ at a node \f$x \in
 * [-1,1]\f$.
 *
 * \arg x The node \f$x \in [-1,1]\f$
 * \arg h The parameter \f$\sigma \in \mathbb{R}_+\f$
 *
 * \return The value of the pherical Gaussian kernel \f$G_\sigma\f$ at the node
 * \f$x\f$
 */
inline double gaussianKernel(const double x, const double sigma)
{
   return exp(2*sigma*(x-1));
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
  double **p;                  /**< The array containing the parameter sets for
                                    the kernels.                              */
  int *m;                      /**< The array containing the cut-off degrees
                                    \f$M\f$.                                  */
  int **ld;                    /**< The array containing the numbers of source
                                    and target nodes, \f$L\f$ and \f$D\f$,
                                    respectively.                             */
  int ip;                      /**< Index variable for p                      */
  int im;                      /**< Index variable for m                      */
  int ild;                     /**< Index variable for l                      */
  int ip_max;                  /**< The maximum index for p                   */
  int im_max;                  /**< The maximum index for m                   */
  int ild_max;                 /**< The maximum index for l                   */
  int tc_max;                  /**< The number of testcases                   */
  int m_max;                   /**< The maximum cut-off degree \f$M\f$ for the
                                    current dataset                           */
  int l_max;                   /**< The maximum number of source nodes \f$L\f$
                                    for the current dataset                   */
  int d_max;                   /**< The maximum number of target nodes \f$D\f$
                                    for the current dataset                   */
  int ipp_max;                 /**< */
  int ipp;                     /**< */
  long ld_max_prec;            /**< */
  long l_max_prec;             /**< */
  int tc;                      /**< Index variable for testcases.             */
  int kt;                      /**< The kernel type                           */
  int cutoff;                  /**< */
  double threshold;            /**< */
  int n_max;                   /**< Next greater power of two with respect to
                                    m_max                                     */
  double t_d;                  /**< */
  double t_dp;                 /**< */
  double t_fd;                 /**< */
  double t_f;                  /**< */
  double nfactor;              /**< */
  double temp;                 /**< */
  double err_f;                /**< */
  double err_fd;               /**< */
  double t;                    /**< */
  int precompute = NO;         /**< */
  complex *ptr;                /**< */
  double* steed;               /**< */
  complex *b;                  /**< The weights \f$\left(b_l\right)_{l=0}
                                    ^{L-1}\f$                                 */
  complex *f_hat;              /**< The spherical Fourier coefficients        */
  complex *a;                  /**< The Fourier-Legendre coefficients         */
  double *xi;                  /**< Target nodes                              */
  double *eta;                 /**< Source nodes                              */
  complex *f_m;                /**< Approximate function values               */
  complex *f;                  /**< Exact function values                     */
  complex *prec;               /**< */
  nfsft_plan plan;             /**< NFSFT plan                                */
  nfsft_plan plan_adjoint;     /**< adjoint NFSFT plan                        */
  int i;                       /** */
  int j;                       /** */
  int k;                       /** */
  int n;                       /** */
  int d;                       /** */
  int l;                       /** */
  int use_nfsft;               /** */
  int use_nfft;                /** */
  int use_fpt;                 /** */
  int nsymbols;                /** */
  long index;                  /** */
  int rinc;                    /** */
  //FILE *file_tex;              /** */
  //FILE *file_dat;              /** */
  FILE *file_gaussian;         /** */
  char filename_tex[100];      /** */
  char filename_dat[100];      /** */
  char filename_gaussian[100]; /** */
  double constant;             /** */

  /* Read number of testcases. */
  fscanf(stdin,"testcases=%d\n",&tc_max);
  fprintf(stdout,"%d\n",tc_max);
  /* Printf out number of testcases. */
  //fprintf(stdout,"Number of testcases: %d\n\n",tc_max);

  /* Process each testcase. */
  for (tc = 0; tc < tc_max; tc++)
  {
    //fprintf(stdout,"Testcase %d:\n",tc);

    //sprintf(filename_dat,"testcase%d.dat",tc);
    //file_dat = fopen(filename_dat,"w");

    /* Check if fast transform shall be used. */
    fscanf(stdin,"nfsft=%d\n",&use_nfsft);
    fprintf(stdout,"%d\n",use_nfsft);
    if (use_nfsft != NO)
      {
      //fprintf(stdout,"  NFSFT = yes\n");
      /* Check if the NFFT shall be used. */
      fscanf(stdin,"nfft=%d\n",&use_nfft);
      fprintf(stdout,"%d\n",use_nfsft);
      //fprintf(stdout,"  NFFT = %d\n",use_nfft);
       if (use_nfft != NO)
      {
        //fprintf(stdout,"  NFFT = yes\n");
        /* Read the cut-off parameter. */
        fscanf(stdin,"cutoff=%d\n",&cutoff);
        fprintf(stdout,"%d\n",cutoff);
        //fprintf(stdout,"  Cutoff = %d\n",cutoff);
      }
      else
      {
        //fprintf(stdout,"  NFFT = no\n");
        cutoff = 3;
      }
      /* Check if the FPT shall be used. */
      fscanf(stdin,"fpt=%d\n",&use_fpt);
      fprintf(stdout,"%d\n",use_fpt);
      //fprintf(stdout,"  FPT = %d\n",use_fpt);
      /* Read the threshold. */
      fscanf(stdin,"threshold=%lf\n",&threshold);
      fprintf(stdout,"%lf\n",threshold);
      //fprintf(stdout,"  Threshold = %E\n",threshold);
      }
    else
    {
      /* Set dummy values. */
      cutoff = 3;
      threshold = 1000000000000.0;
      //fprintf(stdout,"  NFSFT = no\n");
    }

    /* Initialize bandwidth bound. */
    m_max = 0;
    /* Initialize source node bound. */
    l_max = 0;
    /* Initialize target node bound. */
    d_max = 0;

    ld_max_prec = 0;
    l_max_prec = 0;

    /* Read kernel type. One of KT_ABEL_POISSON, KT_SINGULARITY, KT_LOC_SUPP
     * or KT_GAUSSIAN. */
    fscanf(stdin,"kernel=%d\n",&kt);
    fprintf(stdout,"%d\n",kt);

    /* Print out the kernel type. */
    //fprintf(stdout,"  Kernel type: %d\n",kt);

    /* Read the number of parameter sets. */
    fscanf(stdin,"parameter_sets=%d\n",&ip_max);
    fprintf(stdout,"%d\n",ip_max);

    /* Allocate memory for pointers to parameter sets. */
    p = (double**) malloc(ip_max*sizeof(double*));

    /* We now read in the parameter sets. */
    //fprintf(stdout,"  Parameter sets: %d\n",ip_max);

    fscanf(stdin,"parameters=%d\n",&ipp_max);
    fprintf(stdout,"%d\n",ipp_max);
    //fprintf(stdout,"  Parameters=%d\n",ipp_max);

    for (ip = 0; ip < ip_max; ip++)
    {
      p[ip] = (double*) malloc(ipp_max*sizeof(double));
      for (ipp = 0; ipp < ipp_max; ipp++)
      {
        /* Read parameter. */
        fscanf(stdin,"%lf\n",&p[ip][ipp]);
        fprintf(stdout,"%lf\n",p[ip][ipp]);
        /* Print out parameter. */
        //fprintf(stdout,"    %lf\n",p[ip][ipp]);
      }
    }

    /* Read number of cut-off degrees. */
    fscanf(stdin,"bandwidths=%d\n",&im_max);
    m = (int*) malloc(im_max*sizeof(int));

    /* Print out number of cut-off degrees. */
    fprintf(stdout,"%d\n",im_max);
    //fprintf(stdout,"  Bandwidths: %d\n",im_max);

    /* Read cut-off degrees. */
    for (im = 0; im < im_max; im++)
    {
      /* Read cut-off degree. */
      fscanf(stdin,"%d\n",&m[im]);
      fprintf(stdout,"%d\n",m[im]);
      m_max = MAX(m_max,m[im]);
      /* Print out cut-off degree. */
      //fprintf(stdout,"    M = %d\n",m[im]);
    }

    /* Read number of node specifications. */
    fscanf(stdin,"node_sets=%d\n",&ild_max);
    fprintf(stdout,"%d\n",ild_max);
    ld = (int**) malloc(ild_max*sizeof(int*));

    /* Print out number of node specifications. */
    //fprintf(stdout,"  Nodes: %d\n",ild_max);

      ld_max_prec = 0;
      l_max_prec = 0;
    for (ild = 0; ild < ild_max; ild++)
    {
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
      /* Print out parameters. */
      //fprintf(stdout,"    L = %d, D = %d, compare = %d",ld[ild][0],ld[ild][1],
      //  ld[ild][2]);
      if (ld[ild][2] == YES)
      {
        /* Read whether the precomputed version shall also be used. */
        fscanf(stdin,"precomputed=%d\n",&ld[ild][3]);
        fprintf(stdout,"%d\n",ld[ild][3]);
        /* Read the number of repetitions over which measurements are averaged. */
        fscanf(stdin,"repetitions=%d\n",&ld[ild][4]);
        fprintf(stdout,"%d\n",ld[ild][4]);
        //fprintf(stdout,", precomputed = %d, repetitions = %d",ld[ild][3],ld[ild][4]);
          if (ld[ild][3] == YES)
        {
          ld_max_prec = MAX(ld_max_prec,ld[ild][0]*ld[ild][1]);
          l_max_prec = MAX(l_max_prec,ld[ild][0]);
          precompute = YES;
        }
      }
      else
      {
        /* Set default value for the number of repetitions. */
        ld[ild][4] = 1;
      }
      //fprintf(stdout,"\n");
    }
    //fclose(file_dat);

    /* Print out the maximum cut-off degree. */
    //fprintf(stdout,"  Maximum M = %d\n",m_max);
    /* Print out the maximum number of source nodes. */
    //fprintf(stdout,"  Maximum L = %d\n",l_max);
    /* Print out the maximum number of target nodes. */
    //fprintf(stdout,"  Maximum D = %d\n",d_max);

    //n_max = 1<<ngpt(m_max);

    /* Allocate data structures. */
    b = (complex*) malloc(l_max*sizeof(complex));
    eta = (double*) malloc(2*l_max*sizeof(double));
    f_hat = (complex*) malloc(NFSFT_F_HAT_SIZE(m_max)*sizeof(complex));
    a = (complex*) malloc((m_max+1)*sizeof(complex));
    xi = (double*) malloc(2*d_max*sizeof(double));
    f_m = (complex*) malloc(d_max*sizeof(complex));
    f = (complex*) malloc(d_max*sizeof(complex));
    if (precompute == YES)
    {
      prec = (complex*) malloc(ld_max_prec*sizeof(complex));
    }

    /* Generate random source nodes and weights. */
    for (l = 0; l < l_max; l++)
    {
      b[l] = /*1.0;*/drand48() - 0.5;
      eta[2*l] = /*0.0;*/drand48() - 0.5;
      eta[2*l+1] = /*0.25;*/0.5*drand48();
      /*b[l] = 1.0;
      nu[2*l] = 0.0;
      nu[2*l+1] = 0.25;*/
    }

    /* Generate random target nodes. */
    for (d = 0; d < d_max; d++)
    {
      xi[2*d] = /*(d/(double)d_max)-0.5;*/drand48() - 0.5;
      xi[2*d+1] = /*0.25;*/0.5*drand48();
      /*xi[2*d] = (d/(double)d_max)-0.5;
      xi[2*d+1] = 0.25;*/
    }

    //sprintf(filename_tex,"testcase%d.tex",tc);
    //fprintf(stderr,"reached: m_max = %d, threshold = %lf\n",m_max,threshold);
    //fflush(stderr);

    nfsft_precompute(m_max,threshold,0U);
    //fprintf(stderr,"reached2!\n");
    //fflush(stderr);

    for (ip = 0; ip < ip_max; ip++)
    {
      //fprintf(stdout,"  Parameter set %d: ",ip);
      switch (kt)
      {
        case KT_ABEL_POISSON:
          //fprintf(stdout," h = %lf\n",p[ip][0]);
          break;
        case KT_SINGULARITY:
          //fprintf(stdout," h = %lf\n",p[ip][0]);
          break;
        case KT_LOC_SUPP:
          //fprintf(stdout," h = %lf, lambda = %lf\n",p[ip][0],p[ip][1]);
          break;
        case KT_GAUSSIAN:
          //fprintf(stdout," rho = %lf\n",p[ip][0]);
          break;
      }

      /* Kernel coeffcients up to m_max */
      switch (kt)
      {
         case KT_ABEL_POISSON:
          for (k = 0; k <= m_max; k++)
          {
            //fprintf(stderr,"a[%d] = %2.8lf\n",k,SYMBOL_ABEL_POISSON(k,p[ip][0]));
            a[k] = SYMBOL_ABEL_POISSON(k,p[ip][0]);
          }
          break;
        case KT_SINGULARITY:
          for (k = 0; k <= m_max; k++)
          {
            a[k] = SYMBOL_SINGULARITY(k,p[ip][0]);
          }
          break;
        case KT_LOC_SUPP:
          for (k = 0; k <= m_max; k++)
          {
            if (k == 0)
            {
              a[k] = 1.0;
            }
            else if (k == 1)
            {
              a[k] = ((p[ip][1]+1+p[ip][0])/(p[ip][1]+2.0))*a[k-1];
            }
            else
            {
              a[k] = (1.0/(k+p[ip][1]+1))*((2*k-1)*p[ip][0]*a[k-1] - (k-p[ip][1]-2)*a[k-2]);
            }
          }
          break;
        case KT_GAUSSIAN:
          /*steed = (double*) malloc(m_max*sizeof(double));
          gsl_sf_bessel_il_scaled_array(m_max,2*p[ip][0],steed);
          for (k = 0; k <= m_max; k++)
          {
            a[k] = 4*PI*steed[k];
          }
          free(steed);*/
          sprintf(filename_gaussian,"gaussian%.0f.dat",p[ip][0]);
          //fprintf(stderr,"filename = %s\n",filename_gaussian);
          file_gaussian = fopen(filename_gaussian,"r");
          if (file_gaussian != NULL)
          {
            fscanf(file_gaussian,"%d\n",&nsymbols);
            for (k = 0; k <= MIN(nsymbols,m_max); k++)
            {
              fscanf(file_gaussian,"%lf\n",&a[k]);
              //fprintf(stderr,"a[%d] = %.16E\n",k,a[k]);
            }
            for (k = nsymbols+1; k <= m_max; k++)
            {
              a[k] = 0.0;
            }
          }
          else
          {
            //fprintf(stderr,"Couldn't open file %s for reading!\n",filename_gaussian);
          }
          break;
      }

      for (k = 0; k <= m_max; k++)
      {
        a[k] *= (2*k+1)/(4*PI);
      }

      for (ild = 0; ild < ild_max; ild++)
      {
        //fprintf(stdout,"    L = %d, D = %d, ld_max_prec = %d, l_max_prec = %d\n",
        //ld[ild][0],ld[ild][1],ld_max_prec,l_max_prec);
        if (ld[ild][2] != NO)
        {
          /* Check if direct algorithm with precomputation should be tested. */
          if (ld[ild][3] != NO)
          {
            t_dp = 0.0;
            for (i = 0; i < ld[ild][4]; i++)
            {
              ptr = prec;
              rinc = l_max_prec-ld[ild][0];
              for (d = 0; d < ld[ild][1]; d++)
              {
                for (l = 0; l < ld[ild][0]; l++)
                {
                  temp = innerProduct(2*PI*eta[2*l],2*PI*eta[2*l+1],2*PI*xi[2*d],2*PI*xi[2*d+1]);
                  switch (kt)
                  {
                    case KT_ABEL_POISSON:
                      *ptr++ = poissonKernel(temp,p[ip][0]);
                     break;
                   case KT_SINGULARITY:
                     *ptr++ = singularityKernel(temp,p[ip][0]);
                     break;
                   case KT_LOC_SUPP:
                     *ptr++ = locSuppKernel(temp,p[ip][0],p[ip][1]);
                     break;
                   case KT_GAUSSIAN:
                     *ptr++ = gaussianKernel(temp,p[ip][0]);
                      break;
                 }
               }
               ptr += rinc;
             }
             ptr = prec;
             rinc = l_max_prec-ld[ild][0];
             if (kt == KT_LOC_SUPP)
             {
               constant = ((p[ip][1]+1)/(2*PI*pow(1-p[ip][0],p[ip][1]+1)));
               t = second();
               for (d = 0; d < ld[ild][1]; d++)
               {
                 f[d] = 0.0;
                 for (l = 0; l < ld[ild][0]; l++)
                 {
                   f[d] += b[l]*(*ptr++);
                 }
                 f[d] *= constant;
                 ptr += rinc;
               }
               t_dp += second() - t;
             }
             else
             {
               t = second();
               for (d = 0; d < ld[ild][1]; d++)
               {
                 f[d] = 0.0;
                 for (l = 0; l < ld[ild][0]; l++)
                 {
                   f[d] += b[l]*(*ptr++);
                 }
                 ptr += rinc;
               }
               t_dp += second() - t;
             }
           }
           t_dp = t_dp/((double)ld[ild][4]);
           //printf("      t_dp = %f\n",t_dp);
           }
         else
         {
           t_dp = -1.0;
         }

         t_d = 0.0;
         for (i = 0; i < ld[ild][4]; i++)
         {
           switch (kt)
           {
              case KT_ABEL_POISSON:
               t = second();
               for (d = 0; d < ld[ild][1]; d++)
               {
                  f[d] = 0.0;
                  for (l = 0; l < ld[ild][0]; l++)
                 {
                   temp = innerProduct(2*PI*eta[2*l],2*PI*eta[2*l+1],2*PI*xi[2*d],2*PI*xi[2*d+1]);
                   f[d] += b[l]*poissonKernel(temp,p[ip][0]);
                 }
               }
               t_d += second() - t;
               break;
             case KT_SINGULARITY:
               t = second();
               for (d = 0; d < ld[ild][1]; d++)
               {
                 f[d] = 0.0;
                 for (l = 0; l < ld[ild][0]; l++)
                 {
                   temp = innerProduct(2*PI*eta[2*l],2*PI*eta[2*l+1],2*PI*xi[2*d],2*PI*xi[2*d+1]);
                   f[d] += b[l]*singularityKernel(temp,p[ip][0]);
                 }
               }
               t_d += second() - t;
               break;
             case KT_LOC_SUPP:
               constant = ((p[ip][1]+1)/(2*PI*pow(1-p[ip][0],p[ip][1]+1)));
               t = second();
               for (d = 0; d < ld[ild][1]; d++)
               {
                 f[d] = 0.0;
                 for (l = 0; l < ld[ild][0]; l++)
                 {
                   temp = innerProduct(2*PI*eta[2*l],2*PI*eta[2*l+1],2*PI*xi[2*d],2*PI*xi[2*d+1]);
                   f[d] += b[l]*locSuppKernel(temp,p[ip][0],p[ip][1]);
                 }
                 f[d] *= constant;
               }
               t_d += second() - t;
               break;
             case KT_GAUSSIAN:
               t = second();
               for (d = 0; d < ld[ild][1]; d++)
               {
                 f[d] = 0.0;
                 for (l = 0; l < ld[ild][0]; l++)
                 {
                   temp = innerProduct(2*PI*eta[2*l],2*PI*eta[2*l+1],2*PI*xi[2*d],2*PI*xi[2*d+1]);
                   f[d] += b[l]*gaussianKernel(temp,p[ip][0]);
                 }
               }
               t_d = second() - t;
               break;
            }
         }
         t_d = t_d/((double)ld[ild][4]);
           //printf("      t_d = %f\n",t_d);
       }
      else
      {
        t_d = -1.0;
        t_dp = -1.0;
      }

      err_fd = -1.0;
      err_f = -1.0;
      t_fd = -1.0;
      t_f = -1.0;

      //file_dat = fopen(filename_dat,"a");
      //fprintf(file_dat,"\n");
      //fclose(file_dat);

      for (im = 0; im < im_max; im++)
      {
        //fprintf(stderr,"      M = %d:",m[im]);

        /* Init transform plans. */
        nfsft_init_guru(&plan_adjoint,m[im],ld[ild][0],
          ((use_nfft!=0)?(0U):(NFSFT_USE_NDFT)) |
          ((use_fpt!=0)?(0U):(NFSFT_USE_DPT)), cutoff);
        nfsft_init_guru(&plan,m[im],ld[ild][1],
          ((use_nfft!=0)?(0U):(NFSFT_USE_NDFT)) |
          ((use_fpt!=0)?(0U):(NFSFT_USE_DPT)), cutoff);
        plan_adjoint.f_hat = f_hat;
        plan_adjoint.x = eta;
        plan_adjoint.f = b;
        plan.f_hat = f_hat;
        plan.x = xi;
        plan.f = f_m;

        if (use_nfsft == BOTH)
        {
          t_fd = 0.0;
          for (i = 0; i < ld[ild][4]; i++)
          {
            t = second();
            ndsft_adjoint(&plan_adjoint);
           /* Multiplication with diagonal matrix. */
           for (k = 0; k <= m[im]; k++)
           {
             for (n = -k; n <= k; n++)
             {
               f_hat[NFSFT_INDEX(k,n,&plan_adjoint)] *= a[k];
             }
           }
          ndsft_trafo(&plan);
          t_fd += second() - t;
          if (ld[ild][2] != NO)
          {
            err_fd = error_l_infty_1_complex(f, f_m, ld[ild][1], b, ld[ild][0]);
            //printf("\terr_fd = %le\n",err_fd);
          }
        }
          t_fd = t_fd/((double)ld[ild][4]);
        //printf("\tt_fd = %f",t_fd);
      }

      if (use_nfsft != NO)
      {
        t_f = 0.0;
      }
      else
      {
        t_fd = 0.0;
      }
      for (i = 0; i < ld[ild][4]; i++)
        {
         /* Adjoint transform */
         t = second();
         if (use_nfsft != NO)
         {
           nfsft_adjoint(&plan_adjoint);
         }
         else
        {
          ndsft_adjoint(&plan_adjoint);
         }

         /* Multiplication with diagonal matrix. */
         for (k = 0; k <= m[im]; k++)
        {
          for (n = -k; n <= k; n++)
          {
            f_hat[NFSFT_INDEX(k,n,&plan_adjoint)] *= a[k];
           }
         }

        /* Forward transform */
        if (use_nfsft != NO)
        {
          nfsft_trafo(&plan);
        }
        else
        {
           ndsft_trafo(&plan);
         }
         if (use_nfsft != NO)
         {
           t_f += second() - t;
         }
         else
         {
           t_fd += second() - t;
         }

         if (ld[ild][2] != NO)
         {
           if (use_nfsft != NO)
           {
             err_f = error_l_infty_1_complex(f, f_m, ld[ild][1], b, ld[ild][0]);
             //printf("\terr_f = %le\n",err_f);
           }
           else
           {
             err_fd = error_l_infty_1_complex(f, f_m, ld[ild][1], b, ld[ild][0]);
             ///printf("\terr_fd = %le\n",err_fd);
           }
         }
       }
       if (use_nfsft != NO)
       {
         t_f = t_f/((double)ld[ild][4]);
         //printf("\tt_f = %f",t_f);
       }
       else
       {
         t_fd = t_fd/((double)ld[ild][4]);
         //printf("\tt_f = %f",t_fd);
       }

       /*for (d = 0; d < ld[ild][1]; d++)
       {
         fprintf(stderr,"%+5.16f, %+5.16f, %+.3E\n",creal(f_m[d]),creal(f[d]),creal(f_m[d]-f[d]));
         fflush(stderr);
       }*/
      /*fprintf(stderr,"\n");
      for (l = 0; l < ld[ild][0]; l++)
      {
        fprintf(stderr,"%+5.16f\n",creal(b[l]));
        fflush(stderr);
      }*/

      /*for (k = 0; k < ld[ild][1]; k++)
      {
        fprintf(stderr,"f[%d] = %le + I*%le \tf_m[%d] = %le + I*%le\n",
          k,creal(f[k]),cimag(f[k]),k,creal(f_m[k]),cimag(f_m[k]));
      }*/

      //fprintf(stderr,"||b||_1 = %E, L = %d\n",norm_complex_1(b,ld[ild][0]),ld[ild][0]);

      //file_tex = fopen(filename_tex,"a");
      //file_dat = fopen(filename_dat,"a");

      /*fprintf(file_tex,"%6d & %6d & %.1E & %.1E & %.1E & %.1E & %.1E\\\\\n",
      ld[ild][0],ld[ild][1],t_d,t_dp,t_fd,t_f,error_l_infty_1_complex(f, f_m,
      b, ld[ild][1], ld[ild][0]));*/
      fprintf(stdout,"%e\n%e\n%e\n%e\n%e\n%e\n",t_d,t_dp,t_fd,t_f,
        err_fd,err_f);
      //fclose(file_dat);

      /* Finalize plans */
      nfsft_finalize(&plan_adjoint);
      nfsft_finalize(&plan);
    }
  }
}

nfsft_forget();

    if (precompute == YES)
    {
      free(prec);
    }
    free(f);
    free(f_m);
    free(xi);
    free(eta);
    free(a);
    free(f_hat);
    free(b);

    for (ild = 0; ild < ild_max; ild++)
    {
      free(ld[ild]);
    }
    free(ld);

    free(m);

    for (ip = 0; ip < ip_max; ip++)
    {
      free(p[ip]);
    }
    free(p);

    //fprintf(stdout,"\n");
  }

  return EXIT_SUCCESS;
}
