/* 
   Fastsum - Fast summation of spherical radial basis functions

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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#else
#  error Need config.h
#endif

/* Auxilliary headers */
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

/* Library headers. */
#include "nfsft.h"
#include "util.h"
#include "../../nfft/utils.h"

/** The symbol of the Abel-Poisson kernel */
#define SYMBOL_ABEL_POISSON(k,h) pow(h,k)
/** The symbol of the singularity kernel */
#define SYMBOL_SINGULARITY(k,h) (2.0/(2*k+1))*pow(h,k)

/* Kernel types */
/** Abel-Poisson kernel */
#define KT_ABEL_POISSON (0)
/** Singularity kernel */
#define KT_SINGULARITY  (1)
/** Locally supported kernel */
#define KT_LOC_SUPP     (2)
/** Gaussian kernel */
#define KT_GAUSSIAN     (3)

/** Computes the maximum of two integers. */
inline int max(const int a, const int b)
{
  return (a>b)?a:b;
}

/** Computes the inner product on \f$\mathbb{S}^2\f$. */
inline double ip(const double phi1, const double theta1, const double phi2,
                 const double theta2)
{
  return cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2);
}

/** Evaluates the Poisson kernel. */
inline double poissonKernel(const double x, const double h)
{
 	return (1.0/(4*PI))*(1-h*h)/pow(sqrt(1-2*h*x+h*h),3);
}

/** Evaluates the singularity kernel. */
inline double singularityKernel(const double x, const double h)
{
  return (1.0/(2*PI))/sqrt(1-2*h*x+h*h);
}

/** Evaluates the locally supported kernel. */
inline double locSuppKernel(const double x, const double h, const double lambda)
{
 	return (x<=h)?(0.0):(pow((x-h),lambda));
}

/** Evaluates the Gaussian kernel. */
inline double gaussianKernel(const double x, const double rho)
{
 	return exp(2*rho*(x-1));
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
  /** An array containing the parameter sets for the current kernel. */
  double **p;
  /** An array containing the parameters \f$M\f$. */
  int *m;
  /** An array containing the parameters \f$L\f$ and \f$D\f$. */
  int **ld;
  /** Index variable for p */
  int ip;
  /** Index variable for m */
  int im;
  /** Index variable for l */
  int ild;
  /** Maximum index for p */
  int ip_max;
  /** Maximum index for m */
  int im_max;
  /** Maximum index for l */
  int ild_max;

  /** Number of testcases */
  int tc_max;
  /** Maximum \f$M\f$ for the current dataset */
  int m_max;
  /** Maximum \f$L\f$ for the current dataset */
  int l_max;
  /** Maximum \f$D\f$ for the current dataset */
  int d_max;
  
  /** Index variable for testcases. */
  int tc;
  /** The kernel type */
  int kt;
  
  /** Next greater power of two with respect to m_max */
  int n_max;
 	/** 
   * The relative error
   * 
   * \f[
       \frac{\left\|f-f_M\right\|_{\infty}}{\left\|\mathbf{b}\right\|}
   * \f]
   */
	 double err;

	 double t_d, t_f;
  double nfactor;
	
  /** Weights \f$\left(b_l\right)_{l=0}^{L-1}\f$ */
 	complex *b;
  /** Fourier coefficients */
  complex **f_hat;
  /** Symbool coefficients */
 	complex *a;
  /** Target nodes */
 	double *xi;
  /** Source nodes */
 	double *nu;
  /** Approximate function values */
 	complex *f_m;
  /** Exact function values */
  complex *f;
  /** NFSFT plan */
 	nfsft_plan plan;
  /** adjoint NFSFT plan */
  nfsft_plan plan_adjoint;
  
 	int j,k,n;
  
	 FILE *file_tex;
  FILE *file_dat;
  char filename_tex[100];
  char filename_dat[100];
	 	
  /* Read number of testcases. */
  fscanf(stdin,"%d",&tc_max);
  
  fprintf(stdout,"Number of testcases: %d\n\n",tc_max);
  
  /* Process testcases. */
  for (tc = 0; tc < tc_max; tc++)
  {
    fprintf(stdout,"Testcase %d:\n",tc);
    /* Initialize bandwidth bound. */
    m_max = 0;
    /* Initialize source node bound. */
    l_max = 0;
    /* Initialize target node bound. */
    d_max = 0;
    
    /* Read kernel type. One of KT_ABEL_POISSON, KT_SINGULARITY, KT_LOC_SUPP 
     * or KT_GAUSSIAN. */
    fscanf(stdin,"%d",&kt);

    fprintf(stdout,"  Kernel type: %d\n",kt);    
    
    fscanf(stdin,"%d",&ip_max);
    p = (double**) malloc(ip_max*sizeof(double*));

    fprintf(stdout,"  Parameter sets: %d\n",ip_max);    
    
    switch (kt)
    {
      case KT_ABEL_POISSON:
      case KT_SINGULARITY:
      case KT_GAUSSIAN:
        for (ip = 0; ip < ip_max; ip++)
        {  
          p[ip] = (double*) malloc(1*sizeof(double));
          fscanf(stdin,"%lf",&p[ip][0]);
          if (kt == KT_GAUSSIAN)
          {
            fprintf(stdout,"    rho = %lf\n",p[ip][0]);    
          }
          else
          {
            fprintf(stdout,"    h = %lf\n",p[ip][0]);    
          }
        }
        break;
      case KT_LOC_SUPP:
        for (ip = 0; ip < ip_max; ip++)
        {  
          p[ip] = (double*) malloc(2*sizeof(double));
          fscanf(stdin,"%lf",&p[ip][0]);
          fscanf(stdin,"%lf",&p[ip][1]);
          fprintf(stdout,"    h = %lf, lambda = %lf\n",p[ip][0],p[ip][1]);    
        }
        break;
    };
    
    /* Read number of bandwidths. */
    fscanf(stdin,"%d",&im_max);
    m = (int*) malloc(im_max*sizeof(int));

    fprintf(stdout,"  Bandwidths: %d\n",im_max);    
    
    /* Read bandwidths. */
    for (im = 0; im < im_max; im++)
    {  
      fscanf(stdin,"%d",&m[im]);
      m_max = max(m_max,m[im]);
      fprintf(stdout,"    M = %d\n",m[im]);    
    }
    
    /* Read number of nodes specifications. */
    fscanf(stdin,"%d",&ild_max);
    ld = (int**) malloc(ild_max*sizeof(int*));
    
    fprintf(stdout,"  Nodes: %d\n",ild_max);    

    for (ild = 0; ild < ild_max; ild++)
    {  
      ld[ild] = (int*) malloc(2*sizeof(int));
      fscanf(stdin,"%d",&ld[ild][0]);
      l_max = max(l_max,ld[ild][0]);
      fscanf(stdin,"%d",&ld[ild][1]);
      d_max = max(d_max,ld[ild][1]);
      fprintf(stdout,"    L = %d, D = %d\n",ld[ild][0],ld[ild][1]);    
    }
    
    fprintf(stdout,"  Maximum M = %d\n",m_max);    
    fprintf(stdout,"  Maximum L = %d\n",l_max);    
    fprintf(stdout,"  Maximum D = %d\n",d_max);    
    
    continue;
    
    n_max = 1<<ngpt(m_max);	

    /** Allocate data structures. */
    b = (complex*) malloc(l_max*sizeof(complex));
    nu = (double*) malloc(2*l_max*sizeof(double));
    f_hat = (complex**) malloc((2*m_max+1)*sizeof(complex*));
    for (n = -m_max; n <= m_max; n++)
    {
      f_hat[n+m_max] = (complex*) malloc((n_max+1)*sizeof(complex));
    }  
    a = (complex*) malloc((m_max+1)*sizeof(complex));
    xi = (double*) malloc(2*d_max*sizeof(double));
    f_m = (complex*) malloc(d_max*sizeof(complex));
    f = (complex*) malloc(d_max*sizeof(complex));
    
    /* Target nodes */
    /*for (j = 0; j < d_max; j++)
    {
      xi[2*j] = drand48()-0.5;
      xi[2*j+1] = 0.5*drand48();
    }*/

    /* Source nodes. */
    /*for (k = 0; k < l_max; k++)
    {
      b[k] = drand48();
      nu[2*k] = drand48()-0.5;
      nu[2*k+1] = 0.5*drand48();
      fprintf(stderr,"Source node: (%f,%f)\n",nu[2*k],nu[2*k+1]);
    }

    sprintf(filename,"summation%d.tex",kt);
    sprintf(filename2,"summation%d.dat",kt);
    file = fopen(filename,"w");
    file2 = fopen(filename2,"w");
    if (kt == 0 || kt == 1)
    {  
      fprintf(file,"\\begin{tabular}{l|l|l|l|l|l|l|l}\n");
      fprintf(file,"$h$ & $M$ & $K$ & $N$ & $t_{\\text{fast}}$ & $t_{\\text{slow}}$ & $\\text{err}_{\\infty}$ & $\\text{err}_{2}$\\\\\\hline\n");
    }
    else if (kt == 2)
    {
      fprintf(file,"\\begin{tabular}{l|l|l|l|l|l|l|l|l}\n");
      fprintf(file,"$h$ & $\\lambda$ & $M$ & $K$ & $N$ & $t_{\\text{fast}}$ & $t_{\\text{slow}}$ & $\\text{err}_{\\infty}$ & $\\text{err}_{2}$\\\\\\hline\n");
    }

    nfsft_precompute(m_max,2000,0U);
    
    for (ih = 0; ih < ih_max; ih++)
    {
      if (kt == 2)
      {
        ilambda = ih;
      }*/
      /* Kernel coeffcients up to m_max */
      /*for (k = 0; k <= m_max; k++)
      {
        switch (kt)
        {
          case 0:
            a[k] = SYMBOL_ABEL_POISSON(k,h[ih]);            
            break;
          case 1:
            a[k] = SYMBOL_SINGULARITY(k,h[ih]);
            break;
          case 2:
            if (k == 0)
            {
              a[k] = 1.0;
            }
            else if (k == 1)
            {
              a[k] = ((lambda[ilambda]+1+h[ih])/(lambda[ilambda]+2.0))*a[k-1];
            }
            else
            {
              a[k] = (1.0/(k+lambda[ilambda]+1))*((2*k-1)*h[ih]*a[k-1] - (k-lambda[ilambda]-2)*a[k-2]);
            }
            break;                
        }
      }

      //fprintf(stderr,"h[%d] = %f\n",ih,h[ih]);
      fprintf(stderr,"lambda[%d] = %f\n",ilambda,lambda[ilambda]);
      fflush(stderr);
      
      for (k = 0; k <= m_max; k++)
      {
        a[k] *= (2*k+1)/(4*PI);
      }
      
      for (il = 0; il < il_max; il++)
      {       
        id = il;
        
        if (d[id] <= DD_MAX)
        {
          t_d = second();
          for (j = 0; j < d[id]; j++)
          {
            f2[j] = 0.0;
            for (k = 0; k < l[il]; k++)
            {
              switch (kt)
              {
                case 0:
                  f2[j] += b[k]*poissonKernel(2*PI*nu[2*k],2*PI*nu[2*k+1],2*PI*xi[2*j],2*PI*xi[2*j+1],h[ih]);
                  break;
                case 1:
                  f2[j] += b[k]*singularityKernel(2*PI*nu[2*k],2*PI*nu[2*k+1],2*PI*xi[2*j],2*PI*xi[2*j+1],h[ih]);
                  break;
                case 2:
                  f2[j] += b[k]*locSupKernel(2*PI*nu[2*k],2*PI*nu[2*k+1],2*PI*xi[2*j],2*PI*xi[2*j+1],h[ih],lambda[il]);
                  break;                
              }
            }
            if (kt == 2)
            {
              f2[j] *= ((lambda[ilambda]+1)/(2*PI*pow(1-h[ih],lambda[ilambda]+1)));
            }
            //fprintf(stderr,"f(%f,%f) = %f\n",xi[2*j],xi[2*j+1],f2[j]);
            //fflush(stderr);
          } 
          t_d = second() - t_d;
        }
          
        for (im = 0; im < im_max; im++)
        {
          if (kt == 0 || kt == 1)
          {
            fprintf(stderr,"L = %d, D = %d, h = %lf, M = %d\n",l[il],d[id],h[ih],m[im]);
          }
          else if (kt == 2)
          {
            fprintf(stderr,"L = %d, D = %d, h = %lf, lambda = %lf, M = %d\n",l[il],d[id],h[ih],lambda[ilambda],m[im]);
          }
          fflush(stderr);*/
          
          /* Init transform plans. */
/*          plan_adjoint = nfsft_init(m[im],l[il],f_hat,nu,b,0U);
          plan = nfsft_init(m[im],d[id],f_hat,xi,f,0U);*/
          
          /* Adjoint transform */
//          t_f = second();
//          nfsft_adjoint(plan_adjoint);
            
          /* Multiplication with diagonal matrix. */
/*          for (k = 0; k <= m[im]; k++)
          {
            for (n = -k; n <= k; n++)
            {
              f_hat[n+m[im]][k] *= a[k];
            }
          }*/
            
          /* Forward transform */
//          nfsft_trafo(plan);
//          t_f = second() - t_f;

          /* Finalize plans */
//          nfsft_finalize(plan_adjoint);
//          nfsft_finalize(plan);
            
          //for (j = 0; j < d[id]; j++)
          //{
          //  fprintf(stderr,"f_a(%f,%f) = %f\n",xi[2*j],xi[2*j+1],f[j]);
          //  fflush(stderr);
          //}
          
          /*if (d[id] <= DD_MAX)
          {
            if (kt == 0 || kt == 1)
            {
              fprintf(file,"%.2f & %3d & %6d & %6d & %5.2f & %5.2f & %.4E & %4.E\\\\\n",h[ih],m[im],l[il],d[id],t_d,t_f,error_complex_inf(f, f2, d[id]),error_complex_2(f, f2, d[id]));
              fprintf(file2,"%.2f %3d %6d %6d %5.2f %5.2f %.4E %4.E\\\\\n",h[ih],m[im],l[il],d[id],t_d,t_f,error_complex_inf(f, f2, d[id]),error_complex_2(f, f2, d[id]));
            }
            else if (kt == 2)
            {
              fprintf(file,"%.2f & %2f & %3d & %6d & %6d & %5.2f & %5.2f & %.4E & %4.E\\\\\n",h[ih],lambda[ilambda],m[im],l[il],d[id],t_d,t_f,error_complex_inf(f, f2, d[id]),error_complex_2(f, f2, d[id]));
              fprintf(file2,"%.2f %2f %3d %6d %6d %5.2f %5.2f %.4E %4.E\\\\\n",h[ih],lambda[ilambda],m[im],l[il],d[id],t_d,t_f,error_complex_inf(f, f2, d[id]),error_complex_2(f, f2, d[id]));
            }
          }
          else
          {
            if (kt == 0 || kt == 1)
            {
              fprintf(file,"%.2f & %3d & %6d & %6d & -- & %5.2f & -- & --\\\\\n",h[ih],m[im],l[il],d[id],t_f);
              fprintf(file2,"%.2f %3d %6d %6d -1.0 %5.2f -1E0 -1E0\\\\\n",h[ih],m[im],l[il],d[id],t_f);
            }
            else if (kt == 2)
            {
              fprintf(file,"%.2f & %2f & %3d & %6d & %6d & -- & %5.2f & -- & --\\\\\n",h[ih],lambda[ilambda],m[im],l[il],d[id],t_f);
              fprintf(file2,"%.2f %2f %3d %6d %6d -1.0 %5.2f -1E0 -1E0\\\\\n",h[ih],lambda[ilambda],m[im],l[il],d[id],t_f);
            }
          }
        }
      }
      fprintf(file,"\\hline\n");
    }
    
    fprintf(file,"\\end{tabular}\n");
    fclose(file);
    fclose(file2);
    
    nfsft_forget();*/
    
    free(f);
    free(f_m);
    free(xi);
    free(nu);
    free(a);
    for (n = -m_max; n <= m_max; n++)
    {
      free(f_hat[n+m_max]);
    }   
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
    
    fprintf(stdout,"\n");
  }  
  
  return EXIT_SUCCESS;
}
