/* 
   accuracy - Fast spherical convolution example for the NFSFT

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

#ifdef STDC_HEADERS
#  include <stdlib.h>
#  include <stdio.h>
#  include <math.h>
#else
#  error Need ANSI-C headers
#endif

/* Auxilliary headers */
#include <complex.h>
#include "nfsft.h"
#include "util.h"
#include "../../nfft/utils.h"
#include "time.h"

#define M_MIN 32
#define M_MAX 256
#define M_STRIDE 32

#define H_MIN 0.7
#define H_MAX 0.9
#define H_STRIDE 0.2

#define D_MIN 1000
#define D_MAX 10000
#define D_STRIDE 1000

#define DD_MAX 40000

#define L_MIN 1000
#define L_MAX 10000
#define L_STRIDE 1000

#define SYMBOL_ABEL_POISSON(k,h) pow(h,k)
#define SYMBOL_SINGULARITY(k,h) (2.0/(2*k+1))*pow(h,k)

inline double ip(double phi1,double theta1,double phi2,double theta2)
{
  return cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2);
}

inline double poissonKernel(double phi1,double theta1,double phi2,double theta2,
                            double h)
{
  double t = ip(phi1,theta1,phi2,theta2);
 	return (1.0/(4*PI))*(1-h*h)/pow(sqrt(1-2*h*t+h*h),3);
}

inline double singularityKernel(double phi1,double theta1,double phi2,double theta2,
                            double h)
{
  double t = ip(phi1,theta1,phi2,theta2);
 	return (1.0/(2*PI))/sqrt(1-2*h*t+h*h);
}

inline double locSupKernel(double phi1,double theta1,double phi2,double theta2,
                           double h, double lambda)
{
  double t = ip(phi1,theta1,phi2,theta2);
 	return (t<=h)?(0.0):(pow((t-h),lambda));
}

inline int maxv(int *x, int n)
{
  int i;
  int r = 0;
  
  for (i = 0; i < n; i++)
  {
    if (x[i] > r)
    {
      r = x[i];
    }
  }
  return r;
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
  double *h;
  double *lambda;
  int *m;
  int *l;
  int *d;
  int ih;
  int ilambda;
  int im;
  int il;
  int id;
  int ds;
  int kt;
  int ih_max;
  int ilambda_max;
  int im_max;
  int il_max;
  int id_max;
  int ds_max;
  int m_max;
  int l_max;
  int d_max;
  
  /** Next greater power of two with respect to M_MAX */
  int N_MAX;
  /* Bandwidth */
  //int M;
  /** Next greater power of two with respect to M */
  //int N;
 	/** Error */
	 double err;
	 double t_d, t_f;
  double nfactor;
	
 	complex *b;
  complex **f_hat;
 	complex *a;
 	double *xi;
 	double *nu;
 	complex *f, *f2;
 	nfsft_plan plan;
  nfsft_plan plan_adjoint;
 	int j,k,n;
	 FILE *file;
  FILE *file2;
  char filename[100];
  char filename2[100];
	 	
  fscanf(stdin,"%d",&ds_max);
  for (ds = 0; ds < ds_max; ds++)
  {
    fscanf(stdin,"%d",&kt);
    fscanf(stdin,"%d",&ih_max);
    h = (double*) malloc(ih_max*sizeof(double));

    if (kt == 2)
    {
      ilambda_max = ih_max;
      lambda = (double*) malloc(ilambda_max*sizeof(double));
    }
    
    for (j = 0; j < ih_max; j++)
    {  
      fscanf(stdin,"%lf",&h[j]);
      if (kt == 2)
      {
        fscanf(stdin,"%lf",&lambda[j]);
      }
    }

    fscanf(stdin,"%d",&im_max);
    m = (int*) malloc(im_max*sizeof(int));

    for (j = 0; j < im_max; j++)
    {  
      fscanf(stdin,"%d",&m[j]);
    }
    
    fscanf(stdin,"%d",&il_max);
    id_max = il_max;
    l = (int*) malloc(il_max*sizeof(int));
    d = (int*) malloc(id_max*sizeof(int));
    
    for (j = 0; j < il_max; j++)
    {  
      fscanf(stdin,"%d",&l[j]);
      fscanf(stdin,"%d",&d[j]);
    }
    
    m_max = maxv(m,im_max);
    l_max = maxv(l,il_max);
    d_max = maxv(d,id_max);
    N_MAX = 1<<ngpt(m_max);	

    /** Allocate data structures. */
    b = (complex*) malloc(l_max*sizeof(complex));
    nu = (double*) malloc(2*l_max*sizeof(double));
    f_hat = (complex**) malloc((2*m_max+1)*sizeof(complex*));
    for (n = -m_max; n <= m_max; n++)
    {
      f_hat[n+m_max] = (complex*) malloc((N_MAX+1)*sizeof(complex));
    }  
    a = (complex*) malloc((m_max+1)*sizeof(complex));
    xi = (double*) malloc(2*d_max*sizeof(double));
    f = (complex*) malloc(d_max*sizeof(complex));
    f2 = (complex*) malloc(d_max*sizeof(complex));
    
    /* Target nodes */
    for (j = 0; j < d_max; j++)
    {
      xi[2*j] = drand48()-0.5;
      xi[2*j+1] = 0.5*drand48();
      /*xi[2*j] = ((double)j/(d_max)) -0.5;
      xi[2*j+1] = 0.25;
      fprintf(stderr,"Target node: (%f,%f)\n",xi[2*j],xi[2*j+1]);*/
    }

    /* Source nodes. */
    for (k = 0; k < l_max; k++)
    {
      b[k] = drand48();
      nu[2*k] = drand48()-0.5;
      nu[2*k+1] = 0.5*drand48();
      /*b[k] = 1.0;
      nu[2*k] = 0.0;
      nu[2*k+1] = 0.25;
      fprintf(stderr,"Source node: (%f,%f)\n",nu[2*k],nu[2*k+1]);*/
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
      }
      /* Kernel coeffcients up to m_max */
      for (k = 0; k <= m_max; k++)
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
              a[k] = /*2*PI*(pow(1.0-h[ih],lambda[ilambda]+1.0)/(lambda[ilambda]+1.0));*/1.0;
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

      fprintf(stderr,"h[%d] = %f\n",ih,h[ih]);
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
            fprintf(stderr,"f(%f,%f) = %f\n",xi[2*j],xi[2*j+1],f2[j]);
            fflush(stderr);
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
          fflush(stderr);
          
          /* Init transform plans. */
          plan_adjoint = nfsft_init(l[il],m[im],nu,f_hat,b,0U);
          plan = nfsft_init(d[id],m[im],xi,f_hat,f,0U);
          
          /* Adjoint transform */
          t_f = second();
          nfsft_adjoint(plan_adjoint);
          //fprintf(stderr,"1");
          //fflush(stderr);
            
          /* Multiplication with diagonal matrix. */
          for (k = 0; k <= m[im]; k++)
          {
            for (n = -k; n <= k; n++)
            {
              f_hat[n+m[im]][k] *= a[k];
            }
          }
          //fprintf(stderr,"2");
          //fflush(stderr);
            
          /* Forward transform */
          nfsft_trafo(plan);
          t_f = second() - t_f;

          /* Finalize plans */
          nfsft_finalize(plan_adjoint);
          nfsft_finalize(plan);
          //fprintf(stderr,"3");
          //fflush(stderr);
            
          for (j = 0; j < d[id]; j++)
          {
            fprintf(stderr,"f_a(%f,%f) = %f\n",xi[2*j],xi[2*j+1],f[j]);
            fflush(stderr);
          }
          
          if (d[id] <= DD_MAX)
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
    
    nfsft_forget();
    
    free(f);
    free(f2);
    free(xi);
    free(nu);
    free(a);
    for (n = -m_max; n <= m_max; n++)
    {
      free(f_hat[n+m_max]);
    }   
    free(f_hat);
    free(b);
    
    free(h);
    free(m);
    free(l);
    free(d);
    if (kt == 2)
    {
      free(lambda);
    }
  }  
  
  return EXIT_SUCCESS;
}
