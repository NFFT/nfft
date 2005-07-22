/* 
   accuracy2 - Accuracy test for the NFSFT

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
#  include <string.h>
#else
#  error Need ANSI-C headers
#endif

/* Auxilliary headers */
#include <complex.h>
#include <time.h>
#include <fftw3.h>
#include "nfsft.h"
#include "util.h"
#include "../../nfft/utils.h"

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

#define NO 0
#define YES 1

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
  /** An array containing the parameters \f$M\f$. */
  int *m;
  /** The parameter \f$D\f$. */
  int D;
  /** Index variable for m */
  int im;
  /** Maximum index for m */
  int im_max;
  
  /** Number of testcases */
  int tc_max;
  /** Maximum \f$M\f$ for the current dataset */
  int m_max;
  
  /** Index variable for testcases. */
  int tc;
  int cutoff;
  double threshold;
  
  /** Next greater power of two with respect to m_max */
  int n_max;

  double err2, err_inf;

  double nfactor;
  double temp;
 	double t;
 	complex *ptr;
	
  /** Fourier coefficients */
  complex **f_hat;
  /** Fourier coefficients */
  complex **f_hat2;
  complex **f_hatr;
  /** Target nodes */
 	double *xi;
  /** Approximate function values */
 	complex *f;
 	complex *f2;
 	complex *fr;
  /** NFSFT plan */
 	nfsft_plan plan;
  /** adjoint NFSFT plan */
  //nfsft_plan plan_adjoint;
  
 	int i,j,k,n,d,l,use_nfft;
 	long index;
  
  FILE *file_tex;
  FILE *file_dat;
  char filename_tex[100];
  char filename_dat[100];
	 	
  /* Read number of testcases. */
  fscanf(stdin,"testcases=%d\n",&tc_max);
  
  fprintf(stdout,"Number of testcases: %d\n\n",tc_max);
  
  /* Process testcases. */
  for (tc = 0; tc < tc_max; tc++)
  {
    fprintf(stdout,"Testcase %d:\n",tc);

    fscanf(stdin,"nfft=%d\n",&use_nfft);
    fprintf(stdout,"  NFFT = %d\n",use_nfft);
		 	if (use_nfft != NO)
    {
      fprintf(stdout,"  NFFT = yes\n");
      fscanf(stdin,"cutoff=%d\n",&cutoff);
      fprintf(stdout,"  Cutoff = %d\n",cutoff);
	   }
    else
    {
      fprintf(stdout,"  NFFT = no\n");
      cutoff = 3;
    }     
    fscanf(stdin,"threshold=%lf\n",&threshold);
    fprintf(stdout,"  Threshold = %E\n",threshold);

    /* Initialize bandwidth bound. */
    m_max = 0;
        
    /* Read number of bandwidths. */
    fscanf(stdin,"bandwidths=%d\n",&im_max);
    m = (int*) malloc(im_max*sizeof(int));

    fprintf(stdout,"  Bandwidths: %d\n",im_max);    
    
    /* Read bandwidths. */
    for (im = 0; im < im_max; im++)
    {  
      fscanf(stdin,"M=%d\n",&m[im]);
      m_max = max(m_max,m[im]);
      fprintf(stdout,"    M = %d\n",m[im]);    
    }
    
    fprintf(stdout,"  Maximum M = %d\n",m_max);    

    /* Read number of nodes specifications. */
    fscanf(stdin,"D=%d\n",&D);
    
    fprintf(stdout,"  D = %d\n",D);    
        
    n_max = 1<<ngpt(m_max);	

    /** Allocate data structures. */
    f_hat = (complex**) malloc((2*m_max+1)*sizeof(complex*));
    f_hat2 = (complex**) malloc((2*m_max+1)*sizeof(complex*));
    f_hatr = (complex**) malloc((2*m_max+1)*sizeof(complex*));
    for (n = -m_max; n <= m_max; n++)
    {
      f_hat[n+m_max] = (complex*) malloc((n_max+1)*sizeof(complex));
      f_hat2[n+m_max] = (complex*) malloc((n_max+1)*sizeof(complex));
      f_hatr[n+m_max] = (complex*) malloc((n_max+1)*sizeof(complex));
    }  
    xi = (double*) malloc(2*D*sizeof(double));
    f = (complex*) malloc(D*sizeof(complex));
    f2 = (complex*) malloc(D*sizeof(complex));
    fr = (complex*) malloc(D*sizeof(complex));
    
    srand48(time(NULL));
		    
    /* Generate random nodes. */
    for (d = 0; d < D; d++)
    {
      xi[2*d] = drand48() - 0.5;
      xi[2*d+1] = 0.5*drand48();
      /*xi[2*d] = (d/(double)d_max)-0.5;
      xi[2*d+1] = 0.25;*/
    }

    /* Generate random Fourier coefficients. */
    for (n = -m_max; n <= m_max; n++)
    {
      for (k = 0; k <= n_max; k++)
      {
        f_hatr[n+m_max][k] = (drand48()-0.5) + I*(drand48()-0.5);
      }  
    }  

    /* Generate random function samples. */
    for (d = 0; d < D; d++)
    {
      fr[d] = (drand48()-0.5) + I*(drand48()-0.5);
    }  
    
    //sprintf(filename_tex,"testcase%d.tex",tc);
    sprintf(filename_dat,"testcase%d.dat",tc);
    //file_tex = fopen(filename_tex,"w");
    file_dat = fopen(filename_dat,"w");
  		fclose(file_dat);
    
    nfsft_precompute(m_max,threshold,0U);
    
    for (im = 0; im < im_max; im++)
    {
      /* Init transform plan. */
      plan = nfsft_init_guru(m[im],D,&f_hat[m_max-m[im]],xi,f,
                             (use_nfft!=0)?(0U):(NFSFT_USE_NDFT),cutoff);
      fprintf(stderr,"      M = %d:\n ",m[im]);
          
      /* Test foward transform. */    
      copyc_hat(&f_hat[m_max-m[im]],&f_hatr[m_max-m[im]],m[im]);      
      ndsft_trafo(plan);
        
      copyc(f2,f,D);
      
      copyc_hat(&f_hat[m_max-m[im]],&f_hatr[m_max-m[im]],m[im]);
      err_inf = norm_f_hat_1(&f_hat[m_max-m[im]],m[im]);
      nfsft_trafo(plan);
      
      updatec_xpay(f, -1.0, f2, D);
      err2 = norm_complex_2(f,D)/norm_complex_2(f2,D);
      err_inf = norm_complex_inf(f,D)/err_inf;
      
      nfsft_finalize(plan);
            
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
									
	     //fprintf(stderr,"||b||_1 = %E, L = %d\n",norm_complex_1(b,ld[ild][0]),ld[ild][0]);

      file_dat = fopen(filename_dat,"a");
          
      /*fprintf(file_tex,"%6d & %6d & %.1E & %.1E & %.1E & %.1E & %.1E\\\\\n",
      ld[ild][0],ld[ild][1],t_d,t_dp,t_fd,t_f,error_complex_inf(f, f_m, 
      ld[ild][1])/norm_complex_1(b,ld[ild][0]));*/

      fprintf(file_dat,"%3d %.4E %.4E\n",m[im],err2,err_inf);
     	fclose(file_dat);
    }
        
    nfsft_forget();
    
    free(f);
    free(f2);
    free(fr);
    free(xi);
    for (n = -m_max; n <= m_max; n++)
    {
      free(f_hat[n+m_max]);
      free(f_hat2[n+m_max]);
      free(f_hatr[n+m_max]);
    }   
    free(f_hat);
    free(f_hat2);
    free(f_hatr);
  }
  
  return EXIT_SUCCESS;
}
