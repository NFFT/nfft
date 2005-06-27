/* 
   Healpix - Healpix example for the NFSFT

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
#include "infsft.h"
#include "util.h"
#include "time.h"
#include "api.h"

//#define K_MIN 3
//#define K_MAX 360

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
 	complex *f;
  complex *f_temp;
  complex **f_hat;
  complex **f_hat_orig;
  complex **f_hat_temp;
  complex *f_orig;
  double *x;
  double *x_phi, *x_theta;
  double *x2;
  nfsft_plan plan, plan2;
  int n,m,k,d,j,l,i, M, /*k_min, k_max,*/ d_phi, d_theta, N, n_d;
  FILE *file;
  char filename[100];
  double a,b,c,e;
  int D;
  int D_PHI, D_THETA;
  int D2;
  infsft_plan iplan;
  
  D_PHI = 720;
  D_THETA = 360;
  D2 = D_PHI*D_THETA;
 
  if (argc == 1)
  {
    /*k_min = K_MIN;
    k_max = K_MAX; */
    M = 4;
  }
  else if (argc == 2)
  {
    //sscanf(argv[1],"%d",&k_min);
    //sscanf(argv[2],"%d",&k_max);
    sscanf(argv[1],"%d",&M);
  }
  else
  {
    fprintf(stderr,"Healpix - Healpix example for NFSFT\n");
    fprintf(stderr,"Usage: EGM96 M\n");
    return -1;
  }
  
  M = 360;
  
  /* Read Healpix nodes. */  
  sprintf(filename,"healpix.dat");  
  file = fopen(filename,"r");
  
  if (file != NULL)
  {
    fscanf(file,"%d\n",&n_d);
    
    D = 12*n_d *n_d;
    N = 1<<ngpt(M);
    
    fprintf(stderr,"D = %d\n",D);
    fprintf(stderr,"M = %d\n",M);
    fprintf(stderr,"N = %d\n",N);

    /** Allocate data structures. */
    f = (complex*) calloc(D,sizeof(complex));
    f_temp = (complex*) calloc(D2,sizeof(complex));
    f_orig = (complex*) calloc(D,sizeof(complex));
    x = (double*) calloc(2*D,sizeof(double));
    x_phi = (double*) calloc(D_PHI,sizeof(double));
    x_theta = (double*) calloc(D_THETA,sizeof(double));
    x2 = (double*) calloc(2*D2,sizeof(double));
    f_hat = (complex**) malloc((2*M+1)*sizeof(complex*));
    f_hat_temp = (complex**) malloc((2*M+1)*sizeof(complex*));
    f_hat_orig = (complex**) malloc((2*M+1)*sizeof(complex*));
    for (n = -M; n <= M; n++)
    {
      f_hat[n+M] = (complex*) calloc((N+1),sizeof(complex));
      f_hat_temp[n+M] = (complex*) calloc((N+1),sizeof(complex));
      f_hat_orig[n+M] = (complex*) calloc((N+1),sizeof(complex));
    }  
    
    for (d = 0; d < D_PHI; d++)
    {
      x_phi[d] = (d/((double)D_PHI)) - 0.5;
    }

    for (d = 0; d < D_THETA; d++)
    {
      x_theta[d] = (d/((double)(2*D_THETA-2)));
    }
    
    d = 0;
    for (n = 0; n < D_PHI; n++)
    {
      for (k = 0; k < D_THETA; k++)
      {
        x2[2*d] = x_phi[n];
        x2[2*d+1] = x_theta[k];
        d++;
      }
    }
    
    /* Read Healpix nodes */
    for (d = 0; d < D; d++)
    {
  				fscanf(file,"%lf %lf\n",&x[2*d+1],&x[2*d]);
      //fprintf(stderr,"%e %e\n",x[2*d+1],x[2*d]);
    }
    
    fclose(file);

    /* Read EGM96 data. */
    file = fopen("../egm96/egm96.dat","r");
    while (fscanf(file,"%d %d %le %le %le %le\n",&j,&m,&a,&b,&c,&e) != EOF)
    {
      if (j <= M && j > 2)
      {
        //a *= sqrt(2.0*j+1.0); 
        //b *= sqrt(2.0*j+1.0); 
        f_hat_orig[m+M][j] = 0.5*(a + I*b);
        f_hat_orig[-m+M][j] = 0.5*(a - I*b);
      }
    }
    fclose(file);
    /*for (k = 0; k <= M; k++)
    {
      for (n = -k; n <= k ; n++)
      {
        f_hat_orig[n+M][k] = (k<=32)?drand48():0.0;//drand48() + I * drand48();
      }  
    }*/  
    
    copyc_hat(f_hat,f_hat_orig,M);

    nfsft_precompute(M,1000,0U);
    
    plan = nfsft_init(D,M,x,f_hat,f,NFSFT_NORMALIZED);

    nfsft_trafo(plan);
    
    copyc(f_orig,f,D);

    /* Init guess */
    /*for (n = -k_max; n <= k_max; n++)
    {
      memset(f_hat[n+k_max],0U,(N+1)*sizeof(complex));
    }*/  
    
    iplan = infsft_make_plan();
    infsft_init_guru(iplan, plan, NFSFT_CGNR_E);
    
    for (d = 0; d < D; d++)
    {
      iplan->given_f[d] = f_orig[d];
      //iplan->w[d] = wf[d]/(2*M+2);
    }
    
    for (n = -M; n <= M; n++)
    {
      for (k = abs(n); k <= M; k++)
      {
        iplan->f_hat_iter[n+M][k] = 1/sqrt(2*k+1);
      }
    }
    
    //vpr_c_hat(iplan->f_hat_iter,M,"f_hat_iter");

    /* inverse trafo */
    infsft_before_loop(iplan);
    //vpr_c(f_orig,D,"f_orig");
    //vpr_c(iplan->r_iter,D,"r_iter");

    fprintf(stderr,"D2 = %d\n",D2);
    fprintf(stderr,"M = %d\n",M);
    fflush(stderr);
    plan2 = nfsft_init(D2,M,x2,f_hat_temp,f_temp,NFSFT_NORMALIZED);

    copyc_hat(f_hat_temp,f_hat_orig,M);
    nfsft_trafo(plan2);
    sprintf(filename,"it_orig.dat");  
    file = fopen(filename,"w");
    for (d = 0; d < D2; d++)
    {
      fprintf(file,"%16E\n",cabs(f_temp[d]));
    }
    fclose(file);

    copyc_hat(f_hat_temp,iplan->f_hat_iter,M);
    nfsft_trafo(plan2);
    sprintf(filename,"it%02d.dat",0);  
    file = fopen(filename,"w");
    for (d = 0; d < D2; d++)
    {
      fprintf(file,"%16E\n",cabs(f_temp[d]));
    }
    fclose(file);
    
    fprintf(stderr,"%3d: %.3E %.3E,\n",0,sqrt(iplan->dot_r_iter),error_complex_inf(iplan->r_iter,f_orig,iplan->direct_plan->D));
    for(l=0;l<10;l++)
    { 
      infsft_loop_one_step(iplan);
      //vpr_c(iplan->r_iter,D,"r_iter");

      copyc_hat(f_hat_temp,iplan->f_hat_iter,M);
      nfsft_trafo(plan2);
      sprintf(filename,"it%02d.dat",l+1);  
      file = fopen(filename,"w");
      for (d = 0; d < D2; d++)
      {
        fprintf(file,"%16E\n",cabs(f_temp[d]));
      }
      fclose(file);

      fprintf(stderr,"%3d: %.3E %.3E\n",l+1,sqrt(iplan->dot_r_iter),error_complex_inf(iplan->r_iter,f_orig,iplan->direct_plan->D));
    }
    
    infsft_finalize(iplan);  
    nfsft_finalize(plan);  
    nfsft_finalize(plan2);  
    
    nfsft_forget();
    
    for (n = -M; n <= M; n++)
    {
      free(f_hat[n+M]);
      free(f_hat_temp[n+M]);
      free(f_hat_orig[n+M]);
    }   
    free(f_hat);
    free(f_hat_temp);
    free(f_hat_orig);
    free(x);
    free(f);
    free(f_temp);
    free(f_orig);
  }      
  
  return EXIT_SUCCESS;
}
