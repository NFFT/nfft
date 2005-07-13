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

#define NO   (0)
#define YES  (1)
#define BOTH (2)

/** Computes the inner product on \f$\mathbb{S}^2\f$. */
inline double innerProduct(const double phi1, const double theta1, 
                           const double phi2, const double theta2)
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
 	return exp(rho*(x-1));
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
	int dl_max_prec;
  
  /** Index variable for testcases. */
  int tc;
  /** The kernel type */
  int kt;
  int cutoff;
  double threshold;
  
  /** Next greater power of two with respect to m_max */
  int n_max;
 	/** 
   * The relative error
   * 
   * \f[
       \frac{\left\|f-f_M\right\|_{\infty}}{\left\|\mathbf{b}\right\|}
   * \f]
   */
  double t_d, t_dp, t_fd, t_f;
  double nfactor;
  double temp;
	double err;
	int precompute = NO;
	
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
  complex *prec;
  /** NFSFT plan */
 	nfsft_plan plan;
  /** adjoint NFSFT plan */
  nfsft_plan plan_adjoint;
  
 	int j,k,n,d,l,use_nfsft,use_nfft,nsymbols;
	int index;
	int rinc;
  
  FILE *file_tex;
  FILE *file_dat;
  FILE *file_gaussian;
  char filename_tex[100];
  char filename_dat[100];
  char filename_gaussian[100];
	double constant;
	 	
  /* Read number of testcases. */
  fscanf(stdin,"testcases=%d\n",&tc_max);
  
  fprintf(stdout,"Number of testcases: %d\n\n",tc_max);
  
  /* Process testcases. */
  for (tc = 0; tc < tc_max; tc++)
  {
    fprintf(stdout,"Testcase %d:\n",tc);

    fscanf(stdin,"nfsft=%d\n",&use_nfsft);
    if (use_nfsft != NO)
  	{
      fprintf(stdout,"  NFSFT = yes\n");
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
  	}
	  else
	  {
      cutoff = 3;
			threshold = 1000000000000.0;
      fprintf(stdout,"  NFSFT = no\n");
	  }

    /* Initialize bandwidth bound. */
    m_max = 0;
    /* Initialize source node bound. */
    l_max = 0;
    /* Initialize target node bound. */
    d_max = 0;
    
    /* Read kernel type. One of KT_ABEL_POISSON, KT_SINGULARITY, KT_LOC_SUPP 
     * or KT_GAUSSIAN. */
    fscanf(stdin,"kernel=%d\n",&kt);

    fprintf(stdout,"  Kernel type: %d\n",kt);    
    
    fscanf(stdin,"parameter_sets=%d\n",&ip_max);
    p = (double**) malloc(ip_max*sizeof(double*));

    fprintf(stdout,"  Parameter sets: %d\n",ip_max);    
    
    switch (kt)
    {
      case KT_ABEL_POISSON:
      case KT_SINGULARITY:
        for (ip = 0; ip < ip_max; ip++)
        {  
          p[ip] = (double*) malloc(1*sizeof(double));
          fscanf(stdin,"h=%lf\n",&p[ip][0]);
          fprintf(stdout,"    h = %lf\n",p[ip][0]);    
        }
			  break;
      case KT_GAUSSIAN:
        for (ip = 0; ip < ip_max; ip++)
        {  
          p[ip] = (double*) malloc(1*sizeof(double));
          fscanf(stdin,"sigma=%lf\n",&p[ip][0]);
          fprintf(stdout,"    sigma = %lf\n",p[ip][0]);    
        }
        break;
      case KT_LOC_SUPP:
        for (ip = 0; ip < ip_max; ip++)
        {  
          p[ip] = (double*) malloc(2*sizeof(double));
          fscanf(stdin,"h=%lf ",&p[ip][0]);
          fscanf(stdin,"lambda=%lf\n",&p[ip][1]);
          fprintf(stdout,"    h = %lf, lambda = %lf\n",p[ip][0],p[ip][1]);    
        }
        break;
    };
    
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
    
    /* Read number of nodes specifications. */
    fscanf(stdin,"node_sets=%d\n",&ild_max);
    ld = (int**) malloc(ild_max*sizeof(int*));
    
    fprintf(stdout,"  Nodes: %d\n",ild_max);    

  	dl_max_prec = 0;
    for (ild = 0; ild < ild_max; ild++)
    {  
      ld[ild] = (int*) malloc(4*sizeof(int));
      fscanf(stdin,"L=%d ",&ld[ild][0]);
      l_max = max(l_max,ld[ild][0]);
      fscanf(stdin,"D=%d ",&ld[ild][1]);
      d_max = max(d_max,ld[ild][1]);
      fscanf(stdin,"compare=%d ",&ld[ild][2]);
      fscanf(stdin,"precomputed=%d\n",&ld[ild][3]);
			if (ld[ild][2] == YES && ld[ild][3] == YES)
			{
			  precompute = YES;
				dl_max_prec = max(dl_max_prec,ld[ild][0]*ld[ild][1]);
			}
      fprintf(stdout,"    L = %d, D = %d, Compare = %d, Precomputed = %d\n",ld[ild][0],ld[ild][1],ld[ild][2],ld[ild][2]);    
    }
    
    fprintf(stdout,"  Maximum M = %d\n",m_max);    
    fprintf(stdout,"  Maximum L = %d\n",l_max);    
    fprintf(stdout,"  Maximum D = %d\n",d_max);    
    
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
    if (precompute == YES)
		{
  		fprintf(stderr,"reached! %d\n",d_max*l_max);
	  	fflush(stderr);
		  prec = (complex*) malloc(dl_max_prec*sizeof(complex));
		}
    
		
    /* Generate random source nodes and weights. */
    for (l = 0; l < l_max; l++)
    {
      b[l] = /*1.0;*/drand48() - 0.5;
      nu[2*l] = /*0.0;*/drand48() - 0.5;
      nu[2*l+1] = /*0.25;*/0.5*drand48();
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

    sprintf(filename_tex,"testcase%d.tex",tc);
    sprintf(filename_dat,"testcase%d.dat",tc);
    file_tex = fopen(filename_tex,"w");
    file_dat = fopen(filename_dat,"w");
    fprintf(file_dat,"kernel=%d\n",kt);
    fprintf(file_dat,"nfsft=%d\n",use_nfsft);
		if (use_nfsft != NO)
		{
      fprintf(file_dat,"nfft=%d\n",use_nfsft);
      fprintf(file_dat,"cutoff=%d\n",cutoff);
      fprintf(file_dat,"threshold=%lf\n",threshold);
		}
    fprintf(file_tex,"\\begin{tabular}{l|l|l|l|l|l|l|l}\n");
    fprintf(file_tex,"$L$ & $D$ & direct algorithm & with precomputation & fast summation, NDSFT & fast summation, NFSFT & error $E_{\\infty}$\\\\\\hline\n");
    fclose(file_tex);
		fclose(file_dat);
    
    nfsft_precompute(m_max,threshold,0U);
    
    for (ip = 0; ip < ip_max; ip++)
    {
      fprintf(stdout,"  Parameter set %d: ",ip);
      switch (kt)
      {
        case KT_ABEL_POISSON:
          fprintf(stdout," h = %lf\n",p[ip][0]);
          break;
        case KT_SINGULARITY:
          fprintf(stdout," h = %lf\n",p[ip][0]);
          break;
        case KT_LOC_SUPP:
          fprintf(stdout," h = %lf, lambda = %lf\n",p[ip][0],p[ip][1]);
          break;
        case KT_GAUSSIAN:
          fprintf(stdout," rho = %lf\n",p[ip][0]);
          break;
      }          
      
      /* Kernel coeffcients up to m_max */
			switch (kt)
			{
				case KT_ABEL_POISSON:
					for (k = 0; k <= m_max; k++)
					{        
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
          sprintf(filename_gaussian,"gaussian%.0f.dat",p[ip][0]);
					fprintf(stderr,"filename = %s\n",filename_gaussian);
          file_gaussian = fopen(filename_gaussian,"r");
				  if (file_gaussian != NULL)
					{
					  fscanf(file_gaussian,"%d\n",&nsymbols);
						for (k = 0; k <= min(nsymbols,m_max); k++)
						{
  					  fscanf(file_gaussian,"%lf\n",&a[k]);
  					  fprintf(stderr,"a[%d] = %E\n",k,a[k]);
						}
						for (k = nsymbols+1; k <= m_max; k++)
						{
						  a[k] = 0.0;
						}
					}
					else
					{
					  fprintf(stderr,"Couldn't open file %s for reading!\n",filename_gaussian);
					}
					break;                
			}
			
      for (k = 0; k <= m_max; k++)
      {
        a[k] *= (2*k+1)/(4*PI);
      }
      
      for (ild = 0; ild < ild_max; ild++)
      {               
        fprintf(stdout,"    L = %d, D = %d\n",ld[ild][0],ld[ild][1]);
        if (ld[ild][2] != NO)
        {
				  /* Check if direct algorithm with precomputation should be tested. */
  				if (ld[ild][3] != NO)
	  			{
            for (d = 0; d < ld[ild][1]; d++)
            {
              for (l = 0; l < ld[ild][0]; l++)
              {
                temp = innerProduct(2*PI*nu[2*l],2*PI*nu[2*l+1],2*PI*xi[2*d],2*PI*xi[2*d+1]);
								switch (kt)
								{
									case KT_ABEL_POISSON:
										prec[d*l_max+l] = poissonKernel(temp,p[ip][0]);
										break;
									case KT_SINGULARITY:
										prec[d*l_max+l] = singularityKernel(temp,p[ip][0]);
										break;
									case KT_LOC_SUPP:
										prec[d*l_max+l] = locSuppKernel(temp,p[ip][0],p[ip][1]);
										break;  
									case KT_GAUSSIAN:	              
										prec[d*l_max+l] = gaussianKernel(temp,p[ip][0]);
										break;
								}
							}
						}
						rinc = l_max-ld[ild][0];
						index = 0;
						if (kt == KT_LOC_SUPP)
						{
						  constant = ((p[ip][1]+1)/(2*PI*pow(1-p[ip][0],p[ip][1]+1)));
              t_dp = second();
              for (d = 0; d < ld[ild][1]; d++)
              {
                f[d] = 0.0;
                for (l = 0; l < ld[ild][0]; l++)
					  		{
						  	  f[d] += b[l]*prec[index++];
							  }
                f[d] *= constant;
							  index += rinc;
							}
              t_dp = second() - t_dp;
						}
						else
						{
              t_dp = second();
              for (d = 0; d < ld[ild][1]; d++)
              {
                f[d] = 0.0;
                for (l = 0; l < ld[ild][0]; l++)
					  		{
						  	  f[d] += b[l]*prec[index++];
							  }
							  index += rinc;
							}
              t_dp = second() - t_dp;
						}
						printf("t_dp = %f\n",t_dp);
		  		}
					else
					{
					  t_dp = -1.0;
					}
					
					switch (kt)
					{
						case KT_ABEL_POISSON:
							t_d = second();
							for (d = 0; d < ld[ild][1]; d++)
							{
								f[d] = 0.0;
								for (l = 0; l < ld[ild][0]; l++)
								{
									temp = innerProduct(2*PI*nu[2*l],2*PI*nu[2*l+1],2*PI*xi[2*d],2*PI*xi[2*d+1]);
									f[d] += b[l]*poissonKernel(temp,p[ip][0]);
								}
							} 
							t_d = second() - t_d;
							break;
						case KT_SINGULARITY:
							t_d = second();
							for (d = 0; d < ld[ild][1]; d++)
							{
								f[d] = 0.0;
								for (l = 0; l < ld[ild][0]; l++)
								{
									temp = innerProduct(2*PI*nu[2*l],2*PI*nu[2*l+1],2*PI*xi[2*d],2*PI*xi[2*d+1]);
    							f[d] += b[l]*singularityKernel(temp,p[ip][0]);
								}
							} 
							t_d = second() - t_d;
							break;
						case KT_LOC_SUPP:
						  constant = ((p[ip][1]+1)/(2*PI*pow(1-p[ip][0],p[ip][1]+1)));
							t_d = second();
							for (d = 0; d < ld[ild][1]; d++)
							{
								f[d] = 0.0;
								for (l = 0; l < ld[ild][0]; l++)
								{
									temp = innerProduct(2*PI*nu[2*l],2*PI*nu[2*l+1],2*PI*xi[2*d],2*PI*xi[2*d+1]);
    							f[d] += b[l]*locSuppKernel(temp,p[ip][0],p[ip][1]);
								}
                f[d] *= constant;
							} 
							t_d = second() - t_d;
							break;  
						case KT_GAUSSIAN:	              
							t_d = second();
							for (d = 0; d < ld[ild][1]; d++)
							{
								f[d] = 0.0;
								for (l = 0; l < ld[ild][0]; l++)
								{
									temp = innerProduct(2*PI*nu[2*l],2*PI*nu[2*l+1],2*PI*xi[2*d],2*PI*xi[2*d+1]);
    							f[d] += b[l]*gaussianKernel(temp,p[ip][0]);
								}
							} 
							t_d = second() - t_d;
						  break;
					}  

  				printf("t_d = %f\n",t_d);
        }
				else
				{
				  t_d = -1.0;
					t_dp = -1.0;
				}
        
        
        for (im = 0; im < im_max; im++)
        {
          fprintf(stderr,"      M = %d: ",m[im]);
          
          /* Init transform plans. */
          plan_adjoint = nfsft_init_guru(m[im],ld[ild][0],f_hat,nu,b,(use_nfft!=0)?(0U):(NFSFT_USE_NDFT),cutoff);
          plan = nfsft_init_guru(m[im],ld[ild][1],f_hat,xi,f_m,(use_nfft!=0)?(0U):(NFSFT_USE_NDFT),cutoff);
          
					if (use_nfsft == BOTH)
					{
            t_fd = second();
            ndsft_adjoint(plan_adjoint);
						/* Multiplication with diagonal matrix. */
						for (k = 0; k <= m[im]; k++)
						{
							for (n = -k; n <= k; n++)
							{
								f_hat[n+m[im]][k] *= a[k];
							}
						}
            ndsft_trafo(plan);
            t_fd = second() - t_fd;
					}
					else
					{
					  t_fd = -1.0;
					}
					
          /* Adjoint transform */
          t_f = second();
 		      if (use_nfsft != NO)
					{
            nfsft_adjoint(plan_adjoint);
					}
					else
					{
            ndsft_adjoint(plan_adjoint);
					}
            
          /* Multiplication with diagonal matrix. */
          for (k = 0; k <= m[im]; k++)
          {
            for (n = -k; n <= k; n++)
            {
              f_hat[n+m[im]][k] *= a[k];
            }
          }
            
          /* Forward transform */
	 				if (use_nfsft != 0)
				  {
            nfsft_trafo(plan);
					}
					else
					{
            ndsft_trafo(plan);
					}
          t_f = second() - t_f;

          /* Finalize plans */
          nfsft_finalize(plan_adjoint);
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
					
					if (ld[ild][2] != NO)
					{
            err = error_complex_inf(f, f_m, ld[ild][1])/norm_complex_1(b,ld[ild][0]);
  					fprintf(stderr,"err = %E\n",error_complex_inf(f, f_m, ld[ild][1])/norm_complex_1(b,ld[ild][0]));
					}
					else
					{
					  err = -1.0;
					}
					
	        //fprintf(stderr,"||b||_1 = %E, L = %d\n",norm_complex_1(b,ld[ild][0]),ld[ild][0]);

          file_tex = fopen(filename_tex,"a");
          file_dat = fopen(filename_dat,"a");
          
          fprintf(file_tex,"%6d & %6d & %.1E & %.1E & %.1E & %.1E & %.1E\\\\\n",ld[ild][0],ld[ild][1],t_d,t_dp,t_fd,t_f,error_complex_inf(f, f_m, ld[ild][1])/norm_complex_1(b,ld[ild][0]));
          if (ld[ild][2] != 0)
          {
            if (kt == KT_ABEL_POISSON || kt == KT_SINGULARITY || kt == KT_GAUSSIAN)
            {
              fprintf(file_dat,"%.2f %3d %6d %6d %5.2f %5.2f %.4E\n",p[ip][0],m[im],ld[ild][0],ld[ild][1],t_d,t_f,error_complex_inf(f, f_m, ld[ild][1])/norm_complex_1(b,ld[ild][0]));
            }
            else if (kt == KT_LOC_SUPP)
            {
              fprintf(file_dat,"%.2f %2f %3d %6d %6d %5.2f %5.2f %.4E\n",p[ip][0],p[ip][1],m[im],ld[ild][0],ld[ild][1],t_d,t_f,error_complex_inf(f, f_m, ld[ild][1])/norm_complex_1(b,ld[ild][0]));
            }
          }
          else
          {
            if (kt == KT_ABEL_POISSON || kt == KT_SINGULARITY || kt == KT_GAUSSIAN)
            {
              fprintf(file_dat,"%.2f %3d %6d %6d -1.0 %5.2f -1E0\n",p[ip][0],m[im],ld[ild][0],ld[ild][1],t_f);
            }
            else if (kt == KT_LOC_SUPP)
            {
              fprintf(file_dat,"%.2f %2f %3d %6d %6d -1.0 %5.2f -1E0 -1E0\n",p[ip][0],p[ip][1],m[im],ld[ild][0],ld[ild][1],t_f);
            }
          }
  		    fclose(file_tex);
	      	fclose(file_dat);
        }
      }
      fprintf(file_tex,"\\hline\n");
    }
    
    file_tex = fopen(filename_tex,"a");
    file_dat = fopen(filename_dat,"a");
    fprintf(file_tex,"\\end{tabular}\n");
    fclose(file_tex);
    fclose(file_dat);
    
    nfsft_forget();
    
    if (precompute == YES)
		{
		  free(prec);
		}
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
