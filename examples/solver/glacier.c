#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "util.h"
#include "nfft3.h"

double my_weight(double z,double a,double b,double c)
{
    return pow(0.25-z*z,b)/(c+pow(fabs(z),2*a));
}

void glacier(int N,int M)
{  
  int j,k,k0,k1,l,my_N[2],my_n[2];
  double tmp_y;
  nfft_plan p;
  infft_plan ip;
  FILE* fp;
   
  /* initialise p */
  my_N[0]=N; my_n[0]=nfft_next_power_of_2(N);
  my_N[1]=N; my_n[1]=nfft_next_power_of_2(N);
  nfft_init_guru(&p, 2, my_N, M, my_n, 6, 
		 PRE_PHI_HUT| PRE_FULL_PSI|
		 MALLOC_X| MALLOC_F_HAT| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /* initialise ip, specific */
  infft_init_advanced(&ip,&p, CGNE| PRECOMPUTE_DAMP);

  /* init nodes */
  fp=fopen("input_data.dat","r");
  for(j=0;j<p.M_total;j++)
  {
      fscanf(fp,"%le %le %le",&p.x[2*j+0],&p.x[2*j+1],&tmp_y);
      ip.y[j]=tmp_y;
  }
  fclose(fp);
  
  /* precompute psi */
  if(p.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&p);
  
  /* initialise damping factors */
  if(ip.flags & PRECOMPUTE_DAMP)
    for(k0=0;k0<p.N[0];k0++)
      for(k1=0;k1<p.N[1];k1++)
        ip.w_hat[k0*p.N[1]+k1]=
	    my_weight(((double)(k0-p.N[0]/2))/p.N[0],0.5,3,0.001)*
	    my_weight(((double)(k1-p.N[1]/2))/p.N[1],0.5,3,0.001);
  
  /* init some guess */
  for(k=0;k<p.N_total;k++)
      ip.f_hat_iter[k]=0; 

  /* inverse trafo */  
  infft_before_loop(&ip);
  for(l=0;l<40;l++)
    { 
      fprintf(stderr,"Residual ||r||=%e,\n",sqrt(ip.dot_r_iter));
      infft_loop_one_step(&ip);
    }

  for(k=0;k<p.N_total;k++)
    printf("%le %le\n",creal(ip.f_hat_iter[k]),cimag(ip.f_hat_iter[k]));

  infft_finalize(&ip);  
  nfft_finalize(&p);  
}



int main(int argc, char **argv)
{
  glacier(atoi(argv[1]),atoi(argv[2]));

  return 1;
}
