#include "short_dim_2d.h"

/* computes a 2d ndft by 1d nfft along the dimension 0 times
   1d ndft along dimension 1
*/
void short_nfft_trafo(nfft_plan* this_plan, nfft_plan* plan_1d)
{
  int j,k0;
  double omega,sinx,cosx;
  fftw_complex* save_f_hat;

  save_f_hat=plan_1d->f_hat;

  for(j=0;j<this_plan->M;j++)
    {
      this_plan->f[j][0]= 0.0;
      this_plan->f[j][1]= 0.0;
      
      plan_1d->x[j] = this_plan->x[this_plan->d * j + 1]; 
    }

  for(k0=0;k0<this_plan->N[0];k0++) /* for shorties */
    {
      plan_1d->f_hat = this_plan->f_hat + k0*this_plan->N[1];
   
      nfft_trafo(plan_1d);
      
      for(j=0;j<this_plan->M;j++)
	{
	  omega = ((double)(k0 - this_plan->N[0]/2)) * this_plan->x[this_plan->d * j + 0];
	  cosx = cos(2*PI*omega); sinx = sin(2*PI*omega);

	  this_plan->f[j][0] += (plan_1d->f[j][0]*cosx + plan_1d->f[j][1]*sinx);
	  this_plan->f[j][1] += (plan_1d->f[j][1]*cosx - plan_1d->f[j][0]*sinx);
	}
    }
  
 plan_1d->f_hat=save_f_hat;
}

/* computes a 2d ndft^H by 1d nfft^H along the dimension 0 times
   1d ndft^H along dimension 1
*/
void short_nfft_adjoint(nfft_plan* this_plan, nfft_plan* plan_1d)
{
  int j,k0;
  double omega,sinx,cosx;
  fftw_complex* save_f_hat;
 
  save_f_hat=plan_1d->f_hat;

  for(j=0;j<this_plan->M;j++)
    {     
      plan_1d->x[j] = this_plan->x[this_plan->d * j + 1]; 
    }

  for(k0=0;k0<this_plan->N[0];k0++) /* for shorties */
    {
      //printf("%d\n",k0); fflush(stdout);
      for(j=0;j<this_plan->M;j++)
	{
	  omega = ((double)(k0 - this_plan->N[0]/2)) * this_plan->x[this_plan->d * j + 0];
	  cosx = cos(2*PI*omega); sinx = sin(2*PI*omega);
	  
	  plan_1d->f[j][0] = (this_plan->f[j][0]*cosx - this_plan->f[j][1]*sinx);
	  plan_1d->f[j][1] = (this_plan->f[j][1]*cosx + this_plan->f[j][0]*sinx);
	}
      //printf("%d\n",k0); fflush(stdout);

      plan_1d->f_hat = this_plan->f_hat + k0*this_plan->N[1];

      //vpr_c(plan_1d->f_hat,10,"f_hat");
     
      nfft_adjoint(plan_1d);
      //vpr_c(plan_1d->f_hat,10,"f_hat");
    }

  plan_1d->f_hat=save_f_hat;
}


/* computes a 3d ndft by 1d nfft along the dimension 0 times
   2d ndft along dimension 1,2
*/
void short_nfft_trafo_3d_1(nfft_plan* this_plan, nfft_plan* plan_1d)
{
  int j,k0;
  double omega,sinx,cosx;
  fftw_complex* save_f_hat;

  save_f_hat=plan_1d->f_hat;

  for(j=0;j<this_plan->M;j++)
    {
      this_plan->f[j][0]= 0.0;
      this_plan->f[j][1]= 0.0;
      
      plan_1d->x[j] = this_plan->x[this_plan->d * j + 1]; 
    }

  for(k0=0;k0<this_plan->N[0];k0++) /* for shorties */
    {
      plan_1d->f_hat = this_plan->f_hat + k0*this_plan->N[1];
   
      nfft_trafo(plan_1d);
      
      for(j=0;j<this_plan->M;j++)
	{
	  omega = ((double)(k0 - this_plan->N[0]/2)) * this_plan->x[this_plan->d * j + 0] + 
	    ((double)(k0 - this_plan->N[0]/2)) * this_plan->x[this_plan->d * j + 0];
	  cosx = cos(2*PI*omega); sinx = sin(2*PI*omega);

	  this_plan->f[j][0] += (plan_1d->f[j][0]*cosx + plan_1d->f[j][1]*sinx);
	  this_plan->f[j][1] += (plan_1d->f[j][1]*cosx - plan_1d->f[j][0]*sinx);
	}
    }
  
 plan_1d->f_hat=save_f_hat;
}

/* computes a 3d ndft by 2d nfft along the dimension 0,1 times
   1d ndft along dimension 2
*/
void short_nfft_trafo_3d_2(nfft_plan* this_plan, nfft_plan* plan_1d)
{
  int j,k0;
  double omega,sinx,cosx;
  fftw_complex* save_f_hat;

  save_f_hat=plan_1d->f_hat;

  for(j=0;j<this_plan->M;j++)
    {
      this_plan->f[j][0]= 0.0;
      this_plan->f[j][1]= 0.0;
      
      plan_1d->x[j] = this_plan->x[this_plan->d * j + 1]; 
    }

  for(k0=0;k0<this_plan->N[0];k0++) /* for shorties */
    {
      plan_1d->f_hat = this_plan->f_hat + k0*this_plan->N[1];
   
      nfft_trafo(plan_1d);
      
      for(j=0;j<this_plan->M;j++)
	{
	  omega = ((double)(k0 - this_plan->N[0]/2)) * this_plan->x[this_plan->d * j + 0];
	  cosx = cos(2*PI*omega); sinx = sin(2*PI*omega);

	  this_plan->f[j][0] += (plan_1d->f[j][0]*cosx + plan_1d->f[j][1]*sinx);
	  this_plan->f[j][1] += (plan_1d->f[j][1]*cosx - plan_1d->f[j][0]*sinx);
	}
    }
  
 plan_1d->f_hat=save_f_hat;
}






