#include "short_dim_2d.h"

/* computes a 2d ndft by 1d nfft along the dimension 1 times
   1d ndft along dimension 0
*/
void short_nfft_trafo(nfft_plan* this_plan, nfft_plan* plan_1d)
{
  int j,k0;
  double omega,sinx,cosx;

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
}

/* computes a 2d ndft by 1d nfft along the dimension 1 times
   1d ndft along dimension 0
*/
void short_nfft_trafo_horner(nfft_plan* this_plan, nfft_plan* plan_1d, fftw_complex* exp_omega1)
{
  int j,k0;
  double omega,sinx,cosx,temp;

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
	  temp = this_plan->f[j][0];
	  this_plan->f[j][0] = exp_omega1[j][0]*temp - exp_omega1[j][1]*this_plan->f[j][1];
	  this_plan->f[j][1] = exp_omega1[j][1]*temp + exp_omega1[j][0]*this_plan->f[j][1];

	  this_plan->f[j][0] += plan_1d->f[j][0];
	  this_plan->f[j][1] += plan_1d->f[j][1];
	}
    }

  for(j=0;j<this_plan->M;j++)
    { 
      omega = ((double)(1 - ceil(this_plan->N[0]/2.0))) * this_plan->x[this_plan->d * j + 0];
      cosx = cos(2*PI*omega); sinx = sin(2*PI*omega);

      temp=this_plan->f[j][0];
      this_plan->f[j][0] = cosx*temp - sinx*this_plan->f[j][1];
      this_plan->f[j][1] = sinx*temp + cosx*this_plan->f[j][1];
    }
}

/* computes a 3d ndft by 1d nfft along the dimension 2 times
   2d ndft along dimension 0,1
*/
void short_nfft_trafo_3d_1(nfft_plan* this_plan, nfft_plan* plan_1d)
{
  int j,k0,k1;
  double omega,sinx,cosx;

  for(j=0;j<this_plan->M;j++)
    {
      this_plan->f[j][0]= 0.0;
      this_plan->f[j][1]= 0.0;
      
      plan_1d->x[j] = this_plan->x[this_plan->d * j + 2]; 
    }

  for(k0=0;k0<this_plan->N[0];k0++) /* for shorties */
    for(k1=0;k1<this_plan->N[1];k1++) 
      {
	plan_1d->f_hat = this_plan->f_hat + (k0*this_plan->N[1]+k1)*this_plan->N[2];
	
	nfft_trafo(plan_1d);
	
	for(j=0;j<this_plan->M;j++)
	  {
	    omega = ((double)(k0 - this_plan->N[0]/2)) * this_plan->x[this_plan->d * j + 0]
	      +     ((double)(k1 - this_plan->N[1]/2)) * this_plan->x[this_plan->d * j + 1];
	    
	    cosx = cos(2*PI*omega); sinx = sin(2*PI*omega);
	    
	    this_plan->f[j][0] += (plan_1d->f[j][0]*cosx + plan_1d->f[j][1]*sinx);
	    this_plan->f[j][1] += (plan_1d->f[j][1]*cosx - plan_1d->f[j][0]*sinx);
	  }
      }
}

/* computes a 3d ndft by 1d nfft along the dimension 2 times
   2d ndft along dimension 0,1
*/
void short_nfft_trafo_3d_1_horner(nfft_plan* this_plan, nfft_plan* plan_1d)
{
  int j,k0,k1;
  double omega,sinx,cosx;

  for(j=0;j<this_plan->M;j++)
    {
      this_plan->f[j][0]= 0.0;
      this_plan->f[j][1]= 0.0;
      
      plan_1d->x[j] = this_plan->x[this_plan->d * j + 2]; 
    }

  for(k0=0;k0<this_plan->N[0];k0++) /* for shorties */
    for(k1=0;k1<this_plan->N[1];k1++) 
      {
	plan_1d->f_hat = this_plan->f_hat + (k0*this_plan->N[1]+k1)*this_plan->N[2];
	
	nfft_trafo(plan_1d);
	
	for(j=0;j<this_plan->M;j++)
	  {
	    omega = ((double)(k0 - this_plan->N[0]/2)) * this_plan->x[this_plan->d * j + 0]
	      +     ((double)(k1 - this_plan->N[1]/2)) * this_plan->x[this_plan->d * j + 1];
	    
	    cosx = cos(2*PI*omega); sinx = sin(2*PI*omega);
	    
	    this_plan->f[j][0] += (plan_1d->f[j][0]*cosx + plan_1d->f[j][1]*sinx);
	    this_plan->f[j][1] += (plan_1d->f[j][1]*cosx - plan_1d->f[j][0]*sinx);
	  }
      }
}

/* computes a 3d ndft by 2d nfft along the dimension 1,2 times
   1d ndft along dimension 0
*/
void short_nfft_trafo_3d_2(nfft_plan* this_plan, nfft_plan* plan_2d)
{
  int j,k0;
  double omega,sinx,cosx;

  for(j=0;j<this_plan->M;j++)
    {
      this_plan->f[j][0]= 0.0;
      this_plan->f[j][1]= 0.0;
      
      plan_2d->x[2*j+0] = this_plan->x[this_plan->d * j + 1];
      plan_2d->x[2*j+1] = this_plan->x[this_plan->d * j + 2];
    }

  for(k0=0;k0<this_plan->N[0];k0++) /* for shorties */
    {
      plan_2d->f_hat = this_plan->f_hat + k0*this_plan->N[1]*this_plan->N[2];
   
      nfft_trafo(plan_2d);
      
      for(j=0;j<this_plan->M;j++)
	{
	  omega = ((double)(k0 - this_plan->N[0]/2)) * this_plan->x[this_plan->d * j + 0];
	  cosx = cos(2*PI*omega); sinx = sin(2*PI*omega);

	  this_plan->f[j][0] += (plan_2d->f[j][0]*cosx + plan_2d->f[j][1]*sinx);
	  this_plan->f[j][1] += (plan_2d->f[j][1]*cosx - plan_2d->f[j][0]*sinx);
	}
    }
}






