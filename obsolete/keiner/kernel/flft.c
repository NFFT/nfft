#include "flft.h"

void flft(const int M, const int t, const int n, complex *f_hat, 
          struct U_type ***U, const struct nfsft_transform_wisdom *tw)
{
  /** Next greater power of two with respect to M since t=ceil(log2(M)) */
  const int N = 1<<t;
  /** Level index tau */
  int tau;
  /** Index of first block at current level */  
  int firstl;
  /** Index of last block at current level */
  int lastl;
  /** Block index l */
  int l;
  /** Length of polynomial coefficient arrays at next level */
  int plength;
  /** Current matrix U_{n,tau,l} */
  struct U_type act_U;
  /** Polynomial array length for stabilization */
  int plength_stab;
  /** */
  double gamma;
  /** Loop counter */
  int j;
  
  /* Initialize working arrays. */
  memset(tw->work,0U,((N+1)<<1)*sizeof(complex));
  memset(tw->ergeb,0U,((N+1)<<1)*sizeof(complex));

  /* Set first n Fourier-coefficients explicitely to zero. */
  memset(f_hat,0U,n*sizeof(complex));
  memset(&(f_hat[M+1]),0U,(N-M)*sizeof(complex));
  
  /* First step */ 
  for (j = 0; j < N; j++) 
  {
    tw->work[2*j] = f_hat[j];
  }
  
  /* Use three-term recurrence to map last coefficient a_N to a_{N-1} and 
   * a_{N-2}. */
  tw->work[2*(N-2)] += tw->gamma[ROW(n)+N]*f_hat[N];
  tw->work[2*(N-1)] = f_hat[N-1] + tw->beta[ROW(n)+N]*f_hat[N];
  tw->work[2*(N-1)+1] = tw->alpha[ROW(n)+N]*f_hat[N];
  
  /* Compute the remaining steps. */
  plength = 4;
  for (tau = 1; tau < t; tau++)
  {    
    /* Compute first l. */
    firstl = FIRST_L;
    /* Compute last l. */
    lastl = LAST_L;
    
    /* Compute the multiplication steps. */
    for (l = firstl; l <= lastl; l++)
    {  
      /* Initialize second half of coefficient arrays with zeros. */
      memset(&tw->vec1[plength/2],0U,(plength/2)*sizeof(complex));
      memset(&tw->vec2[plength/2],0U,(plength/2)*sizeof(complex));
      memset(&tw->vec3[plength/2],0U,(plength/2)*sizeof(complex));
      memset(&tw->vec4[plength/2],0U,(plength/2)*sizeof(complex));
      
      /* Copy coefficients into first half. */
      memcpy(tw->vec1,&(tw->work[(plength/2)*(4*l+0)]),(plength/2)*sizeof(complex));
      memcpy(tw->vec2,&(tw->work[(plength/2)*(4*l+1)]),(plength/2)*sizeof(complex));
      memcpy(tw->vec3,&(tw->work[(plength/2)*(4*l+2)]),(plength/2)*sizeof(complex));
      memcpy(tw->vec4,&(tw->work[(plength/2)*(4*l+3)]),(plength/2)*sizeof(complex));     
      
      /* Get matrix U_{n,tau,l} */
      act_U = U[n][tau][l];
      
      /* Check if step is stable. */
      if (act_U.stable)
      {
        /* Multiply third and fourth polynomial with matrix U. */
        multiplyU(tw->vec3, tw->vec4, act_U, tau, n, plength*l+1, tw, 
                  tw->gamma[ROWK(n)+plength*l+1-n+1]);        
        for (j = 0; j < plength; j++)
        {
          tw->work[plength*2*l+j] = tw->vec1[j] += tw->vec3[j];
          tw->work[plength*(2*l+1)+j] = tw->vec2[j] += tw->vec4[j];
        }          
      }
      else
      {	
        /* Stabilize. */
        plength_stab = pow2(t);
                
        memset(&tw->vec3[plength/2],0U,(plength_stab-plength/2)*sizeof(complex));
        memset(&tw->vec4[plength/2],0U,(plength_stab-plength/2)*sizeof(complex));
        
        multiplyU(tw->vec3, tw->vec4, act_U, t-1, n, 1, tw, 0.0);
        
        for (j = 0; j < plength_stab; j++)
        {
          tw->ergeb[plength_stab+j] += tw->vec4[j];
        }
        
        /* Don't change result. */
        memcpy(&(tw->work[plength*2*l]),tw->vec1,plength*sizeof(complex));
        memcpy(&(tw->work[plength*(2*l+1)]),tw->vec2,plength*sizeof(complex)); 
      }
    }
    plength = plength<<1;
  } 
  
  /* */
  for (j = 0; j < 2*(N+1); j++)
  {
    tw->ergeb[j] += tw->work[j];
  }  
  
  /* The last step */ 
  gamma = tw->gamma_m1[n];
  
  if (n%2 == 0)
  {
    if (n == 0)
    {
      f_hat[0] = gamma*(tw->ergeb[0]+tw->ergeb[N+1]*0.5);
      f_hat[1] = gamma*(tw->ergeb[1]+(tw->ergeb[N]+tw->ergeb[N+2]*0.5));
      f_hat[N-1] = gamma*(tw->ergeb[N-1]+tw->ergeb[N+N-2]*0.5);
      f_hat[N] = gamma*(tw->ergeb[N+N-1]*0.5);
      for (j = 2; j < N-1; j++)
      {
        f_hat[j] = gamma*(tw->ergeb[j]+(tw->ergeb[j+N-1]+tw->ergeb[j+N+1])*0.5);
      } 
    }
    else
    {
      f_hat[0] = gamma*(tw->ergeb[0]+tw->ergeb[N]-tw->ergeb[N+1]*0.5);
      f_hat[1] = gamma*(tw->ergeb[1]+tw->ergeb[N+1]-(tw->ergeb[N]+tw->ergeb[N+2]*0.5));
      f_hat[N-1] = gamma*(tw->ergeb[N-1]+tw->ergeb[N+N-1]-tw->ergeb[N+N-2]*0.5);
      f_hat[N] = gamma*(tw->ergeb[N+N]-tw->ergeb[N+N-1]*0.5);
      for (j = 2; j < N-1; j++)
      {
        f_hat[j] = gamma*(tw->ergeb[j]+tw->ergeb[j+N]-(tw->ergeb[j+N-1]+tw->ergeb[j+N+1])*0.5);
      } 
    }
  }
  else
  {
    f_hat[N] = 0.0;
    for (j = 0; j < N; j++)
    {
      f_hat[j] = /*---*/ - gamma*(tw->ergeb[j]+tw->ergeb[j+N]);
    }
  }
}


void flft_adjoint(const int M, const int t, const int n, complex *f_hat,  
                  struct U_type ***U, const struct nfsft_transform_wisdom *tw)
{
  int plength, /*nsteps, */firstl, lastl, tau;
  int N, N1;
  double alpha, beta, gamma, gammaconst;
  int i,k,l,j;
  int coeff_index;
  struct U_type act_U;
  int rindex = ROWK(n);
  double *ngamma = &(tw->gamma[ROWK(n)]);
  int /*degree_stab, */plength_stab, tau_stab;
  
  /* Calculate auxilliary values. */
  N = pow2(t);
  N1 = N + 1;
  
  /* Initialize working arrays. */
  memset(tw->work,0U,(N1<<1)*sizeof(complex));
  memset(tw->ergeb,0U,(N1<<1)*sizeof(complex));
  
  /* The final step */ 
  gamma = tw->gamma_m1[n];
    
  /* First half consists always of coefficient vector multiplied by I_{N+1}, 
   * i.e. a copy of this vector. */
  for (i = 0; i <= N; i++)
  {
    tw->work[i] = gamma*f_hat[i]; 
  }
  
  /* Distinguish by n for the second half. */
  if (n%2 == 0)
  {
    if (n == 0)
    {
      /* Second half is T_{N+1}^T */
      tw->work[N+1+0] = gamma * f_hat[1];
      for (i = 1; i < N; i++)
      {
        tw->work[N+1+i] = gamma*0.5*(f_hat[i-1] + f_hat[i+1]);
      } 
      tw->work[N+1+N] = 0.5*gamma*f_hat[N-1];     
    }
    else
    {
      /* Second half is I_{N+1} - T_{N+1}^T */
      tw->work[N+1+0] = gamma * (f_hat[0] - f_hat[1]);
      for (i = 1; i < N; i++)
      {
        tw->work[N+1+i] = gamma * (f_hat[i] - 0.5*(f_hat[i-1] + f_hat[i+1]));
      } 
      tw->work[N+1+N] = gamma * (f_hat[N] - 0.5*f_hat[N-1]);     
    }
  }
  else
  {
    /* Second half is I_{N+1} */    
    for (i = 0; i <= N; i++)
    {
      tw->work[N+1+i] = /*---*/ - gamma*f_hat[i];
    }
  }  
  
  memmove(&tw->work[N],&tw->work[N+1],N*sizeof(complex));
  memset(&tw->work[2*N],0U,2*sizeof(complex));
  
  memcpy(tw->old,tw->work,2*N*sizeof(complex));
  
  /* Compute the remaining steps. */
  plength = N;
  
  for (tau = t-1; tau >= 1; tau--)
  {    
    /* Compute first l. */
    //firstl = pow2((int)log2(n)-(n==N?1:0)-tau-1);
    firstl = FIRST_L;    
    /* Compute last l. */
    lastl = LAST_L;
    
    /* Compute the multiplication steps. */
    for (l = firstl; l <= lastl; l++)
    {  
      /* Initialize second half of coefficient arrays with zeros. */
      memcpy(tw->vec3,&(tw->work[(plength/2)*(4*l+0)]),plength*sizeof(complex));
      memcpy(tw->vec4,&(tw->work[(plength/2)*(4*l+2)]),plength*sizeof(complex));     

      memcpy(&tw->work[(plength/2)*(4*l+1)],&(tw->work[(plength/2)*(4*l+2)]),(plength/2)*sizeof(complex));
           
      /* Get matrix U_(2^tau-1)^n() */
      act_U = U[n][tau][l];
      
      /* Check if step is stable. */
      if (act_U.stable)
      {
        /* Multiply third and fourth polynomial with matrix U. */
        gammaconst = ngamma[plength*l+1-n+1];        
        multiplyU_adjoint(tw->vec3, tw->vec4, act_U, tau, n, plength*l+1, tw, gammaconst);
        memcpy(&(tw->vec3[plength/2]),tw->vec4,(plength/2)*sizeof(complex));
        
        for (j = 0; j < plength; j++)
        {
          tw->work[plength*(4*l+2)/2+j] = tw->vec3[j];
        }

      }
      else
      {	
        /* Stabilize. */
        plength_stab = pow2(t);
        tau_stab = t-1;        
        
        memcpy(tw->vec3,tw->old,plength_stab*sizeof(complex));
        memcpy(tw->vec4,&(tw->old[plength_stab]),plength_stab*sizeof(complex));

        multiplyU_adjoint(tw->vec3, tw->vec4, act_U, tau_stab, n, 1, tw, 0.0);
        memcpy(&(tw->vec3[plength/2]),tw->vec4,(plength/2)*sizeof(complex));
        for (j = 0; j < plength; j++)
        {
          tw->work[(plength/2)*(4*l+2)+j] = tw->vec3[j];
        }        
	     }
    }
    plength = plength>>1;    
  }    
  
  /* First step */ 
  memset(f_hat,0U,N1*sizeof(complex));
  for (k = 0; k < N; k++) 
  {
    f_hat[k] = tw->work[2*k];
  }
  
  coeff_index = N - n;
  alpha = tw->alpha[rindex+coeff_index];
  beta = tw->beta[rindex+coeff_index];
  gamma = tw->gamma[rindex+coeff_index];

  f_hat[N] = gamma*tw->work[2*N-4] + beta*tw->work[2*N-2] + alpha*tw->work[2*N-1];
}
