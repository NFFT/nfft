#include "flft.h"
#include "u.h"
#include "util.h"
#include "legendre.h"
#include "api.h"

void flft(const int M, const int t, const int n, complex *const f_hat, 
          struct nfsft_wisdom *const wisdom, int *nstab, int *ntotal)
{
  /** 
   * Next greater power of two with respect to M since 
   * \f$t := \left\lceil \log_2 M \right\rceil\f$ 
   */
  const int N = 1<<t;
  
  /** Level index \f$tau\f$ */
  int tau;
  /** Index of first block at current level */  
  int firstl;
  /** Index of last block at current level */
  int lastl;
  /** Block index \f$l\f$ */
  int l;
  /** Length of polynomial coefficient arrays at next level */
  int plength;
  /** Polynomial array length for stabilization */
  int plength_stab;
  
  /** Multidimensional array of matrices \$U_{n,\tau,l}\$ */
  struct U_type ****const U  = wisdom->U;
  /** Current matrix \f$U_{n,tau,l}\f$ */
  struct U_type act_U;
  
  /** */
  double mygamma;

  /** Loop counter */
  int j;
  
  complex *work_ptr;
  complex *f_hat_ptr;
  
  int ntilde = min(n,N-2);
  int Mtilde = min(M,N-1);
  
  /* Initialize working arrays. */
  memset(wisdom->work,0U,((N+1)<<1)*sizeof(complex));
  memset(wisdom->ergeb,0U,((N+1)<<1)*sizeof(complex));

  /* Set first n Fourier coefficients explicitly to zero. */
  memset(f_hat,0U,n*sizeof(complex));
  memset(&(f_hat[M+1]),0U,(N-M)*sizeof(complex));
  
  work_ptr = wisdom->work;
  f_hat_ptr = f_hat;
  /* First step */ 
  for (j = 0; j < N; j++) 
  {
    *work_ptr = *f_hat_ptr++;
    work_ptr += 2;
  }
  
  /* Use three-term recurrence to map last coefficient a_N to a_{N-1} and 
   * a_{N-2}. */
  wisdom->work[2*(N-2)]   += gamma_al(N-1,n)*f_hat[N];
  wisdom->work[2*(N-1)]   += beta_al(N-1,n)*f_hat[N];
  wisdom->work[2*(N-1)+1]  = alpha_al(N-1,n)*f_hat[N];
  
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
      /* Copy vectors to multiply into working arrays zero-padded to twice the length. */
      memcpy(wisdom->vec3,&(wisdom->work[(plength/2)*(4*l+2)]),(plength/2)*sizeof(complex));
      memcpy(wisdom->vec4,&(wisdom->work[(plength/2)*(4*l+3)]),(plength/2)*sizeof(complex));     
      memset(&wisdom->vec3[plength/2],0U,(plength/2)*sizeof(complex));
      memset(&wisdom->vec4[plength/2],0U,(plength/2)*sizeof(complex));
      
      /* Copy coefficients into first half. */
      memcpy(&(wisdom->work[(plength/2)*(4*l+2)]),&(wisdom->work[(plength/2)*(4*l+1)]),(plength/2)*sizeof(complex));
      memset(&(wisdom->work[(plength/2)*(4*l+1)]),0U,(plength/2)*sizeof(complex));
      memset(&(wisdom->work[(plength/2)*(4*l+3)]),0U,(plength/2)*sizeof(complex));
      
      /* Get matrix U_{n,tau,l} */
      act_U = U[n][tau][l][0];
      
      *ntotal = *ntotal + 1;
      
      /* Check if step is stable. */
      if (act_U.stable)
      {
        /*fprintf("-----\n");
        fprintf("gamma = %f\n",wisdom->gamma[ROWK(n)+plength*l+1-n+1]);
        fprintf("gamma(n-1) = %f\n",gamma_al(n,plength*l));
        fprintf("gamma(n) = %f\n",gamma_al(n,plength*l+1));
        fprintf("gamma(n+1) = %f\n",gamma_al(n,plength*l+2));
        fprintf("-----\n");*/
        
        /* Multiply third and fourth polynomial with matrix U. */
        multiplyU(wisdom->vec3, wisdom->vec4, act_U, tau, n, plength*l+1, wisdom, 
                  gamma_al(plength*l+1,n));
        if (gamma_al(plength*l+1,n) != 0.0)
        {  
          for (j = 0; j < plength; j++)
          {
            wisdom->work[plength*2*l+j] += wisdom->vec3[j];
          }          
        }
        for (j = 0; j < plength; j++)
        {
          wisdom->work[plength*(2*l+1)+j] += wisdom->vec4[j];
        }          
      }
      else
      {	        
        /* Stabilize. */
        *nstab = *nstab + 1;

        /* Get U suitable for current N. */
        if (wisdom->flags & NFSFT_BW_WINDOW)
        {
        }
        else
        {
          act_U = U[n][tau][l][t-tau-1];
        }

        /* The lengh of the polynomials */
        plength_stab = 1<<t;
                
        /* Set rest of vectors explicitely to zero */
        memset(&wisdom->vec3[plength/2],0U,(plength_stab-plength/2)*sizeof(complex));
        memset(&wisdom->vec4[plength/2],0U,(plength_stab-plength/2)*sizeof(complex));
        
        /* Multiply only fourth polynomial with U for stabilization since 
         * gamma_1^n is zero. Add result to final result vector. */
        multiplyU(wisdom->vec3, wisdom->vec4, act_U, t-1, n, 1, wisdom, 0.0);
        for (j = 0; j < plength_stab; j++)
        {
          wisdom->ergeb[plength_stab+j] += wisdom->vec4[j];
        }
      }
    }
    /* Double length of polynomials. */
    plength = plength<<1;
  } 
  
  /* Add the resulting cascade coeffcients to the coeffcients accumulated from 
   * the stabilization steps. */
  for (j = 0; j < 2*(N+1); j++)
  {
    wisdom->ergeb[j] += wisdom->work[j];
  }  
  
  /* The last step. Compute the Chebyshev coeffcients c_k^n from the 
   * polynomials in front of P_0^n and P_1^n. */ 
  mygamma = gamma_al(-1,n);
  if (n%2 == 0)
  {
    if (n == 0)
    {
      f_hat[0] = mygamma*(wisdom->ergeb[0]+wisdom->ergeb[N+1]*0.5);
      f_hat[1] = mygamma*(wisdom->ergeb[1]+(wisdom->ergeb[N]+wisdom->ergeb[N+2]*0.5));
      f_hat[N-1] = mygamma*(wisdom->ergeb[N-1]+wisdom->ergeb[N+N-2]*0.5);
      f_hat[N] = mygamma*(wisdom->ergeb[N+N-1]*0.5);
      for (j = 2; j < N-1; j++)
      {
        f_hat[j] = mygamma*(wisdom->ergeb[j]+(wisdom->ergeb[j+N-1]+wisdom->ergeb[j+N+1])*0.5);
      } 
    }
    else
    {
      f_hat[0] = mygamma*(wisdom->ergeb[0]+wisdom->ergeb[N]-wisdom->ergeb[N+1]*0.5);
      f_hat[1] = mygamma*(wisdom->ergeb[1]+wisdom->ergeb[N+1]-(wisdom->ergeb[N]+wisdom->ergeb[N+2]*0.5));
      f_hat[N-1] = mygamma*(wisdom->ergeb[N-1]+wisdom->ergeb[N+N-1]-wisdom->ergeb[N+N-2]*0.5);
      f_hat[N] = mygamma*(wisdom->ergeb[N+N]-wisdom->ergeb[N+N-1]*0.5);
      for (j = 2; j < N-1; j++)
      {
        f_hat[j] = mygamma*(wisdom->ergeb[j]+wisdom->ergeb[j+N]-(wisdom->ergeb[j+N-1]+wisdom->ergeb[j+N+1])*0.5);
      } 
    }
  }
  else
  {
    f_hat[N] = 0.0;
    for (j = 0; j < N; j++)
    {
      f_hat[j] = /*---*/ - mygamma*(wisdom->ergeb[j]+wisdom->ergeb[j+N]);
    }
  }
}


void flft_adjoint(const int M, const int t, const int n, complex *const f_hat, 
                  struct nfsft_wisdom *const wisdom)
{
  const int N = 1<<t;
  
  /** Level index tau */
  int tau;
  /** Index of first block at current level */  
  int firstl;
  /** Index of last block at current level */
  int lastl;
  /** Block index \f$l\f$ */
  int l;
  /** Length of polynomial coefficient arrays at next level */
  int plength;
  /** Polynomial array length for stabilization */
  int plength_stab;
  
  /** Multidimensional array of matrices \$U_{n,\tau,l}\$ */
  struct U_type ****const U = wisdom->U;
  /** Current matrix U_{n,tau,l} */
  struct U_type act_U;
  
  /** */
  double gamma;
  
  /** Loop counter */
  int j;

  int ntilde = min(n,N-2);;
  int Mtilde = min(M,N-1);
  
  /* Initialize working arrays. */
  memset(wisdom->work,0U,((N+1)<<1)*sizeof(complex));
  memset(wisdom->ergeb,0U,((N+1)<<1)*sizeof(complex));
  
  /* The final step */ 
  gamma = gamma_al(-1,n);//wisdom->gamma[ROW(n)];//gamma_m1[n];
    
  /* First half consists always of coefficient vector multiplied by I_{N+1}, 
   * i.e. a copy of this vector. */
  for (j = 0; j <= N; j++)
  {
    wisdom->work[j] = gamma*f_hat[j]; 
  }
  
  /* Distinguish by n for the second half. */
  if (n%2 == 0)
  {
    if (n == 0)
    {
      /* Second half is T_{N+1}^T */
      wisdom->work[N+0] = gamma * f_hat[1];
      for (j = 1; j < N; j++)
      {
        wisdom->work[N+j] = gamma*0.5*(f_hat[j-1] + f_hat[j+1]);
      } 
      wisdom->work[N+N] = 0.5*gamma*f_hat[N-1];     
    }
    else
    {
      /* Second half is I_{N+1} - T_{N+1}^T */
      wisdom->work[N+0] = gamma * (f_hat[0] - f_hat[1]);
      for (j = 1; j < N; j++)
      {
        wisdom->work[N+j] = gamma * (f_hat[j] - 0.5*(f_hat[j-1] + f_hat[j+1]));
      } 
      wisdom->work[N+N] = gamma * (f_hat[N] - 0.5*f_hat[N-1]);     
    }
  }
  else
  {
    /* Second half is I_{N+1} */    
    for (j = 0; j <= N; j++)
    {
      wisdom->work[N+j] = /*---*/ - gamma*f_hat[j];
    }
  }  
    
  /** Save copy of inpute data for stabilization steps. */
  memcpy(wisdom->old,wisdom->work,2*N*sizeof(complex));
  
  /* Compute the remaining steps. */
  plength = N;
  for (tau = t-1; tau >= 1; tau--)
  {    
    /* Compute first l. */
    firstl = FIRST_L;    
    /* Compute last l. */
    lastl = LAST_L;
    
    /* Compute the multiplication steps. */
    for (l = firstl; l <= lastl; l++)
    {  
      /* Initialize second half of coefficient arrays with zeros. */
      memcpy(wisdom->vec3,&(wisdom->work[(plength/2)*(4*l+0)]),plength*sizeof(complex));
      memcpy(wisdom->vec4,&(wisdom->work[(plength/2)*(4*l+2)]),plength*sizeof(complex));     

      memcpy(&wisdom->work[(plength/2)*(4*l+1)],&(wisdom->work[(plength/2)*(4*l+2)]),(plength/2)*sizeof(complex));
           
      /* Get matrix U_(2^tau-1)^n() */
      act_U = U[n][tau][l][0];
      
      /* Check if step is stable. */
      if (act_U.stable)
      {
        /* Multiply third and fourth polynomial with matrix U. */
        multiplyU_adjoint(wisdom->vec3, wisdom->vec4, act_U, tau, n, plength*l+1, 
                          wisdom,  gamma_al(plength*l+1,n));
        memcpy(&(wisdom->vec3[plength/2]),wisdom->vec4,(plength/2)*sizeof(complex));
        
        for (j = 0; j < plength; j++)
        {
          wisdom->work[plength*(4*l+2)/2+j] = wisdom->vec3[j];
        }

      }
      else
      {	
        /* Stabilize. */
        
        /* Get U suitable for current N. */
        if (wisdom->flags & NFSFT_BW_WINDOW)
        {
        }
        else
        {
          act_U = U[n][tau][l][t-tau-1];
        }
        
        /* Stabilize. */
        plength_stab = 1<<t;
        
        memcpy(wisdom->vec3,wisdom->old,plength_stab*sizeof(complex));
        memcpy(wisdom->vec4,&(wisdom->old[plength_stab]),plength_stab*sizeof(complex));

        multiplyU_adjoint(wisdom->vec3, wisdom->vec4, act_U, t-1, n, 1, wisdom, 0.0);
        memcpy(&(wisdom->vec3[plength/2]),wisdom->vec4,(plength/2)*sizeof(complex));
        for (j = 0; j < plength; j++)
        {
          wisdom->work[(plength/2)*(4*l+2)+j] = wisdom->vec3[j];
        }        
	     }
    }
    /* Half the length of polynomial arrays. */
    plength = plength>>1;    
  }    
  
  /* First step */ 
  memset(f_hat,0U,(N+1)*sizeof(complex));
  for (j = 0; j < N; j++) 
  {
    f_hat[j] = wisdom->work[2*j];
  }  
  f_hat[N] = gamma_al(N-1,n)*wisdom->work[2*N-4] + beta_al(N-1,n)*wisdom->work[2*N-2] +  alpha_al(N-1,n)*wisdom->work[2*N-1];
}


void flft_stab(const int M, const int t, const int n, complex *const f_hat, 
          struct nfsft_wisdom *const wisdom, int *nstab, int *ntotal)
{
  /** 
  * Next greater power of two with respect to M since 
  * \f$t := \left\lceil \log_2 M \right\rceil\f$ 
  */
  const int N = 1<<t;
  
  /** Level index \f$tau\f$ */
  int tau;
  /** Index of first block at current level */  
  int firstl;
  /** Index of last block at current level */
  int lastl;
  /** Block index \f$l\f$ */
  int l;
  /** Length of polynomial coefficient arrays at next level */
  int plength;
  /** Polynomial array length for stabilization */
  int plength_stab;
  
  /** Multidimensional array of matrices \$U_{n,\tau,l}\$ */
  struct U_type ****const U  = wisdom->U;
  /** Current matrix \f$U_{n,tau,l}\f$ */
  struct U_type act_U;
  
  /** */
  double gamma;
  
  /** Loop counter */
  int j;
  
  complex *work_ptr;
  complex *f_hat_ptr;
  
  int ntilde = min(n,N-2);;
  int Mtilde = min(M,N-1);
  
  /* Initialize working arrays. */
  /*memset(wisdom->work,0U,((N+1)<<1)*sizeof(complex));
  memset(wisdom->ergeb,0U,((N+1)<<1)*sizeof(complex));*/
  
  /* Set first n Fourier coefficients explicitly to zero. */
  /*memset(f_hat,0U,n*sizeof(complex));
  memset(&(f_hat[M+1]),0U,(N-M)*sizeof(complex));*/
  
  //work_ptr = wisdom->work;
  //f_hat_ptr = f_hat;
  /* First step */ 
  /*for (j = 0; j < N; j++) 
  {
    *work_ptr = *f_hat_ptr++;
    work_ptr += 2;
  }*/
  
  /* Use three-term recurrence to map last coefficient a_N to a_{N-1} and 
    * a_{N-2}. */
  /*wisdom->work[2*(N-2)]   += wisdom->gamma[ROW(n)+N]*f_hat[N];
  wisdom->work[2*(N-1)]   += wisdom->beta[ROW(n)+N]*f_hat[N];
  wisdom->work[2*(N-1)+1]  = wisdom->alpha[ROW(n)+N]*f_hat[N];*/
  
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
      /* Copy vectors to multiply into working arrays zero-padded to twice the length. */
      /*memcpy(wisdom->vec3,&(wisdom->work[(plength/2)*(4*l+2)]),(plength/2)*sizeof(complex));
      memcpy(wisdom->vec4,&(wisdom->work[(plength/2)*(4*l+3)]),(plength/2)*sizeof(complex));     
      memset(&wisdom->vec3[plength/2],0U,(plength/2)*sizeof(complex));
      memset(&wisdom->vec4[plength/2],0U,(plength/2)*sizeof(complex));*/
      
      /* Copy coefficients into first half. */
      /*memcpy(&(wisdom->work[(plength/2)*(4*l+2)]),&(wisdom->work[(plength/2)*(4*l+1)]),(plength/2)*sizeof(complex));
      memset(&(wisdom->work[(plength/2)*(4*l+1)]),0U,(plength/2)*sizeof(complex));
      memset(&(wisdom->work[(plength/2)*(4*l+3)]),0U,(plength/2)*sizeof(complex));*/
      
      /* Get matrix U_{n,tau,l} */
      act_U = U[n][tau][l][0];
      
      *ntotal = *ntotal + 1;
      
      /* Check if step is stable. */
      if (act_U.stable)
      {
        /* Multiply third and fourth polynomial with matrix U. */
        /*multiplyU(wisdom->vec3, wisdom->vec4, act_U, tau, n, plength*l+1, wisdom, 
                  wisdom->gamma[ROWK(n)+plength*l+1-n+1]);        
        if (wisdom->gamma[ROWK(n)+plength*l+1-n+1] != 0.0)
        {  
          for (j = 0; j < plength; j++)
          {
            wisdom->work[plength*2*l+j] += wisdom->vec3[j];
          }          
        }
        for (j = 0; j < plength; j++)
        {
          wisdom->work[plength*(2*l+1)+j] += wisdom->vec4[j];
        }*/          
      }
      else
      {	
        /* Stabilize. */
        *nstab = *nstab + 1;
        
        /* Get U suitable for current N. */
        //act_U = U[n][tau][l][t-tau-1];
        
        /* The lengh of the polynomials */
        //plength_stab = 1<<t;
        
        /* Set rest of vectors explicitely to zero */
        //memset(&wisdom->vec3[plength/2],0U,(plength_stab-plength/2)*sizeof(complex));
        //memset(&wisdom->vec4[plength/2],0U,(plength_stab-plength/2)*sizeof(complex));
        
        /* Multiply only fourth polynomial with U for stabilization since 
          * gamma_1^n is zero. Add result to final result vector. */
        /*multiplyU(wisdom->vec3, wisdom->vec4, act_U, t-1, n, 1, wisdom, 0.0);
        for (j = 0; j < plength_stab; j++)
        {
          wisdom->ergeb[plength_stab+j] += wisdom->vec4[j];
        }*/
      }
    }
    /* Double length of polynomials. */
    plength = plength<<1;
  } 
  
  /* Add the resulting cascade coeffcients to the coeffcients accumulated from 
    * the stabilization steps. */
  /*for (j = 0; j < 2*(N+1); j++)
  {
    wisdom->ergeb[j] += wisdom->work[j];
  }*/  
  
  /* The last step. Compute the Chebyshev coeffcients c_k^n from the 
    * polynomials in front of P_0^n and P_1^n. */ 
  //gamma = wisdom->gamma[ROW(n)];//_m1[n];
  /*  if (n%2 == 0)
    {
      if (n == 0)
      {
        f_hat[0] = gamma*(wisdom->ergeb[0]+wisdom->ergeb[N+1]*0.5);
        f_hat[1] = gamma*(wisdom->ergeb[1]+(wisdom->ergeb[N]+wisdom->ergeb[N+2]*0.5));
        f_hat[N-1] = gamma*(wisdom->ergeb[N-1]+wisdom->ergeb[N+N-2]*0.5);
        f_hat[N] = gamma*(wisdom->ergeb[N+N-1]*0.5);
        for (j = 2; j < N-1; j++)
        {
          f_hat[j] = gamma*(wisdom->ergeb[j]+(wisdom->ergeb[j+N-1]+wisdom->ergeb[j+N+1])*0.5);
        } 
      }
      else
      {
        f_hat[0] = gamma*(wisdom->ergeb[0]+wisdom->ergeb[N]-wisdom->ergeb[N+1]*0.5);
        f_hat[1] = gamma*(wisdom->ergeb[1]+wisdom->ergeb[N+1]-(wisdom->ergeb[N]+wisdom->ergeb[N+2]*0.5));
        f_hat[N-1] = gamma*(wisdom->ergeb[N-1]+wisdom->ergeb[N+N-1]-wisdom->ergeb[N+N-2]*0.5);
        f_hat[N] = gamma*(wisdom->ergeb[N+N]-wisdom->ergeb[N+N-1]*0.5);
        for (j = 2; j < N-1; j++)
        {
          f_hat[j] = gamma*(wisdom->ergeb[j]+wisdom->ergeb[j+N]-(wisdom->ergeb[j+N-1]+wisdom->ergeb[j+N+1])*0.5);
        } 
      }
    }
    else
    {
      f_hat[N] = 0.0;
      for (j = 0; j < N; j++)
      {*/
        //f_hat[j] = /*---*/ - gamma*(wisdom->ergeb[j]+wisdom->ergeb[j+N]);
      /*}
    }*/
}

