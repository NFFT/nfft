#include "api.h"
#include "util.h"
#include "u.h"
#include "direct.h"
#include "legendre.h"

/* Global structure for precomputed values. */
static struct nfsft_wisdom wisdom = {false};
int nfft_size[2] = {0,0};
int fftw_size[2] = {0,0};

//#define THRESHOLD 1000
//#define THRESHOLD 78.203

int gthreshold;

/**
 * Fast Legendre Function Transform
 *
 * \param 
 */
void flft(struct U_type ***U, int M, int t, int n, complex *f_hat, 
         struct nfsft_transform_wisdom *tw)
{
  /** Level index \f$\tau\f$ */
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
  /** Next greater power of two with respect to M since t=ceil(log2(M)) */
  int N = 1<<t;
  /** N + 1 */
  int N1 = N + 1;
  /** Polynomial array lenbgth for stabilization */
  int plength_stab;
  /** */
  double gamma;
  /** Loop counter */
  int j;
  
  /* Initialize working arrays. */
  memset(tw->work,0U,(N1<<1)*sizeof(complex));
  memset(tw->ergeb,0U,(N1<<1)*sizeof(complex));

  /* Set first n Fourier-coefficients explicitely to zero. They should be 
   * zero right from the beginning anyway. */
  memset(f_hat,0U,n*sizeof(complex));
  
  /* First step */ 
  for (j = 0; j < N; j++) 
  {
    tw->work[2*j] = f_hat[j];
  }
  
  /* Use three-term recurrence to map last coefficient a_N to a_{N-1} and 
   * a_{N-2}. */
  tw->work[2*(N-2)] += wisdom.gamma[ROW(n)+N]*f_hat[N];
  tw->work[2*(N-1)] = f_hat[N-1] + wisdom.beta[ROW(n)+N]*f_hat[N];
  tw->work[2*(N-1)+1] = wisdom.alpha[ROW(n)+N]*f_hat[N];
  
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
        multiplyU(tw->vec3, tw->vec4, act_U, tau, n, plength*l+1, tw, wisdom.gamma[ROWK(n)+plength*l+1-n+1]);        
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
                
        /* Ng = Zweierpotenz um schnelle Polynomultiplikation auszuf¸hren */
        /* Vektor auf die volle Laenge verlaengern */
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
  
  // Copy result.
  for (j = 0; j < 2*N1; j++)
  {
    tw->ergeb[j] += tw->work[j];
  }  
  
  /* The final step */ 
  gamma = wisdom.gamma_m1[n];
  
  /* funktioniert nur fuer diese alphas !!*/
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


/**
 * Fast Legendre Function Transform
 *
 * \param 
 */
void flft_adjoint(struct U_type ***U, int M, int t, int n, complex *f_hat, 
          struct nfsft_transform_wisdom *tw)
{
  int plength, nsteps, firstl, lastl, tau;
  int N, N1;
  double alpha, beta, gamma, gammaconst;
  int i,k,l,j;
  int coeff_index;
  struct U_type act_U;
  int rindex = ROWK(n);
  double *ngamma = &(wisdom.gamma[ROWK(n)]);
  int degree_stab, plength_stab, tau_stab;
  
  /* Calculate auxilliary values. */
  N = pow2(t);
  N1 = N + 1;
  
  /* Initialize working arrays. */
  memset(tw->work,0U,(N1<<1)*sizeof(complex));
  memset(tw->ergeb,0U,(N1<<1)*sizeof(complex));
  
  /* The final step */ 
  gamma = wisdom.gamma_m1[n];
    
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
  alpha = wisdom.alpha[rindex+coeff_index];
  beta = wisdom.beta[rindex+coeff_index];
  gamma = wisdom.gamma[rindex+coeff_index];

  f_hat[N] = gamma*tw->work[2*N-4] + beta*tw->work[2*N-2] + alpha*tw->work[2*N-1];
}


/**
 * Converts Chebyshev coefficients to Fourier coefficients.
 */
void cheb2exp(complex *f_hat, complex **DATA,int M,int N)
{
  int k,n;
  complex *data;
  complex last, act;
  complex *f_hat_p;
  complex *f_hat_n;
  int l1,l2,u1,u2;
  int dim, dimh;
  int rowz;
  int colz;
  dimh = N + 1;
  dim = dimh<<1;
  rowz = (dimh)*dim;
  colz = dimh;

  if (M%2 == 0)
  {
    l1 = -M;
    u1 = M;
    l2 = -M+1;
    u2 = M-1;
  }
  else
  {
    l2 = -M;
    u2 = M;
    l1 = -M+1;
    u1 = M-1;
  }
  
  memset(f_hat,0U,dim*dim*sizeof(complex));
  
  /* Process even terms. */
  for (n = l1; n <= u1; n += 2)
  {
    data = DATA[n+M];
    
    f_hat[rowz+n*dim+colz] = *(data++);
    
    f_hat_p = &f_hat[rowz+n*dim+colz+1];
    f_hat_n = &f_hat[rowz+n*dim+colz-1];
    
    for(k = 1; k <= N; k++)
    {
      *f_hat_p++ = (*(data))/2.0;
      *f_hat_n-- = (*(data++))/2.0;
    }
  }
  
  /* Process odd terms. */
  for (n = l2; n <= u2; n += 2)
  {
    data = DATA[n+M];
    
    f_hat[rowz+n*dim+colz] = *(data++);
    
    f_hat_p = &f_hat[rowz+n*dim+colz+1];
    f_hat_n = &f_hat[rowz+n*dim+colz-1];
    
    for(k = 1; k <= N; k++)
    {
      *f_hat_p++ = (*(data))/2.0;
      *f_hat_n-- = (*(data++))/2.0;
    }
    
    /* Incorporate sine term. */
    last = f_hat[rowz+n*dim+colz-N];
    f_hat[rowz+n*dim+colz-N] = I * f_hat[rowz+n*dim+colz-N+1]/2.0;
    for (k = -N+1; k <= N-1; k++)
    {
      act = f_hat[rowz+n*dim+colz+k];
      f_hat[rowz+n*dim+colz+k] = I * (f_hat[rowz+n*dim+colz+k+1] - last)/2.0;       
      last = act;
    }
    f_hat[rowz+n*dim+colz+N] = - I * last/2.0;
  }  
}


/**
 * Converts Chebyshev coefficients to Fourier coefficients.
 */
void cheb2exp_adjoint(complex *f_hat, complex **DATA,int M,int N)
{
  int k,n;
  complex *data;
  complex last, act;
  complex *f_hat_p;
  complex *f_hat_n;
  int l1,l2,u1,u2;
  int dim, dimh;
  int rowz;
  int colz;
  
  dimh = N + 1;
  dim = dimh<<1;
  rowz = (dimh)*dim;
  colz = dimh;
  
  if (M%2 == 0)
  {
    l1 = -M;
    u1 = M;
    l2 = -M+1;
    u2 = M-1;
  }
  else
  {
    l2 = -M;
    u2 = M;
    l1 = -M+1;
    u1 = M-1;
  }
  
  /* Process even terms. */
  for (n = l1; n <= u1; n += 2)
  {   
    data = DATA[n+M];

    memset(data, 0U, (N+1)*sizeof(complex));
    
    *data++ = f_hat[rowz+n*dim+colz];
    
    f_hat_p = &f_hat[rowz+n*dim+colz+1];
    f_hat_n = &f_hat[rowz+n*dim+colz-1];
    
    for(k = 1; k <= N; k++)
    {
      *data++ = (*f_hat_p++ + *f_hat_n--)/2.0;
    }
  }
  
  /* Process odd terms. */
  for (n = l2; n <= u2; n += 2)
  {       
    /* Incorporate sine term. */
    last = f_hat[rowz+n*dim+colz-N];
    f_hat[rowz+n*dim+colz-N] = I * f_hat[rowz+n*dim+colz-N+1]/2.0;
    for (k = -N+1; k <= N-1; k++)
    {
      act = f_hat[rowz+n*dim+colz+k];
      f_hat[rowz+n*dim+colz+k] = I * (f_hat[rowz+n*dim+colz+k+1] - last)/2.0;       
      last = act;
    }
    f_hat[rowz+n*dim+colz+N] = -I * last/2.0;

    data = DATA[n+M];
    memset(data, 0U, (N+1)*sizeof(complex));
    
    *data++ = f_hat[rowz+n*dim+colz];
    
    f_hat_p = &f_hat[rowz+n*dim+colz+1];
    f_hat_n = &f_hat[rowz+n*dim+colz-1];
    
    for(k = 1; k <= N; k++)
    {
      *data++ = (*f_hat_p++ + *f_hat_n--)/2.0;
    }
  }  
}

// --- Wisdom stuff ---

struct nfsft_transform_wisdom* init_transform_wisdom(int t, double threshold)
{  
  int i;
  int ti;
  int N = 1<<t;
  int N1;
  struct nfsft_transform_wisdom *tw;
  
  N1 = N+1;
  tw = malloc(sizeof(struct nfsft_transform_wisdom));
  wisdom.transform_wisdoms[t] = tw;
  tw->N = N;
  tw->t = t;
  
  tw->work = (complex*) calloc(3*N1,sizeof(complex));   
  tw->old = malloc(2*N*sizeof(complex));
  tw->ergeb = (complex*) calloc(3*N1,sizeof(complex)); /* hier werden schrittweise die Cheb.-koeffizienten aufgebaut */
  tw->vec1 = (complex*) fftw_malloc(sizeof(complex)*N1);
  tw->vec2 = (complex*) fftw_malloc(sizeof(complex)*N1);
  tw->vec3 = (complex*) fftw_malloc(sizeof(complex)*N1);
  tw->vec4 = (complex*) fftw_malloc(sizeof(complex)*N1);
  tw->a2 = (complex*) fftw_malloc(sizeof(complex)*N1);
  tw->b2 = (complex*) fftw_malloc(sizeof(complex)*N1);
  
  tw->kinds = malloc(2*sizeof(fftw_r2r_kind));
  tw->kinds[0] = FFTW_REDFT01;
  tw->kinds[1] = FFTW_REDFT01;
  tw->kindsr = malloc(2*sizeof(fftw_r2r_kind));
  tw->kindsr[0] = FFTW_REDFT10;
  tw->kindsr[1] = FFTW_REDFT10;
  
  tw->plans_dct3 = (fftw_plan*) fftw_malloc(sizeof(fftw_plan)*(t-1));
  tw->plans_dct2 = (fftw_plan*) fftw_malloc(sizeof(fftw_plan)*(t-1));
  tw->lengths = (int *) malloc((t-1)*sizeof(int));
  for(i=0,ti=4;i<t-1;i++,ti<<=1)
  {
    tw->lengths[i] = ti;
    tw->plans_dct3[i] = fftw_plan_many_r2r(1, &tw->lengths[i], 2, (double*)tw->vec1, NULL, 2, 1,
                                           (double*)tw->vec1, NULL, 2, 1, tw->kinds, 0);
    tw->plans_dct2[i] = fftw_plan_many_r2r(1, &tw->lengths[i], 2, (double*)tw->vec1, NULL, 2, 1,
                                           (double*)tw->vec1, NULL, 2, 1, tw->kindsr, 0);
  }

  tw->U = precomputeU(t, threshold, wisdom.alpha, wisdom.beta, wisdom.gamma);
    
  return tw;
}

void forget_transform_wisdom(struct nfsft_transform_wisdom *tw)
{
  int ts,mu;
  int i,i_j,Ni,nleg;
  
  for (nleg = 0; nleg <= tw->N; nleg++)
  {
    if (nleg > 1)
    { 
      mu = (int) (log((double)nleg)/log(2.0));
    }  
    
    for (ts = 1; ts < tw->t; ts++)
	   {
	     if (nleg > 1) 
      {
        i_j = pow2(mu-ts-1);
      }  
	     else 
      {
        i_j = 0;
      }  
      
	     Ni = (int) ceil((double)tw->N/(double)pow2(ts+1)) - i_j;
      
      for (i = 0; i < Ni; i++)
	     {
	       fftw_free(tw->U[nleg][ts][i].m1);
	       fftw_free(tw->U[nleg][ts][i].m2);
	       fftw_free(tw->U[nleg][ts][i].m3);
	       fftw_free(tw->U[nleg][ts][i].m4);
      }
      
      fftw_free(tw->U[nleg][ts]);
    }
    
    fftw_free(tw->U[nleg]);
  }
  
  fftw_free(tw->U);
  
  fftw_free(tw->work);
  fftw_free(tw->old);
  fftw_free(tw->ergeb);
  fftw_free(tw->vec1);
  fftw_free(tw->vec2);
  fftw_free(tw->vec3);
  fftw_free(tw->vec4);
  fftw_free(tw->a2);
  fftw_free(tw->b2);
  
  for(i = 0; i < tw->t-1; i++)
  {
    fftw_destroy_plan(tw->plans_dct3[i]);
    fftw_destroy_plan(tw->plans_dct2[i]);
  }  
  free(tw->plans_dct3);
  free(tw->plans_dct2);
  free(tw->kinds);
  free(tw->kindsr);
  free(tw->lengths);
}

inline void init_wisdom()
{ 
  if (wisdom.initialized == false)
  {
    wisdom.initialized = true;
    wisdom.transform_wisdoms = (struct nfsft_transform_wisdom**) 
      calloc(BWEXP_MAX+1,sizeof(struct nfsft_transform_wisdom*));
    
    alpha_al_all(wisdom.alpha,BW_MAX);
    beta_al_all(wisdom.beta,BW_MAX);
    gamma_al_all(wisdom.gamma,BW_MAX);
    gamma_al_m1_all(wisdom.gamma_m1,BW_MAX);
  }  
}

void nfsft_forget_wisdom()
{
  static int i;
  for (i = 0; i <= BWEXP_MAX; i++)
  {
    if (wisdom.transform_wisdoms[i] != 0)
    {  
      forget_transform_wisdom(wisdom.transform_wisdoms[i]);
    }
  }
  free(wisdom.transform_wisdoms);
  wisdom.initialized = false;
}

/*void export_transform_wisdom(struct nfsft_transform_wisdom *tw, FILE *f)
{
  int ts,mu;
  int i,i_j,Ni,nleg;
  int N,N2;
  
  for (nleg = 0; nleg <= tw->N; nleg++)
  {
    if (nleg > 1)
    { 
      mu = (int) (log2(nleg));
    }  
    
    for (ts = 1; ts < tw->t; ts++)
	   {
	     if (nleg > 1) 
      {
        i_j = pow2(mu-ts-1);
      }  
	     else 
      {
        i_j = 0;
      }  
      
	     Ni = (int) ceil((double)tw->N/(double)pow2(ts+1)) - i_j;
      
      N = pow2(ts);
	     N2 = 2*N;
      
      for (i = 0; i < Ni; i++)
	     {
        fwrite(tw->U[nleg][ts][i].m1,sizeof(double),N2+1,f);
        fwrite(tw->U[nleg][ts][i].m1,sizeof(double),N2+1,f);
        fwrite(tw->U[nleg][ts][i].m1,sizeof(double),N2+1,f);
        fwrite(tw->U[nleg][ts][i].m1,sizeof(double),N2+1,f);
      }
    }
  }
  
  fftw_free(tw->work);
  fftw_free(tw->ergeb);
  fftw_free(tw->vec1);
  fftw_free(tw->vec2);
  fftw_free(tw->vec3);
  fftw_free(tw->vec4);
  
  for(i = 0; i < tw->t-1; i++)
  {
    fftw_destroy_plan(tw->plans_dct3[i]);
    fftw_destroy_plan(tw->plans_dct2[i]);
  }  
  free(tw->plans_dct3);
  free(tw->plans_dct2);
  free(tw->kinds);
  free(tw->kindsr);
  free(tw->lengths);
}

void nfsft_export_wisdom(const char* filename)
{
  int i;
  int nul;
  int bwexp_max = BWEXP_MAX;
  nul = 0;
  FILE *f = fopen(filename,"w");
  fwrite(&bwexp_max,sizeof(int),1,f);
  for (i = 0; i <= BWEXP_MAX; i++)
  {
    if (wisdom.transform_wisdoms[i] != 0)
    {
      export_transform_wisdom(wisdom.transform_wisdoms[i],f);
    }
    else
    {
      fwrite(&nul,sizeof(int),1,f);      
    }
  }  
  fclose(f);
}

void nfsft_import_wisdom(const char* filename)
{
  FILE *f = fopen(filename,"r");
  int max;
  fread(&max,sizeof(int),1,f);
  printf("BWEXP_MAX=%d\n",max);
  fclose(f);
}*/

// --- Public interface ---

void nfsft_compute_wisdom(int M, int threshold)
{
  int N;
  int t;
  t = (int) ceil(log((double)M)/log(2.0));
  /* Calculate N as next greater power of 2 of the bandwidth M. */
  N = pow2(t);  
  init_wisdom();
  if (M > 1)
  {  
    init_transform_wisdom(t,threshold);
  }
  gthreshold = threshold;
}

nfsft_plan nfsft_init(int d, int m, double *angles, fftw_complex **f_hat, 
                      fftw_complex *f, nfsft_flags flags)
{
  nfsft_plan plan = malloc(sizeof(struct nfsft_plan_s));
  
  plan->D = d;
  plan->M = m;
  plan->angles = angles;
  plan->f_hat = f_hat;
  plan->f = f;
  plan->threshold = gthreshold;
  
  return(plan);
}

void nfsft_finalize(nfsft_plan plan)
{
  free(plan);
}

void ndsft_trafo(nfsft_plan plan)
{
  /** Counter for loops */
  int i;
  
  init_wisdom();

  if (plan->M == 0)
  {
    for (i = 0; i < plan->D; i++)
    {
      plan->f[i] = plan->f_hat[0][0];
    }  
  }
  else
  {
    /* Compute direct transform. */
    ndsft(plan->D, plan->angles, plan->f, plan->M, plan->f_hat, &wisdom);
  }  
}


void ndsft_adjoint(nfsft_plan plan)
{
  /** Counter for loops */
  int i;
  
  init_wisdom();
  
  adjoint_ndsft(plan->D, plan->angles, plan->f, plan->M, plan->f_hat, 
                &wisdom);
}


void nfsft_trafo(nfsft_plan plan)
{
  /** Exponent of next greater power of 2 realtive to M */
  int t;
  /** Next greater power of 2 relative to M */
  int N;
  /** Counter for loops */
  int i, n;
  /** */
  struct nfsft_transform_wisdom *tw;
  /** */
  nfft_plan myplan;

  /* Init global structure. */
  init_wisdom();

  if (plan->M < 3)
  {
    ndsft_trafo(plan);
  }
  else
  {
    t = ngpt(plan->M);
    /** Next greater power of 2 relative to M */
    N = 1<<t;

    /* Ensure that precomputation has been done. */
    if (wisdom.transform_wisdoms[t] == 0)
    {
      tw = init_transform_wisdom(t,gthreshold);
    }  
    else
    {
      tw = wisdom.transform_wisdoms[t];
    }  
    
    /* Compute FLFT. */
    for(n = -plan->M, i = 0; n <= plan->M; n++, i++) 
    {
      flft(tw->U, plan->M, t, abs(n), plan->f_hat[i], tw);
    }
   
    /* Compute NFFT. */
    nfft_size[0] = 2*(N+1);
    nfft_size[1] = 2*(N+1);
    fftw_size[0] = 4*N;
    fftw_size[1] = 4*N;
    nfft_init_specific(&myplan, 2, nfft_size, plan->D, fftw_size, 
                       6, PRE_PHI_HUT| PRE_PSI| MALLOC_F_HAT | FFT_OUT_OF_PLACE, 
                       FFTW_ESTIMATE| FFTW_DESTROY_INPUT);
    /* Assign angle array. */
    myplan.x = plan->angles;
    if (myplan.nfft_flags & PRE_PSI) 
    {  
      nfft_precompute_psi(&myplan); 
    } 
    /* Assign result array. */
    myplan.f = plan->f;

    /* Convert Chebyshev coefficients to Fourier coefficients. */
    cheb2exp(myplan.f_hat, plan->f_hat, plan->M, N); 
    
    /* Execute NFFT. */
    nfft_trafo(&myplan);
    nfft_finalize(&myplan);
  }    
}

void nfsft_adjoint(nfsft_plan plan)
{
  /** Exponent of next greater power of 2 realtive to M */
  int t;
  /** Next greater power of 2 relative to M */
  int N;
  /** Counter for loops */
  int i, n;
  /** */
  struct nfsft_transform_wisdom *tw;
  /** */
  nfft_plan myplan;
  
  /* Init global structure. */
  init_wisdom();
  
  if (plan->M < 3)
  {
    ndsft_adjoint(plan);
  }
  else
  {
    t = ngpt(plan->M);
    /** Next greater power of 2 relative to M */
    N = 1<<t;
    
    /* Ensure that precomputation has been done. */
    if (wisdom.transform_wisdoms[t] == 0)
    {
      tw = init_transform_wisdom(t,gthreshold);
    }  
    else
    {
      tw = wisdom.transform_wisdoms[t];
    }  
    
    /* Calculate adjoint NFFT. */
    nfft_size[0] = 2*(N+1);
    nfft_size[1] = 2*(N+1);
    fftw_size[0] = 4*N;
    fftw_size[1] = 4*N;
    nfft_init_specific(&myplan, 2, nfft_size, plan->D, fftw_size, 
                       6, PRE_PHI_HUT| PRE_PSI| MALLOC_F_HAT | FFT_OUT_OF_PLACE, 
                       FFTW_ESTIMATE| FFTW_DESTROY_INPUT);
    /* Assign angle array. */
    myplan.x = plan->angles;
    /* Assign result array. */
    myplan.f = plan->f;
    if (myplan.nfft_flags & PRE_PSI) 
    {  
      nfft_precompute_psi(&myplan); 
    } 
    
    /* Execute adjoint NFFT. */
    nfft_adjoint(&myplan);
    
    /* Convert Chebyshev coefficients to Fourier coefficients. */
    cheb2exp_adjoint(myplan.f_hat, plan->f_hat, plan->M, N); 
    
    nfft_finalize(&myplan);      
    
    /* Calculate FLFT. */
    for(n = -plan->M, i = 0; n <= plan->M; n++, i++) 
    {
      flft_adjoint(tw->U,plan->M,t,abs(n),plan->f_hat[i],tw);
    }          
  }    
}
