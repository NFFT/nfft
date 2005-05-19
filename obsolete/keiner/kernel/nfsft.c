#include "api.h"
#include "util.h"
#include "u.h"
#include "direct.h"
#include "legendre.h"
#include "flft.h"

/* Global structure for precomputed values. */
static struct nfsft_wisdom wisdom = {false,false};
int nfft_size[2] = {0,0};
int fftw_size[2] = {0,0};

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

/*struct nfsft_transform_wisdom* init_transform_wisdom(int t, double threshold)
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
  tw->ergeb = (complex*) calloc(3*N1,sizeof(complex));
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
  
  tw->alpha = wisdom.alpha;
  tw->beta = wisdom.beta;
  tw->gamma = wisdom.gamma;
  tw->gamma_m1 = wisdom.gamma_m1;
    
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
}*/

inline void precompute_coeffs()
{ 
  if (wisdom.initialized_coeffs == false)
  {
    alpha_al_all(wisdom.alpha,BW_MAX);
    beta_al_all(wisdom.beta,BW_MAX);
    gamma_al_all(wisdom.gamma,BW_MAX);
    gamma_al_m1_all(wisdom.gamma_m1,BW_MAX);
    wisdom.initialized_coeffs = true;    
  }  
}

void nfsft_forget()
{
  int i;
  
  if (wisdom.initialized == true)
  {
    if (wisdom.N >= 4)
    {  
      forgetU(wisdom.U,wisdom.N,wisdom.t);
      
      fftw_free(wisdom.work);
      fftw_free(wisdom.old);
      fftw_free(wisdom.ergeb);
      fftw_free(wisdom.vec1);
      fftw_free(wisdom.vec2);
      fftw_free(wisdom.vec3);
      fftw_free(wisdom.vec4);
      fftw_free(wisdom.a2);
      fftw_free(wisdom.b2);
      
      for(i = 0; i < wisdom.t-1; i++)
      {
        fftw_destroy_plan(wisdom.plans_dct3[i]);
        fftw_destroy_plan(wisdom.plans_dct2[i]);
      }  
      
      free(wisdom.plans_dct3);
      free(wisdom.plans_dct2);
      free(wisdom.kinds);
      free(wisdom.kindsr);
      free(wisdom.lengths);
    }
    wisdom.initialized = false;
  }
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

void nfsft_precompute(int M, int threshold)
{
  int i;
  int ti;

  if (wisdom.initialized == true)
  {
    if (wisdom.N < M || wisdom.threshold != threshold)
    {
      nfsft_forget();
    }
    else
    {
      return;
    }
  }
  
  precompute_coeffs();

  wisdom.t = (int) ceil(log((double)M)/log(2.0));
  wisdom.N = pow2(wisdom.t);

  if (wisdom.N >= 4)
  {  
    wisdom.threshold = threshold;
    wisdom.work  = (complex*) calloc(3*(wisdom.N+1),sizeof(complex));   
    wisdom.old   = (complex*) malloc(2*(wisdom.N)*sizeof(complex));
    wisdom.ergeb = (complex*) calloc(3*(wisdom.N+1),sizeof(complex));
    wisdom.vec1  = (complex*) fftw_malloc(sizeof(complex)*(wisdom.N+1));
    wisdom.vec2  = (complex*) fftw_malloc(sizeof(complex)*(wisdom.N+1));
    wisdom.vec3  = (complex*) fftw_malloc(sizeof(complex)*(wisdom.N+1));
    wisdom.vec4  = (complex*) fftw_malloc(sizeof(complex)*(wisdom.N+1));
    wisdom.a2    = (complex*) fftw_malloc(sizeof(complex)*(wisdom.N+1));
    wisdom.b2    = (complex*) fftw_malloc(sizeof(complex)*(wisdom.N+1));
    wisdom.kinds     = (fftw_r2r_kind*) malloc(2*sizeof(fftw_r2r_kind));
    wisdom.kinds[0]  = FFTW_REDFT01;
    wisdom.kinds[1]  = FFTW_REDFT01;
    wisdom.kindsr    = (fftw_r2r_kind*) malloc(2*sizeof(fftw_r2r_kind));
    wisdom.kindsr[0] = FFTW_REDFT10;
    wisdom.kindsr[1] = FFTW_REDFT10;
    wisdom.plans_dct3 = (fftw_plan*) fftw_malloc(sizeof(fftw_plan)*(wisdom.t-1));
    wisdom.plans_dct2 = (fftw_plan*) fftw_malloc(sizeof(fftw_plan)*(wisdom.t-1));
    wisdom.lengths = (int*) malloc((wisdom.t-1)*sizeof(int));
    for (i = 0, ti = 4; i < wisdom.t-1; i++, ti<<=1)
    {
      wisdom.lengths[i] = ti;
      wisdom.plans_dct3[i] = fftw_plan_many_r2r(1, &wisdom.lengths[i], 2, (double*)wisdom.vec1, NULL, 2, 1,
                                             (double*)wisdom.vec1, NULL, 2, 1, wisdom.kinds, 0);
      wisdom.plans_dct2[i] = fftw_plan_many_r2r(1, &wisdom.lengths[i], 2, (double*)wisdom.vec1, NULL, 2, 1,
                                             (double*)wisdom.vec1, NULL, 2, 1, wisdom.kindsr, 0);
    }  
    wisdom.U = precomputeU(wisdom.t, wisdom.threshold, wisdom.alpha, wisdom.beta, wisdom.gamma);
  }
  
  wisdom.initialized = true;
}

nfsft_plan nfsft_init(int d, int m, double *angles, fftw_complex **f_hat, 
                      fftw_complex *f, nfsft_flags flags)
{
  nfsft_plan plan = malloc(sizeof(struct nfsft_plan_s));
  
  plan->D = d;
  plan->M = m;
  plan->N = 1<<ngpt(plan->M);
  plan->angles = angles;
  plan->f_hat = f_hat;
  plan->f = f;
  
  /* Initialize NFFT. */
  nfft_size[0] = 2*(plan->N+1);
  nfft_size[1] = 2*(plan->N+1);
  fftw_size[0] = 4*plan->N;
  fftw_size[1] = 4*plan->N;
  nfft_init_specific(&plan->plan_nfft, 2, nfft_size, plan->D, fftw_size, 
                     6, PRE_PHI_HUT| PRE_PSI| MALLOC_F_HAT | FFT_OUT_OF_PLACE, 
                     FFTW_ESTIMATE| FFTW_DESTROY_INPUT);
  /* Assign angle array. */
  plan->plan_nfft.x = plan->angles;
  if (plan->plan_nfft.nfft_flags & PRE_PSI) 
  {  
    nfft_precompute_psi(&plan->plan_nfft); 
  } 
  /* Assign result array. */
  plan->plan_nfft.f = plan->f;
  
  return(plan);
}

void nfsft_finalize(nfsft_plan plan)
{
  nfft_finalize(&plan->plan_nfft);  
  free(plan);
}

void ndsft_trafo(nfsft_plan plan)
{
  /** Counter for loops */
  int i;
  
  precompute_coeffs();

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
  precompute_coeffs();
  
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
  //struct nfsft_transform_wisdom *tw;

  /* Init global structure. */
  precompute_coeffs();
  
  if (wisdom.initialized == false || plan->M > wisdom.N)
  {
    fprintf(stderr,"M too big!\n");
  }  

  if (plan->M < 3)
  {
    ndsft_trafo(plan);
  }
  else
  {
    t = ngpt(plan->M);
    /** Next greater power of 2 relative to M */
    N = 1<<t;
    
    /* Compute FLFT. */
    for (n = -plan->M, i = 0; n <= plan->M; n++, i++) 
    {
      flft(plan->M, t, abs(n), plan->f_hat[i], &wisdom);
    }
   
    /* Convert Chebyshev coefficients to Fourier coefficients. */
    cheb2exp(plan->plan_nfft.f_hat, plan->f_hat, plan->M, N); 
    
    /* Execute NFFT. */
    nfft_trafo(&plan->plan_nfft);
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
  //struct nfsft_transform_wisdom *tw;
  /** */
  //nfft_plan myplan;
  
  if (wisdom.initialized == false || plan->M > wisdom.N)
  {
    fprintf(stderr,"M too big!\n");
  }  
  
  /* Init global structure. */
  precompute_coeffs();
  
  if (plan->M < 3)
  {
    ndsft_adjoint(plan);
  }
  else
  {
    t = ngpt(plan->M);
    /** Next greater power of 2 relative to M */
    N = 1<<t;
        
    /* Execute adjoint NFFT. */
    nfft_adjoint(&plan->plan_nfft);
    
    /* Convert Chebyshev coefficients to Fourier coefficients. */
    cheb2exp_adjoint(plan->plan_nfft.f_hat, plan->f_hat, plan->M, N); 
        
    /* Calculate FLFT. */
    for(n = -plan->M, i = 0; n <= plan->M; n++, i++) 
    {
      flft_adjoint(plan->M, t, abs(n), plan->f_hat[i], &wisdom);
    }          
  }    
}
