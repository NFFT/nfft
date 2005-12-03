/*! \file taylor_nfft.c
 *
 * \brief Testing the nfft againt a Taylor expansion based version.
 *
 * \author Stefan Kunis
 *
 * References: Time and memory requirements of the Nonequispaced FFT
 *
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "util.h"
#include "nfft3.h"

typedef struct
{
  nfft_plan p;                          /**< used for fftw and data          */

  int *idx0;                            /**< index of next neighbour of x_j
                                             on the oversampled regular grid */
  double *deltax0;                      /**< distance to the grid point      */
} taylor_plan;

/**
 * Initialisation of a transform plan.
 *
 * \arg ths The pointer to a taylor plan
 * \arg N The multi bandwidth
 * \arg M The number of nodes
 * \arg n The fft length
 * \arg m The order of the Taylor expansion
 *
 * \author Stefan Kunis
 */
void taylor_init(taylor_plan *ths, int N, int M, int n, int m)
{
  /* Note: no nfft precomputation! */
  nfft_init_guru((nfft_plan*)ths, 1, &N, M, &n, m,
                 MALLOC_X| MALLOC_F_HAT| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_ESTIMATE| FFTW_PRESERVE_INPUT);

  ths->idx0=(int*)fftw_malloc(M*sizeof(int));
  ths->deltax0=(double*)fftw_malloc(M*sizeof(double));
}

/**
 * Precomputation of weights and indices in Taylor expansion.
 *
 * \arg ths The pointer to a taylor plan
 *
 * \author Stefan Kunis
 */
void taylor_precompute(taylor_plan *ths)
{
  int j;

  nfft_plan* cths=(nfft_plan*)ths;

  for(j=0;j<cths->M_total;j++)
    {
      ths->idx0[j] = ((int)round((cths->x[j]+0.5)*cths->n[0]) +
                                 cths->n[0]/2)%cths->n[0];
      ths->deltax0[j] = cths->x[j] - (round((cths->x[j]+0.5)*cths->n[0]) /
                                      cths->n[0] - 0.5);
    }
}

/**
 * Destroys a transform plan.
 *
 * \arg ths The pointer to a taylor plan
 *
 * \author Stefan Kunis, Daniel Potts
 */
void taylor_finalize(taylor_plan *ths)
{
  fftw_free(ths->deltax0);
  fftw_free(ths->idx0);

  nfft_finalize((nfft_plan*)ths);
}

/**
 * Executes a Taylor-NFFT, see equation (1.1) in [Guide], computes fast and
 * approximate by means of a Taylor expansion
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 *
 * \arg ths The pointer to a taylor plan
 *
 * \author Stefan Kunis
 */
void taylor_trafo(taylor_plan *ths)
{
  int j,k,l,ll;
  complex *f, *f_hat, *g1;
  double *deltax;
  int *idx;

  nfft_plan *cths=(nfft_plan*)ths;

  for(j=0, f=cths->f; j<cths->M_total; j++)
    *f++ = 0;

  for(k=0; k<cths->n_total; k++)
    cths->g1[k]=0;

  for(k=-cths->N_total/2, g1=cths->g1+cths->n_total-cths->N_total/2,
      f_hat=cths->f_hat; k<0; k++)
    (*g1++)=cpow( - 2*PI*I*k,cths->m)* (*f_hat++);
   
  cths->g1[0]=cths->f_hat[cths->N_total/2];

  for(k=1, g1=cths->g1+1, f_hat=cths->f_hat+cths->N_total/2+1;
      k<cths->N_total/2; k++)
    (*g1++)=cpow( - 2*PI*I*k,cths->m)* (*f_hat++);

  for(l=cths->m-1; l>=0; l--)
    {
      for(k=-cths->N_total/2, g1=cths->g1+cths->n_total-cths->N_total/2;
          k<0; k++)
        (*g1++) /= (-2*PI*I*k);

      for(k=1, g1=cths->g1+1; k<cths->N_total/2; k++)
        (*g1++) /= (-2*PI*I*k);

      fftw_execute(cths->my_fftw_plan1);

      ll=(l==0?1:l);
      for(j=0, f=cths->f, deltax=ths->deltax0, idx=ths->idx0; j<cths->M_total;
          j++)
        (*f++) = ((*f) * (*deltax++) + cths->g2[*idx++]) /ll;
    }
}

/**
 * Compares NDFT, NFFT, and Taylor-NFFT
 *
 * \arg N The bandwidth
 * \arg N The number of nodes
 * \arg n The FFT-size for the NFFT
 * \arg m The cut-off for window function
 * \arg n_taylor The FFT-size for the Taylor-NFFT
 * \arg m_taylor The order of the Taylor approximation
 * \arg test_accuracy Flag for NDFT computation
 *
 * \author Stefan Kunis
 */
void taylor_time_accuracy(int N, int M, int n, int m, int n_taylor,
                          int m_taylor, unsigned test_accuracy)
{
  int j,k,r;

  double t_ndft, t_nfft, t_taylor, t;
  complex *swapndft;

  taylor_plan tp;
  nfft_plan np;

  printf("%d\t%d\t",N, M);

  taylor_init(&tp,N,M,n_taylor,m_taylor);

  nfft_init_guru(&np, 1, &N, M, &n, m,
                 PRE_PHI_HUT| PRE_FG_PSI|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_ESTIMATE| FFTW_DESTROY_INPUT);

  /** share nodes, input, and output vectors */
  np.x=tp.p.x;
  np.f_hat=tp.p.f_hat;
  np.f=tp.p.f;
  
  /** output vector ndft */
  if(test_accuracy)
    swapndft=(complex*)fftw_malloc(M*sizeof(complex));

  /** init pseudo random nodes */
  for(j=0;j<np.M_total;j++)
    np.x[j]=drand48()-0.5;

  /** nfft precomputation */
  taylor_precompute(&tp);

  /** nfft precomputation */
  if(np.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&np);

  /** init pseudo random Fourier coefficients */
  for(k=0;k<np.N_total;k++)
    np.f_hat[k] = drand48() + I* drand48();

  /** NDFT */
  if(test_accuracy)
    {
      SWAP_complex(np.f,swapndft);
      
      t_ndft=0;
      r=0; 
      while(t_ndft<0.01)
        {
          r++;
          t=second();
          ndft_trafo(&np);
          t=second()-t;
          t_ndft+=t;
        }
      t_ndft/=r;

      SWAP_complex(np.f,swapndft);
      printf("%.2e\t",t_ndft);
    }
  else 
    printf("nan\t\t");

  /** NFFT */
  t_nfft=0;
  r=0;
  while(t_nfft<0.01)
    {
      r++;
      t=second();
      nfft_trafo(&np);
      t=second()-t;
      t_nfft+=t;
    }
  t_nfft/=r;

  printf("%.2f\t%d\t%.2e\t",((double)n)/N, m, t_nfft);

  if(test_accuracy)
    printf("%.2e\t",error_l_infty_complex(swapndft, np.f, np.M_total));
  else
    printf("nan\t\t");

  /** TAYLOR NFFT */
  t_taylor=0;
  r=0;
  while(t_taylor<0.01)
    {
      r++;
      t=second();
      taylor_trafo(&tp);
      t=second()-t;
      t_taylor+=t;
    }
  t_taylor/=r;


  printf("%.2f\t%d\t%.2e\t",((double)n_taylor)/N,m_taylor,t_taylor);

  if(test_accuracy)
    printf("%.2e\n",error_l_infty_complex(swapndft, np.f, np.M_total));
  else
    printf("nan\t\n");

  fflush(stdout);
  
  /** finalise */
  if(test_accuracy)
    fftw_free(swapndft);

  nfft_finalize(&np);
  taylor_finalize(&tp);
}

int main(int argc,char **argv)
{
  int l,m,trial,sigma;
  int N=(1U<< 10);

  if(argc<=2)
    {
      fprintf(stderr,"taylor_nfft type first last trials sigma_nfft sigma_taylor.\n");
      return -1;
    }
  
  fprintf(stderr,"Testing the Nfft & a Taylor expansion based version.\n\n");
  fprintf(stderr,"Columns: N, M, t_ndft, sigma_nfft, m_nfft, t_nfft, e_nfft");
  fprintf(stderr,", sigma_taylor, m_taylor, t_taylor, e_taylor\n");

  /* time vs. N=M */
  if(atoi(argv[1])==0)
    {
      fprintf(stderr,"Fixed target accuracy, timings.\n\n");
      for(l=atoi(argv[2]); l<=atoi(argv[3]); l++)
        for(trial=0; trial<atoi(argv[4]); trial++)
          if(l<=10)
            taylor_time_accuracy((1U<< l), (1U<< l), (int)(atof(argv[5])*
                                 (1U<< l)), 6, (int)(atof(argv[6])*(1U<< l)),
                                 6, 1);
          else
            taylor_time_accuracy((1U<< l), (1U<< l), (int)(atof(argv[5])*
                                 (1U<< l)), 6, (int)(atof(argv[6])*(1U<< l)),
                                 6, 0);
    }

  /* error vs. m */
  if(atoi(argv[1])==1)
    {
      fprintf(stderr,"Fixed N=M=%d, error vs. m.\n\n",N);
      for(m=atoi(argv[2]); m<=atoi(argv[3]); m++)
        for(trial=0; trial<atoi(argv[4]); trial++)
          taylor_time_accuracy(N,N, (int)(atof(argv[5])*N), m,
                                    (int)(atof(argv[6])*N), m, 1);
    }

  return 1;
}
