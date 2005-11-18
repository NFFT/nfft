/*! \file taylor_nfft.c
 *
 * \brief Testing the nfft againt a Taylor expansion based version.
 *
 * \author Stefan Kunis
 *
 * References: [AnDa96], i.e., Chris Anderson and Marie Dahleh:
 * Rapid computation on the discrete Fourier transform
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
		 FFTW_MEASURE| FFTW_PRESERVE_INPUT);

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

  double t_ndft, t_nfft, t_taylor, t, tc;
  complex *swapndft;

  taylor_plan tp;
  nfft_plan np;

  printf("%d\t%d\t",N, M);

  taylor_init(&tp,N,M,n_taylor,m_taylor);

  nfft_init_guru(&np, 1, &N, M, &n, m,
                 PRE_PHI_HUT| FG_PSI|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /** share nodes, input, and output vectors */
  np.x=tp.p.x;
  np.f_hat=tp.p.f_hat;
  np.f=tp.p.f;
  
  /** output vector ndft */
  if(test_accuracy)
    swapndft=(complex*)fftw_malloc(M*sizeof(complex));

  /** init pseudo random nodes */
  for(j=0;j<np.M_total;j++)
    np.x[j]=((double)rand())/RAND_MAX-0.5;

  /** nfft precomputation */
  taylor_precompute(&tp);

  /** nfft precomputation */
  if(np.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&np);

  /** init pseudo random Fourier coefficients */
  for(k=0;k<np.N_total;k++)
    np.f_hat[k] = ((double)rand())/RAND_MAX + I* ((double)rand())/RAND_MAX;

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

  printf("%.1f\t%d\t%.2e\t",((double)n)/N, m, t_nfft);

  if(test_accuracy)
    printf("%.2e\t",error_l_infty_complex(swapndft, np.f, np.M_total));
  else
    printf("nan\t\t");

  /** TAYLOR NFFT */
  t_taylor=0;
  r=0;
tc=second();
  while(t_taylor<0.01)
    {
      r++;
      t=second();
      taylor_trafo(&tp);
      t=second()-t;
      t_taylor+=t;
    }
tc=second()-tc;
tc/=r;
  t_taylor/=r;


  printf("%.1f\t%d\t%.2e\t%.2e\t",((double)n_taylor)/N,m_taylor,t_taylor,tc);

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

int main()
{
  int l,m;

  for(l=4;l<15;l++)
    taylor_time_accuracy((1U<< l), (1U<< l), (1U<< (l+2)), 5, (1U<< (l+3)), 5, 1);
  for(l=15;l<20;l++)
    taylor_time_accuracy((1U<< l), (1U<< l), (1U<< (l+2)), 5, (1U<< (l+3)), 5, 0);


printf("polynomial degree N,\nnumber of nodes M,\nnfft oversampling factor sigma,\nnfft truncation parameter m,\n");
printf("taylor nfft oversampling factor D,\ntaylor nfft truncation parameter L,\nerrors e=|F-F_approx|_infty/|f|_1, and\ntimes t in sec.\n\n"); 

printf("N\tM\tsigma\tm\tD\tL\te_nfft\t\te_taylornfft\tt_ndft\t\tt_nfft\t\tt_taylornfft\n");


/*  for(m=0;m<20;m++)
    taylor_time_accuracy((1U<< 12), (1U<< 12), (1U<< (14)), m, (1U<< (15)), m, 1);
*/
printf("\n");

  for(l=4;l<20;l++)
    if(l<13)
      taylor_time_accuracy((1U<< l), (1U<< l), (1U<< (l+2)), 4, (1U<< (l+3)), 5, 1);
    else
      taylor_time_accuracy((1U<< l), (1U<< l), (1U<< (l+2)), 4, (1U<< (l+3)), 5, 0);

  return 1;
}
