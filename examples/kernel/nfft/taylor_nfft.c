/*! \file taylor_nfft.c
 *
 * \brief Testing the nfft againt a Taylor expansion based version.
 *
 * \author Stefan Kunis
 *
 * XXX References Anderson & Dahleh, Candes...
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

  double *deltax0;
  int *idx0;
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
  nfft_init_guru((nfft_plan*)ths, 1, &N, M, &n, m,
		 FFTW_INIT| FFT_OUT_OF_PLACE| MALLOC_X| MALLOC_F_HAT| MALLOC_F,
		 FFTW_MEASURE| FFTW_PRESERVE_INPUT);
  /* Note: no nfft precomputation! */

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
      ths->idx0[j] = ((int)round((cths->x[j]+0.5)*cths->n[0]) + cths->n[0]/2)%cths->n[0];
      ths->deltax0[j] = cths->x[j] - (round((cths->x[j]+0.5)*cths->n[0]) / cths->n[0] - 0.5);
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
 * Executes a NFFT, see equation (1.1) in [Guide], computes fast and
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

  for(k=0; k<cths->n[0]; k++)
    cths->g1[k]=0;

  for(k=-cths->N_total/2, g1=cths->g1+cths->n_total-cths->N_total/2, f_hat=cths->f_hat; k<0; k++)
    (*g1++)=cpow( - 2*PI*I*k,cths->m)* (*f_hat++);
   
  cths->g1[0]=cths->f_hat[cths->N_total/2];

  for(k=1, g1=cths->g1+1, f_hat=cths->f_hat+cths->N_total/2+1; k<cths->N_total/2; k++)
    (*g1++)=cpow( - 2*PI*I*k,cths->m)* (*f_hat++);

  for(l=cths->m-1; l>=0; l--)
    {
      for(k=-cths->N_total/2, g1=cths->g1+cths->n_total-cths->N_total/2; k<0; k++)
        (*g1++) /= (-2*PI*I*k);

      for(k=1, g1=cths->g1+1; k<cths->N_total/2; k++)
        (*g1++) /= (-2*PI*I*k);

      fftw_execute(cths->my_fftw_plan1);

      ll=(l==0?1:l);
      for(j=0, f=cths->f, deltax=ths->deltax0, idx=ths->idx0; j<cths->M_total; j++)
        (*f++) = ((*f) * (*deltax++) + cths->g2[*idx++]) /ll;
    }
}


void taylor_time_accuracy(int N, int M, int n, int m, int n_taylor,
                          int m_taylor, unsigned test_accuracy)
{
  int j,k,l,ll,r;

  double t_ndft, t_nfft, t_taylor, t;
  complex *swapndft;

  taylor_plan tp;
  nfft_plan np;

  taylor_init(&tp,N,M,n_taylor,m_taylor);

  nfft_init_guru(&np, 1, &N, M, &n, m,
                 PRE_PHI_HUT| //PRE_LIN_PSI|
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
  if(np.nfft_flags & PRE_LIN_PSI)
      nfft_precompute_lin_psi(&np);
  if(np.nfft_flags & PRE_PSI)
      nfft_precompute_psi(&np);
  if(np.nfft_flags & PRE_FULL_PSI)
      nfft_precompute_full_psi(&np);

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
    }
  else 
    {
      t_ndft=sqrt(-1);
    }

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

  printf("%d\t%d\t%.1f\t%d\t%.1f\t%d\t",N, M, ((double)n)/N, m, ((double)n_taylor)/N, m_taylor);

  if(test_accuracy)
    printf("%.2e\t",error_l_infty_1_complex(swapndft, np.f, np.M_total,
                    np.f_hat, np.N_total));
  else
    printf("--------\t");

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

  if(test_accuracy)
    printf("%.2e\t",error_l_infty_1_complex(swapndft, np.f, np.M_total,
                    np.f_hat, np.N_total));
  else
    printf("--------\t");
  

  printf("%.2e\t%.2e\t%.2e\n",t_ndft,t_nfft,t_taylor);

  /** finalise */
  if(test_accuracy)
    fftw_free(swapndft);

  nfft_finalize(&np);
  taylor_finalize(&tp);
}

void taylor_simple_test()
{
  int j,k,l;                            /**< index for nodes and freqencies  */
  nfft_plan my_plan, tp;       /**< plan for the nfft               */

  int N,M,m,n_taylor,facl;
  int* idx0;
  double* deltax0;
  int sigma;

  N=4096;
  M=4096;

  /** Taylor nfft */
  sigma=8;
  m=4;

  n_taylor=sigma*N;

  idx0=(int*)fftw_malloc(M*sizeof(int));
  deltax0=(double*)fftw_malloc(M*sizeof(double));

  /** init an one dimensional plan */
  nfft_init_1d(&my_plan, N, M);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
    {
      my_plan.x[j]=((double)rand())/RAND_MAX-0.5;
      idx0[j]=((int)round((my_plan.x[j]+0.5)*n_taylor) + n_taylor/2)%n_taylor;
      deltax0[j]=my_plan.x[j] - (round((my_plan.x[j]+0.5)*n_taylor) / n_taylor - 0.5);
    }

  /** precompute psi, the entries of the matrix B */
    if(my_plan.nfft_flags & PRE_LIN_PSI)
      nfft_precompute_lin_psi(&my_plan);

  if(my_plan.nfft_flags & PRE_PSI)
      nfft_precompute_psi(&my_plan);
      
  if(my_plan.nfft_flags & PRE_FULL_PSI)
      nfft_precompute_full_psi(&my_plan);

  /** init pseudo random Fourier coefficients and show them */
  for(k=0;k<my_plan.N_total;k++)
    my_plan.f_hat[k] = ((double)rand())/RAND_MAX + I* ((double)rand())/RAND_MAX;

  vpr_complex(my_plan.f_hat,10,"given Fourier coefficients, vector f_hat"); 

  /** direct trafo and show the result */
  ndft_trafo(&my_plan);
  vpr_complex(my_plan.f,10,"ndft, vector f"); 

  /** approx. trafo and show the result */
  nfft_trafo(&my_plan);
  vpr_complex(my_plan.f,10,"nfft, vector f");  

  /** Taylor nfft */

  /** init an one dimensional plan */
  nfft_init_guru(&tp, 1, &N, M, &n_taylor, 0,
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_ESTIMATE| FFTW_PRESERVE_INPUT);

  tp.x=my_plan.x;
  tp.f_hat=my_plan.f_hat;
  tp.f=my_plan.f;

  for(j=0; j<M; j++)
    tp.f[j]=0;

  for(k=0; k<n_taylor; k++)
    tp.g1[k]=0;

  for(l=0; l<m; l++)
    {
      facl=((l==0)?1:l*facl);

      for(k=-N/2; k<0; k++)
        tp.g1[n_taylor+k]=cpow( - 2*PI*I*k,l)*tp.f_hat[k+N/2];
      
      tp.g1[0]=tp.f_hat[N/2];

      for(k=1; k<N/2; k++)
        tp.g1[k]=cpow( - 2*PI*I*k,l)*tp.f_hat[k+N/2];

      fftw_execute(tp.my_fftw_plan1);

      for(j=0; j<M; j++)
        tp.f[j] += tp.g2[idx0[j]]*pow(deltax0[j],l)/facl;

      vpr_complex(tp.f,10,"taylor, vector f");  
    } 

  /** finalise the one dimensional plan */
  nfft_finalize(&tp);
  nfft_finalize(&my_plan);
  fftw_free(deltax0);
  fftw_free(idx0);
}

int main()
{
  int l,m;

//  time_accuracy_taylor_nfft_1d((1U<< 4), 100, (1U<< 5), 3, (1U<< 7), 4, 1); exit(-1);

printf("polynomial degree N,\nnumber of nodes M,\nnfft oversampling factor sigma,\nnfft truncation parameter m,\n");
printf("taylor nfft oversampling factor D,\ntaylor nfft truncation parameter L,\nerrors e=|F-F_approx|_infty/|f|_1, and\ntimes t in sec.\n\n"); 

printf("N\tM\tsigma\tm\tD\tL\te_nfft\t\te_taylornfft\tt_ndft\t\tt_nfft\t\tt_taylornfft\n");

//  for(m=1;m<14;m++)
//    time_accuracy_taylor_nfft_1d((1U<< 12), 10000, (1U<< 13), m, (1U<< 16), m, 1);

  for(l=4;l<20;l++)
    if(l<13)
      taylor_time_accuracy((1U<< l), (1U<< l), (1U<< (l+1)), 3, (1U<< (l+3)), 4, 1);
    else
      taylor_time_accuracy((1U<< l), (1U<< l), (1U<< (l+1)), 3, (1U<< (l+3)), 4, 0);

  return 1;
}
