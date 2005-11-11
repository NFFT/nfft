/*! \file flags.c
 *
 * \brief Testing the nfft.
 *
 * \author Stefan Kunis
 *
 * References: Time and Memory Requirements of the Nonequispaced FFT
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "util.h"
#include "nfft3.h"

void flags_cp(nfft_plan *dst, nfft_plan *src)
{
  dst->x=src->x;
  dst->f_hat=src->f_hat;
  dst->f=src->f;
  dst->g1=src->g1;
  dst->g2=src->g2;
  dst->my_fftw_plan1=src->my_fftw_plan1;
  dst->my_fftw_plan2=src->my_fftw_plan2;
}

void time_accuracy(int d, int N, int M, int n, int m, unsigned test_accuracy)
{
  int j,k,r;

  double t_ndft, t;
  complex *swapndft;

  int NN[d], nn[d];

  nfft_plan p;
  nfft_plan p_pre_phi_hut;
  nfft_plan p_fg_psi;
  nfft_plan p_pre_lin_psi;
  nfft_plan p_pre_fg_psi;
  nfft_plan p_pre_psi;
  nfft_plan p_pre_full_psi;

  for(r=0; r<d; r++)
    {
      NN[r]=N;
      nn[r]=n;
    }

  /** output vector ndft */
  if(test_accuracy)
    swapndft=(complex*)fftw_malloc(M*sizeof(complex));

  nfft_init_guru(&p, d, NN, M, nn, m,
                 MALLOC_X| MALLOC_F_HAT| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /** init pseudo random nodes */
  for(j=0; j<p.d*p.M_total; j++)
    p.x[j]=((double)rand())/RAND_MAX-0.5;

  nfft_init_guru(&p_pre_phi_hut, d, NN, M, nn, m, PRE_PHI_HUT,0);
  flags_cp(&p_pre_phi_hut, &p);
  if(p_pre_phi_hut.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&p_pre_phi_hut);

  nfft_init_guru(&p_fg_psi, d, NN, M, nn, m, FG_PSI,0);
  flags_cp(&p_fg_psi, &p);
  if(p_fg_psi.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&p_fg_psi);

  nfft_init_guru(&p_pre_lin_psi, d, NN, M, nn, m, PRE_LIN_PSI,0);
  flags_cp(&p_pre_lin_psi, &p);
  if(p_pre_lin_psi.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&p_pre_lin_psi);

  nfft_init_guru(&p_pre_fg_psi, d, NN, M, nn, m, PRE_FG_PSI,0);
  flags_cp(&p_pre_fg_psi, &p);
  if(p_pre_fg_psi.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&p_pre_fg_psi);

  nfft_init_guru(&p_pre_psi, d, NN, M, nn, m, PRE_PSI,0);
  flags_cp(&p_pre_psi, &p);
  if(p_pre_psi.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&p_pre_psi);

  nfft_init_guru(&p_pre_full_psi, d, NN, M, nn, m, PRE_FULL_PSI,0);
  flags_cp(&p_pre_full_psi, &p);
  if(p_pre_full_psi.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&p_pre_full_psi);

  /** init pseudo random Fourier coefficients */
  for(k=0; k<p.N_total; k++)
    p.f_hat[k] = ((double)rand())/RAND_MAX + I* ((double)rand())/RAND_MAX;

  /** NDFT */
  if(test_accuracy)
    {
      SWAP_complex(p.f,swapndft);
      
      t_ndft=0;
      r=0; 
      while(t_ndft<0.01)
        {
          r++;
          t=second();
          ndft_trafo(&p);
          t=second()-t;
          t_ndft+=t;
        }
      t_ndft/=r;

      SWAP_complex(p.f,swapndft);
    }
  else 
    {
      t_ndft=0;
    }

  /** NFFT */
  nfft_trafo(&p);
  nfft_trafo(&p_pre_phi_hut);
  nfft_trafo(&p_fg_psi);
  nfft_trafo(&p_pre_lin_psi);
  nfft_trafo(&p_pre_fg_psi);
  nfft_trafo(&p_pre_psi);
  nfft_trafo(&p_pre_full_psi);

  if(test_accuracy)
    printf("%.2e\t",error_l_2_complex(swapndft, p.f, p.M_total));
  else
    printf("--------\t");


  printf("%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
         t_ndft,
         p.MEASURE_TIME_t[0],
         p_pre_phi_hut.MEASURE_TIME_t[0],

         p.MEASURE_TIME_t[1],

         p.MEASURE_TIME_t[2],
         p_fg_psi.MEASURE_TIME_t[2],
         p_pre_lin_psi.MEASURE_TIME_t[2],
         p_pre_fg_psi.MEASURE_TIME_t[2],
         p_pre_psi.MEASURE_TIME_t[2],
         p_pre_full_psi.MEASURE_TIME_t[2]);

  /** finalise */
  nfft_finalize(&p_pre_full_psi);
  nfft_finalize(&p_pre_psi);
  nfft_finalize(&p_pre_fg_psi);
  nfft_finalize(&p_pre_lin_psi);
  nfft_finalize(&p_fg_psi);
  nfft_finalize(&p_pre_phi_hut);
  nfft_finalize(&p);

  if(test_accuracy)
    fftw_free(swapndft);
}



int main()
{
  int l,m;

  int d=3;

  for(l=3;l<20;l++)
    if(l<10)
      time_accuracy(d, (1U<< l), (1U<< (d*l)), (1U<< (l+1)), 5, 1);
    else
      time_accuracy(d, (1U<< l), (1U<< (d*l)), (1U<< (l+1)), 5, 0);


  return 1;
}
