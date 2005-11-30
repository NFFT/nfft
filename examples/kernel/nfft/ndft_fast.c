/*! \file ndft_fast.c
 *
 * \brief Testing ndft, Horner-like ndft, and fully precomputed ndft.
 *
 * \author Stefan Kunis
 *
 * References: [BaMi], i.e., Bagchi Mitra XXXXXXXX:
 *
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "util.h"
#include "nfft3.h"

void ndft_horner_trafo(nfft_plan *ths)
{
  int j,k;
  complex *f_hat_k, *f_j;
  complex exp_omega_0, exp_omega;

  for(j=0, f_j=ths->f; j<ths->M_total; j++, f_j++)
    (*f_j) =0;

  for(j=0, f_j=ths->f; j<ths->M_total; j++, f_j++)
    {
      exp_omega_0 = cexp(+2*PI*I*ths->x[j]);
      for(k=0, f_hat_k= ths->f_hat; k<ths->N[0]; k++, f_hat_k++)
        {
          (*f_j)+=(*f_hat_k);
          (*f_j)*=exp_omega_0;
	}
      (*f_j)*=cexp(-PI*I*ths->N[0]*ths->x[j]);
    }
} /* ndft_horner_trafo */

void ndft_pre_full_trafo(nfft_plan *ths, complex *A)
{
  int j,k;
  complex *f_hat_k, *f_j;
  complex *A_local;

  for(j=0, f_j=ths->f; j<ths->M_total; j++, f_j++)
    (*f_j) =0;

  for(j=0, f_j=ths->f, A_local=A; j<ths->M_total; j++, f_j++)
    for(k=0, f_hat_k= ths->f_hat; k<ths->N[0]; k++, f_hat_k++, A_local++)
      (*f_j) += (*f_hat_k)*(*A_local);
} /* ndft_pre_full_trafo */

void ndft_pre_full_init(nfft_plan *ths, complex *A)
{
  int j,k;
  complex *f_hat_k, *f_j, *A_local;

  for(j=0, f_j=ths->f, A_local=A; j<ths->M_total; j++, f_j++)
    for(k=0, f_hat_k= ths->f_hat; k<ths->N[0]; k++, f_hat_k++, A_local++)
      (*A_local) = cexp(-2*PI*I*(k-ths->N[0]/2)*ths->x[j]);

} /* ndft_pre_full_init */

void ndft_time(int N, int M, unsigned test_ndft, unsigned test_pre_full)
{
  nfft_plan np;
  int j,k,r;
  double t, t_ndft, t_horner, t_pre_full, t_nfft;
  complex *A;
  int n=2*N;

  printf("%d\t%d\t",N, M);

  nfft_init_1d(&np, N, M);

  /** init pseudo random nodes */
  for(j=0;j<np.M_total;j++)
    np.x[j]=drand48()-0.5;

  if(test_pre_full)
   {
     A=(complex*)fftw_malloc(N*M*sizeof(complex));
     ndft_pre_full_init(&np, A);
   }

  /** init pseudo random Fourier coefficients */
  for(k=0;k<np.N_total;k++)
    np.f_hat[k] = drand48() + I* drand48();

  /** NDFT */
  if(test_ndft)
    {
      t_ndft=0;
      r=0;
      while(t_ndft<0.1)
        {
          r++;
          t=second();
          ndft_trafo(&np);
          t=second()-t;
          t_ndft+=t;
        }
      t_ndft/=r;

      printf("%.2e\t",t_ndft);
    }
  else
    printf("nan\t\t");

  /** Horner NDFT */
  t_horner=0;
  r=0;
  while(t_horner<0.1)
    {
      r++;
      t=second();
      ndft_horner_trafo(&np);
      t=second()-t;
      t_horner+=t;
    }
  t_horner/=r;

  printf("%.2e\t", t_horner);

  /** Fully precomputed NDFT */
  if(test_pre_full)
    {
      t_pre_full=0;
      r=0;
      while(t_pre_full<0.1)
        {
          r++;
          t=second();
          ndft_pre_full_trafo(&np,A);
          t=second()-t;
          t_pre_full+=t;
        }
      t_pre_full/=r;

      printf("%.2e\t", t_pre_full);
    }
  else
    printf("nan\t\t");

  t_nfft=0;
  r=0;
  while(t_nfft<0.1)
    {
      r++;
      t=second();
      nfft_trafo(&np);
      t=second()-t;
      t_nfft+=t;
    }
  t_nfft/=r;

  printf("%.2e\n", t_nfft);

  fflush(stdout);

  if(test_pre_full)
    fftw_free(A);

  nfft_finalize(&np);
}

int main(int argc,char **argv)
{
  int l,m,trial,sigma;
  int N=(1U<< 10);

  if(argc<=2)
    {
      fprintf(stderr,"ndft_fast type first last trials\n");
      return -1;
    }
  
  fprintf(stderr,"Testing ndft, Horner-like ndft, fully precomputed ndft.\n");
  fprintf(stderr,"Columns: N, M, t_ndft, t_horner, t_pre_full, t_nfft\n\n");

  /* time vs. N=M */
  if(atoi(argv[1])==0)
    {
      for(l=atoi(argv[2]); l<=atoi(argv[3]); l++)
        for(trial=0; trial<atoi(argv[4]); trial++)
          if(l<13)
            ndft_time((1U<< l), (1U<< l), 1, 1);
          else
            if(l<15)
              ndft_time((1U<< l), (1U<< l), 1, 0);
            else
              ndft_time((1U<< l), (1U<< l), 0, 0);
    }

  return 1;
}
