#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "util.h"
#include "nfft3.h"

void simple_test_nsfft(int d, int J, int M)
{
  int j,k;                              /**< index for nodes and freqencies  */
  nsfft_plan p;                   /**< plan for the nfft               */
  complex *swap_sndft;

  nsfft_init(&p, d, J, M, 6, SNDFT);

  swap_sndft=(complex*) fftw_malloc(M*sizeof(complex));

  nsfft_init_random_nodes_coeffs(&p);

  vpr_complex(p.f_hat, 8,"frequencies, vector f_hat, 0,...,7");

  /** direct trafo and show the result */
  nsdft_trafo(&p);
  vpr_complex(p.f,p.M_total,"nsdft, vector f"); 

  /** approx. trafo and show the result */
  nsfft_trafo(&p);
  vpr_complex(p.f,p.M_total,"nsfft, vector f");

  /** direct adjoint and show the result */
  nsdft_adjoint(&p);
  vpr_complex(p.f_hat, 8,"adjoint nsdft, vector f_hat, 0,...,7");

  /** approx. adjoint and show the result */
  nsfft_adjoint(&p);
  vpr_complex(p.f_hat, 8,"adjoint nsfft, vector f_hat, 0,...,7");

  /** finalise the one dimensional plan */
  nsfft_finalize(&p);
}

void accuracy_nsfft(int d, int J, int M, int m)
{
  int j,k;                              /**< index for nodes and freqencies  */
  nsfft_plan p;                         /**< plan for the nfft               */
  complex *swap_sndft_trafo, *swap_sndft_adjoint;

  nsfft_init(&p, d, J, M, m, SNDFT);

  swap_sndft_trafo=(complex*) fftw_malloc(p.M_total*sizeof(complex));
  swap_sndft_adjoint=(complex*) fftw_malloc(p.N_total*sizeof(complex));

  nsfft_init_random_nodes_coeffs(&p);

  /** direct trafo */
  nsdft_trafo(&p);
  
  SWAP_complex(swap_sndft_trafo,p.f);

  /** approx. trafo */
  nsfft_trafo(&p);
  
  printf("%5d\t %+.5E\t",J, 
         error_l_infty_1_complex(swap_sndft_trafo, p.f, p.M_total,
                                 p.f_hat, p.N_total));
  fflush(stdout);

  vrand_unit_complex(p.f, p.M_total);

  /** direct adjoint */
  nsdft_adjoint(&p);
  
  SWAP_complex(swap_sndft_adjoint,p.f_hat);

  /** approx. adjoint */
  nsfft_adjoint(&p);
  
  printf("%+.5E\n", 
         error_l_infty_1_complex(swap_sndft_adjoint, p.f_hat,
                                 p.N_total,
                                 p.f, p.M_total));
  fflush(stdout);

  fftw_free(swap_sndft_adjoint);
  fftw_free(swap_sndft_trafo);

  /** finalise the one dimensional plan */
  nsfft_finalize(&p);
}

void time_nsfft(int d, int J, int M, unsigned test_nsdft, unsigned test_nfft)
{
  int j,k,r;                            /**< index for nodes and freqencies  */
  nsfft_plan p;                         /**< plan for the nsfft              */
  nfft_plan np;                         /**< plan for the nfft               */

  double t,t_nsdft,t_nfft,t_nsfft;

  int N[2]={int_2_pow(J+2),int_2_pow(J+2)};
  int n[2]={2*N[0],2*N[1]};

  /** init */
  nsfft_init(&p, d, J, M, 4, SNDFT);
  nsfft_init_random_nodes_coeffs(&p);

  /* transforms */
  if(test_nsdft)
  {
    t_nsdft=0;
    r=0;
    while(t_nsdft<0.1)
    {
      r++;
      t=second();
      nsdft_trafo(&p);
      t=second()-t;
      t_nsdft+=t;
    }
    t_nsdft/=r;
  }
  else
    t_nsdft=nan("");   

  if(test_nfft)
  {
    nfft_init_guru(&np,d,N,M,n,4, FG_PSI| MALLOC_F_HAT| MALLOC_F| FFTW_INIT, FFTW_MEASURE);
    np.x=p.act_nfft_plan->x;
    if(np.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&np);
    nsfft_cp(&p, &np);

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

    nfft_finalize(&np);
  }
  else
    t_nfft=nan(""); 

  t_nsfft=0;
  r=0;
  while(t_nsfft<0.1)
    {
      r++;
      t=second();
      nsfft_trafo(&p);
      t=second()-t;
      t_nsfft+=t;
    }
  t_nsfft/=r;

  printf("%d\t%.2e\t%.2e\t%.2e\n",
	 J,
         t_nsdft,
	 t_nfft,
	 t_nsfft);

  fflush(stdout);

  /** finalise */
  nsfft_finalize(&p);
}


int main(int argc,char **argv)
{
  int l,m,d,trial,N,K,J,M;
    
  if(argc<=2)
  {
    fprintf(stderr,"simple_test type d [first last trials]\n");
    return -1;
  }

  d=atoi(argv[2]);
  fprintf(stderr,"Testing the nsfft (nfft on the hyperbolic cross).\n");

  if(atoi(argv[1])==0)
  {
    fprintf(stderr,"Computing a %d dimensional nsdft and nsfft\n\n",d);
    simple_test_nsfft(d,5,8);
  }

  if(atoi(argv[1])==1)
  {
    fprintf(stderr,"Testing the accuracy of the nsfft vs. nsdft\n");
    fprintf(stderr,"Columns: d, E_{1,\\infty}(trafo) E_{1,\\infty}(adjoint)\n\n");
    for(J=1; J<10; J++)
      accuracy_nsfft(d, J, 1000, 6);
  }
 
  if(atoi(argv[1])==2)
  {
    fprintf(stderr,"Testing the computation time of the nsdft, nfft, and nsfft\n");
    fprintf(stderr,"Columns: d, J, M, t_nsdft, t_nfft, t_nsfft\n\n");
    for(J=atoi(argv[3]); J<=atoi(argv[4]); J++)
    {
      if(d==2)
	M=(J+4)*int_2_pow(J+1);
      else
	M=6*int_2_pow(J)*(int_2_pow((J+1)/2+1)-1)+int_2_pow(3*(J/2+1));
      
      if(d*(J+2)<=24)
	time_nsfft(d, J, M, 1, 1);
      else
	if(d*(J+2)<=24)
	  time_nsfft(d, J, M, 0, 1);
	else
	  time_nsfft(d, J, M, 0, 0);  
    }
  }

  return 1;
}
