#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "util.h"
#include "nfft3.h"

void accuracy_nsfft(int d, int J, int M, int m)
{
  nsfft_plan p;
  double complex *swap_sndft_trafo, *swap_sndft_adjoint;

  nsfft_init(&p, d, J, M, m, SNDFT);

  swap_sndft_trafo=(double complex*) fftw_malloc(p.M_total*
						 sizeof(double complex));
  swap_sndft_adjoint=(double complex*) fftw_malloc(p.N_total*
						   sizeof(double complex));

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
  int r, N[d], n[d];
  int m, m_nfft, m_nsfft;
  double t, t_nsdft, t_nfft, t_nsfft;

  nsfft_plan p;
  nfft_plan np;

  for(r=0;r<d;r++)
  {
    N[r]=int_2_pow(J+2);
    n[r]=(3*N[r])/2;
    /*n[r]=2*N[r];*/
  }

  /** init */
  m=total_used_memory();
  nsfft_init(&p, d, J, M, 4, SNDFT);
  m_nsfft=total_used_memory()-m;
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
    m=total_used_memory();
    nfft_init_guru(&np,d,N,M,n,6, FG_PSI| MALLOC_F_HAT| MALLOC_F| FFTW_INIT,
		   FFTW_MEASURE);
    m_nfft=total_used_memory()-m;
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
  {
    t_nfft=nan(""); 
    m_nfft=-1;
  }

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

  printf("%d\t%.2e\t%.2e\t%.2e\t%d\t%d\n",
	 J,
         t_nsdft,
	 t_nfft,
	 t_nsfft,
         m_nfft,
	 m_nsfft);

  fflush(stdout);

  /** finalise */
  nsfft_finalize(&p);
}


int main(int argc,char **argv)
{
  int d, J, M;
    
  if(argc<=2)
  {
    fprintf(stderr,"nsfft_test type d [first last trials]\n");
    return -1;
  }

  d=atoi(argv[2]);
  fprintf(stderr,"Testing the nfft on the hyperbolic cross (nsfft).\n");

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
