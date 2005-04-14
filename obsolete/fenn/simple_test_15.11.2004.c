#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/resource.h>
#include <time.h>
#include <string.h>

#include "nfft.h"
#include "utils.h"
#include "short_dim_2d.h"
#include "sparse_nfft.h"


void test_short()
{
  nfft_plan my_plan_2d, my_plan_1d;

  int N[2],n[2];
  int k,j,M;

  N[0]=1; N[1]=12; 
  n[0]=2*N[0]; n[1]=2*N[1];
  M=10;

  nfft_init_specific(&my_plan_2d,2,N,M,n,5,MALLOC_X| MALLOC_F| MALLOC_F_HAT, FFTW_MEASURE);
  nfft_init_specific(&my_plan_1d,1,&N[1],M,&n[1],5,MALLOC_X| MALLOC_F| MALLOC_F_HAT, FFTW_MEASURE);

  /* srand(1); */
  /* init frequencies */
  /*
  for(k=0;k<my_plan_2d.N_L;k++)
    {
      my_plan_2d.f_hat[k][0]=(double)rand()/RAND_MAX;
      my_plan_2d.f_hat[k][1]=(double)rand()/RAND_MAX;
    }
  */
  
  /* init nodes */
  for(j=0;j<my_plan_2d.M;j++) 
    {
      my_plan_2d.x[2*j+0]=((double)rand())/RAND_MAX-0.5;
      my_plan_2d.x[2*j+1]=((double)rand())/RAND_MAX-0.5;
      
      my_plan_2d.f[j][0]=((double)rand())/RAND_MAX-0.5;
      my_plan_2d.f[j][1]=((double)rand())/RAND_MAX-0.5;
    }
  
  /*  ndft_trafo(&my_plan_2d);
      vpr_c(my_plan_2d.f,my_plan_2d.M,"ndft_2d");
      
      short_nfft_trafo(&my_plan_2d,&my_plan_1d);
      vpr_c(my_plan_2d.f,my_plan_1d.M,"nfft_1d x ndft_1d");*/
  
  ndft_adjoint(&my_plan_2d);
  vpr_c(my_plan_2d.f_hat,10,"f_hat, ndft_2d");
      
  short_nfft_adjoint(&my_plan_2d,&my_plan_1d);
  vpr_c(my_plan_2d.f_hat,10,"f_hat, ndft_1d x nfft_1d");

  nfft_finalize(&my_plan_1d);
  nfft_finalize(&my_plan_2d);
}



void test_error2(int J, int m)
{
  int M, my_N;

  sparse_nfft_plan my_sparse_plan;
  fftw_complex* swap_sndft;
  
  my_N=int_2_pow(J+2);
  M=(J+4)*int_2_pow(J+1);

  /* SNFFT vs. SNDFT */

  /* initialise plans */
  sparse_nfft_init(&my_sparse_plan, J, M, m);

  /* memory allocation for swap */
  swap_sndft=(fftw_complex*) fftw_malloc(M*sizeof(fftw_complex));

  /* init nodes and frequencies */
  sparse_init_random_nodes_coeffs(&my_sparse_plan);

  sparse_ndft_trafo(&my_sparse_plan);
  SWAPC(swap_sndft,my_sparse_plan.f);

  sparse_nfft_trafo(&my_sparse_plan);

  printf("%e\t%e\t",E_2_error_c(swap_sndft,my_sparse_plan.f,M),
	 E_infty_error_c(swap_sndft,my_sparse_plan.f,M,my_sparse_plan.f_hat,my_sparse_plan.N_S));
  fflush(stdout);

  fftw_free(swap_sndft);

  sparse_nfft_finalize(&my_sparse_plan);

  /*---------------*/
  printf("\n\n\n");

  /* SNFFT^H vs. SNDFT^H */
  /* initialise plans */
  sparse_nfft_init(&my_sparse_plan, J, M, m);

  /* memory allocation for swap */
  swap_sndft=(fftw_complex*) fftw_malloc(my_sparse_plan.N_S*sizeof(fftw_complex));

  /* init nodes and frequencies */
  sparse_init_random_nodes_values(&my_sparse_plan);

  sparse_ndft_adjoint(&my_sparse_plan);
  SWAPC(swap_sndft,my_sparse_plan.f);

  sparse_nfft_adjoint(&my_sparse_plan);

  vpr_c(swap_sndft,my_sparse_plan.N_S,"snfft^h");
  vpr_c(my_sparse_plan.f_hat,my_sparse_plan.N_S,"snfft^h");
  printf("%e\n",E_2_error_c(swap_sndft,my_sparse_plan.f_hat,my_sparse_plan.N_S));
  fflush(stdout);

  fftw_free(swap_sndft);

  sparse_nfft_finalize(&my_sparse_plan);
}

void test_error(int J, int m)
{
  int M, my_N, r;

  sparse_nfft_plan my_sparse_plan;
  fftw_complex* swap_sndft;
  
  my_N=int_2_pow(J+2);
  M=(J+4)*int_2_pow(J+1);

  /* SNFFT^H vs. SNDFT^H */
  /* initialise plans */
  sparse_nfft_init(&my_sparse_plan, J, M, m);

  /* memory allocation for swap */
  swap_sndft=(fftw_complex*) fftw_malloc(my_sparse_plan.N_S*sizeof(fftw_complex));

  /* init nodes and frequencies */
  sparse_init_random_nodes_values(&my_sparse_plan);

  sparse_ndft_adjoint(&my_sparse_plan);
  SWAPC(swap_sndft,my_sparse_plan.f_hat);

  sparse_nfft_adjoint(&my_sparse_plan);

  for(r=0; r<4*((J+1)/2+1); r++)
    {
      printf("\n%d\n",r);
      vpr_c(swap_sndft+int_2_pow(J)*r,int_2_pow(J),"sndft^h");
      vpr_c(my_sparse_plan.f_hat+int_2_pow(J)*r,int_2_pow(J),"snfft^h");

      getc(stdin);
    }

  vpr_c(swap_sndft+4*((J+1)/2+1)*int_2_pow(J),int_2_pow(2*(J/2)+2),"center");
  vpr_c(my_sparse_plan.f_hat+4*((J+1)/2+1)*int_2_pow(J),int_2_pow(2*(J/2)+2),"center");

  printf("%e\n",E_2_error_c(swap_sndft,my_sparse_plan.f_hat,my_sparse_plan.N_S));
  fflush(stdout);

  fftw_free(swap_sndft);

  sparse_nfft_finalize(&my_sparse_plan);
}

void test_time(int J, int FULL_NFFT, int SPARSE_NDFT)
{
  int M, my_N;
  int m=5;

  double t;
  int N[2];
  int n[2];

  nfft_plan my_full_plan;
  sparse_nfft_plan my_sparse_plan;
  
  my_N=int_2_pow(J+2);
  M=(J+4)*int_2_pow(J+1);
  N[0]=my_N; N[1]=my_N; n[0]=2*my_N; n[1]=2*my_N;

  /* initialise plans */
  sparse_nfft_init(&my_sparse_plan, J, M, m);
  if(FULL_NFFT)
    nfft_init_specific(&my_full_plan,2,N,M,n,m,MALLOC_X| MALLOC_F| MALLOC_F_HAT, FFTW_MEASURE);

  /* init nodes and frequencies */
  sparse_init_random_nodes_coeffs(&my_sparse_plan);
  if(FULL_NFFT)
    {
      copy_sparse_to_full(&my_sparse_plan, &my_full_plan);
      /*test_copy_sparse_to_full(&my_sparse_plan, &my_full_plan);*/

      t=second();
      nfft_trafo(&my_full_plan);
      t=second()-t;
      printf("%e\t",t);
      //vpr_c(my_full_plan.f,my_full_plan.M,"full, f");
      //vpr_c(my_full_plan.f_hat,10,"full, f_hat");
    }
  else
    printf("inf\t");

  /* init nodes and values */
  /*sparse_init_random_nodes_values(&my_sparse_plan);
    if(FULL_NFFT)
    copy_sparse_to_full(&my_sparse_plan, &my_full_plan);*/


  if(SPARSE_NDFT)
    {
      t=second();
      sparse_ndft_trafo(&my_sparse_plan);
      t=second()-t;
      printf("%e\t",t);
      //vpr_c(my_sparse_plan.f,my_sparse_plan.M,"sparse, f");
      //copy_sparse_to_full(&my_sparse_plan, &my_full_plan);
      //vpr_c(my_full_plan.f_hat,10,"sparse, f_hat");
    }
  else
    printf("inf\t");

  t=second();
  sparse_nfft_trafo(&my_sparse_plan);
  t=second()-t;
  printf("%e\n",t);

  //printf("Error (fast sparse-fast full  ): %e\n",E_2_error_c(my_sparse_plan.f,my_full_plan.f,M));

  if(FULL_NFFT)
    nfft_finalize(&my_full_plan); 
  sparse_nfft_finalize(&my_sparse_plan);
}

void test_memory(int J, int FULL_NFFT, int SPARSE_NDFT)
{
  int M, my_N;
  int m=5;

  int N[2];
  int n[2];
  int mem;

  nfft_plan my_full_plan;
  sparse_nfft_plan my_sparse_plan;
  
  my_N=int_2_pow(J+2);
  M=(J+4)*int_2_pow(J+1);
  N[0]=my_N; N[1]=my_N; n[0]=2*my_N; n[1]=2*my_N;

  /* initialise my_full_plan */
  mem=total_used_memory();
  sparse_nfft_init(&my_sparse_plan, J, M, m);
  printf("%d\t",total_used_memory()-mem);

  if(FULL_NFFT)
    {
      mem=total_used_memory();
      nfft_init_specific(&my_full_plan,2,N,M,n,m,MALLOC_X| MALLOC_F| MALLOC_F_HAT, FFTW_MEASURE);
      printf("%d\n",total_used_memory()-mem);      
    }
  else
    printf("inf\n");

  if(FULL_NFFT)
    nfft_finalize(&my_full_plan); 

  sparse_nfft_finalize(&my_sparse_plan);
}

/**********************************************************/
/* main                                                   */
/**********************************************************/
int main(int argc, char *argv[])
{
  int J,JJ,m;

  /* evaluate parameters */
  if (argc<2)
    {
      printf("simple_test J\n");
      exit(-1);
    }

  J=atoi(argv[1]);  

  test_full_to_sparse(J);

  exit(-1);

  /* --------------------------------------------------- */

  test_error2(J,6);

  exit(-1);

  /* test_error */
  for(m=2;m<=10;m++)
    {
      printf("%d\t",m);
      test_error(J,m);
    }

  /* test_time */
  for(JJ=5;JJ<J;JJ++)
    {
      if(JJ<10)
	test_time(JJ,1,1);
      else
	{
	  if(JJ<13)
	    test_time(JJ,0,1);
	  else
	    test_time(JJ,0,0);
	}
      fflush(stdout);
    }

  /* test_memory */
  for(JJ=5;JJ<J;JJ++)
    {
      if(JJ<10)
	test_memory(JJ,1,1);
      else
	{
	  if(JJ<13)
	    test_memory(JJ,0,1);
	  else
	    test_memory(JJ,0,0);
	}
      fflush(stdout);
    }

  return(1);
}
