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


/** test functions during debugging */
void test_time_cos()
{
  const int n=100000000;
  int k;
  double t,y;
  double X[n];

  for(k=0;k<n;k++)
    X[k]=(double)rand()/RAND_MAX;

  t=second();
  y=0;
  for(k=0;k<n;k++)
    {
      y=y+cos(2*PI*k*X[k])+sin(2*PI*k*X[k]);
      //y=y+2*PI*k*X[k]+X[k]*X[k];
    }
  t=second()-t;
  printf("time=%e, result=%e\n",t,y);

}

void test_short()
{
  nfft_plan my_plan_2d, my_plan_1d;
  fftw_complex *exp_omega, *swap1, *swap2;
  double t1,t2,t3;

  int N[2],n[2];
  int k,j,M;

  N[0]=4; N[1]=128; 
  n[0]=2*N[0]; n[1]=2*N[1];
  M=N[0]*N[1];

  nfft_init_specific(&my_plan_2d,2,N,M,n,5,MALLOC_X| MALLOC_F| MALLOC_F_HAT| FFTW_INIT, FFTW_MEASURE);
  nfft_init_specific(&my_plan_1d,1,&N[1],M,&n[1],5,MALLOC_X| MALLOC_F| FFTW_INIT, FFTW_MEASURE);

  exp_omega = (fftw_complex*) fftw_malloc(my_plan_2d.M*sizeof(fftw_complex));
  swap1 = (fftw_complex*) fftw_malloc(my_plan_2d.M*sizeof(fftw_complex));
  swap2 = (fftw_complex*) fftw_malloc(my_plan_2d.M*sizeof(fftw_complex));

  srand(1);
  /* init frequencies */
  for(k=0;k<my_plan_2d.N_L;k++)
    {
      my_plan_2d.f_hat[k][0]=(double)rand()/RAND_MAX;
      my_plan_2d.f_hat[k][1]=(double)rand()/RAND_MAX;
    }
  
  /* init nodes */
  for(j=0;j<my_plan_2d.M;j++) 
    {
      my_plan_2d.x[2*j+0]=((double)rand())/RAND_MAX-0.5;
      my_plan_2d.x[2*j+1]=((double)rand())/RAND_MAX-0.5;
    }
  
  t1=second();
  ndft_trafo(&my_plan_2d);
  t1=second()-t1;
  //vpr_c(my_plan_2d.f,my_plan_2d.M,"ndft_2d");
  SWAPC(swap1,my_plan_2d.f);
  
  t2=second();
  short_nfft_trafo(&my_plan_2d,&my_plan_1d);
  t2=second()-t2;
  //vpr_c(my_plan_2d.f,my_plan_1d.M,"nfft_1d x ndft_1d");
  SWAPC(swap2,my_plan_2d.f);

  for(j=0;j<M;j++)
    {
      exp_omega[j][0] = cos(2*PI*my_plan_2d.x[2 * j + 0]);
      exp_omega[j][1] = sin(2*PI*my_plan_2d.x[2 * j + 0]);
    }

  t3=second();
  short_nfft_trafo_horner(&my_plan_2d,&my_plan_1d,exp_omega);
  t3=second()-t3;
  //vpr_c(my_plan_2d.f,my_plan_1d.M,"nfft_1d x ndft_1d");

  printf("%e\t%e\t\t%e\t%e\n",
	 E_2_error_c(swap2, swap1, M),
	 E_infty_error_c(swap2, swap1, M, my_plan_2d.f_hat, my_plan_2d.N_L),
	 E_2_error_c(my_plan_2d.f, swap1, M),
	 E_infty_error_c(my_plan_2d.f, swap1, M, my_plan_2d.f_hat, my_plan_2d.N_L));

  printf("%e\t\t%e\t\t%e\n",t1,t2,t3);

  fftw_free(swap2);
  fftw_free(swap1);
  fftw_free(exp_omega);
  nfft_finalize(&my_plan_1d);
  nfft_finalize(&my_plan_2d);
}

void test_snfft_3d(int J, int m)
{
  int M, my_N;
  int N[3];
  int n[3];
  double t;

  sparse_nfft_plan_3d my_sparse_plan;

  my_N=int_2_pow(J+2);
  N[0]=my_N; n[0]=2*my_N;
  N[1]=my_N; n[1]=2*my_N;
  N[2]=my_N; n[2]=2*my_N;
  
  //M=6*int_2_pow(J)*(int_2_pow((J+1)/2+1)-1)+int_2_pow(3*(J/2+1));
  M=10;

  /* NSFFT vs. SNDFT */

  /* initialise plans */

  t=second();
  sparse_nfft_init_3d(&my_sparse_plan, J, M, m, SNDFT);
  t=second()-t;
  printf("init done: %f sec\n",t);

  /* init nodes and frequencies */
  sparse_init_random_nodes_coeffs_3d(&my_sparse_plan);
  
  t=second();
  sparse_nfft_trafo_3d(&my_sparse_plan);
  t=second()-t;
  printf("snfft done: %f sec\n",t);

  vpr_c(my_sparse_plan.f,M,"snfft:f"); //fflush(stdout);

  t=second();
  sparse_ndft_trafo_3d(&my_sparse_plan);
  t=second()-t;
  printf("sndft done: %f sec\n",t);

  vpr_c(my_sparse_plan.f,M,"sndft:f");

  sparse_nfft_finalize_3d(&my_sparse_plan);
}


/** test error snfft vs. sndft
 *  d=2
 *  with respect to m (and J)
 */
void test_error(int J, int m)
{
  int M, my_N;

  sparse_nfft_plan my_sparse_plan;
  fftw_complex* swap_sndft;
  
  my_N=int_2_pow(J+2);
  M=(J+4)*int_2_pow(J+1);

  /* SNFFT vs. SNDFT */

  /* initialise plans */
  sparse_nfft_init(&my_sparse_plan, J, M, m, SNDFT);

  /* memory allocation for swap */
  swap_sndft=(fftw_complex*) fftw_malloc(M*sizeof(fftw_complex));

  /* init nodes and frequencies */
  sparse_init_random_nodes_coeffs(&my_sparse_plan);
  fflush(stdout);

  sparse_ndft_trafo(&my_sparse_plan);

  SWAPC(swap_sndft,my_sparse_plan.f);

  sparse_nfft_trafo(&my_sparse_plan);

  printf("%e\t%e\n",E_2_error_c(swap_sndft,my_sparse_plan.f,M),
	 E_infty_error_c(swap_sndft,my_sparse_plan.f,M,my_sparse_plan.f_hat,my_sparse_plan.N_S));
  fflush(stdout);

  fftw_free(swap_sndft);

  sparse_nfft_finalize(&my_sparse_plan);
}

/** test error snfft vs. sndft and nfft vs. sndft
 *  d=3
 *  with respect to m (and J)
 */
void test_error_3d(int J, int m, int FULL_NFFT)
{
  int M, my_N;
  int N[3];
  int n[3];

  sparse_nfft_plan_3d my_sparse_plan;
  fftw_complex* swap_sndft;

  nfft_plan my_full_plan;

  my_N=int_2_pow(J+2);
  N[0]=my_N; n[0]=2*my_N;
  N[1]=my_N; n[1]=2*my_N;
  N[2]=my_N; n[2]=2*my_N;
  
  M=6*int_2_pow(J)*(int_2_pow((J+1)/2+1)-1)+int_2_pow(3*(J/2+1));

  /* memory allocation for swap */
  swap_sndft=(fftw_complex*) fftw_malloc(M*sizeof(fftw_complex));

  /* initialise plans */
  sparse_nfft_init_3d(&my_sparse_plan, J, M, m, SNDFT);

  /* init nodes and frequencies */
  sparse_init_random_nodes_coeffs_3d(&my_sparse_plan);
  
  sparse_ndft_trafo_3d(&my_sparse_plan);
  
  SWAPC(swap_sndft,my_sparse_plan.f);

  sparse_nfft_trafo_3d(&my_sparse_plan);

  if(FULL_NFFT)
    {
      nfft_init_specific(&my_full_plan,3,N,M,n,m, MALLOC_X| MALLOC_F| MALLOC_F_HAT| FFTW_INIT, FFTW_MEASURE);
      copy_sparse_to_full_3d(&my_sparse_plan, &my_full_plan);
      nfft_trafo(&my_full_plan);
      
      printf("%e\t%e\t\t%e\t%e\n",
	     E_2_error_c(my_sparse_plan.f, swap_sndft, M),
	     E_infty_error_c(my_sparse_plan.f, swap_sndft, M, my_sparse_plan.f_hat, my_sparse_plan.N_S),
	     E_2_error_c(my_full_plan.f, swap_sndft, M),
	     E_infty_error_c(my_full_plan.f, swap_sndft, M, my_sparse_plan.f_hat, my_sparse_plan.N_S));

      nfft_finalize(&my_full_plan);
    }
  else
    {
      printf("%e\t%e\t\tnan\tnan\n",
	     E_2_error_c(my_sparse_plan.f, swap_sndft, M),
	     E_infty_error_c(my_sparse_plan.f, swap_sndft, M, my_sparse_plan.f_hat, my_sparse_plan.N_S));
    }
  fflush(stdout);

  sparse_nfft_finalize_3d(&my_sparse_plan);
  fftw_free(swap_sndft);
}

/** test time snfft vs. sndft vs. nfft
 *  d=2
 *  with respect to J (and m)
 */
void test_time(int J, int m, int FULL_NFFT, int SPARSE_NDFT)
{
  int M, my_N;

  double t;
  int N[2];
  int n[2];

  nfft_plan my_full_plan;
  sparse_nfft_plan my_sparse_plan;
  
  my_N=int_2_pow(J+2);
  N[0]=my_N; n[0]=2*my_N;
  N[1]=my_N; n[1]=2*my_N;

  M=(J+4)*int_2_pow(J+1);
 
  /* initialise plans */
  if(SPARSE_NDFT)
    sparse_nfft_init(&my_sparse_plan, J, M, m, SNDFT);
  else
    sparse_nfft_init(&my_sparse_plan, J, M, m, 0);
 
  /* init nodes and frequencies */
  sparse_init_random_nodes_coeffs(&my_sparse_plan);

  t=second();
  sparse_nfft_trafo(&my_sparse_plan);
  t=second()-t;
  printf("%e\t",t);
  fflush(stdout);

  if(SPARSE_NDFT)
    {
      t=second();
      sparse_ndft_trafo(&my_sparse_plan);
      t=second()-t;
      printf("%e\t",t);
    }
  else
    printf("nan\t");
  fflush(stdout);
  
  if(FULL_NFFT)
    {
      nfft_init_specific(&my_full_plan,2,N,M,n,m,MALLOC_X| MALLOC_F| MALLOC_F_HAT| FFTW_INIT, FFTW_MEASURE);
      copy_sparse_to_full(&my_sparse_plan, &my_full_plan);

      t=second();
      nfft_trafo(&my_full_plan);
      t=second()-t;
      printf("%e\n",t);

      nfft_finalize(&my_full_plan); 
    }
  else
    printf("nan\n");
  fflush(stdout);

  sparse_nfft_finalize(&my_sparse_plan);
}

/** test time snfft vs. sndft vs. nfft
 *  d=3
 *  with respect to J (and m)
 */
void test_time_3d(int J, int m, int FULL_NFFT, int SPARSE_NDFT)
{
  int M, my_N;
  double t;

  int N[3],n[3];

  nfft_plan my_full_plan;
  sparse_nfft_plan_3d my_sparse_plan;

  my_N=int_2_pow(J+2);
  N[0]=my_N; n[0]=2*my_N;
  N[1]=my_N; n[1]=2*my_N;
  N[2]=my_N; n[2]=2*my_N;
  
  M=6*int_2_pow(J)*(int_2_pow((J+1)/2+1)-1)+int_2_pow(3*(J/2+1));

  /* SNFFT, NFFT vs. SNDFT */

  /* initialise plans */
  if(SPARSE_NDFT)
    sparse_nfft_init_3d(&my_sparse_plan, J, M, m, SNDFT);
  else
    sparse_nfft_init_3d(&my_sparse_plan, J, M, m, 0);

  /* init nodes and frequencies */
  sparse_init_random_nodes_coeffs_3d(&my_sparse_plan);

  t=second();
  sparse_nfft_trafo_3d(&my_sparse_plan);
  t=second()-t;
  printf("%e\t",t);
  fflush(stdout);	
  
  if(SPARSE_NDFT)
    {
      t=second();
      sparse_ndft_trafo_3d(&my_sparse_plan);
      t=second()-t;
      printf("%e\t",t);
    }
  else
    printf("nan\t");
  fflush(stdout);

  if(FULL_NFFT)
    {
      nfft_init_specific(&my_full_plan,3,N,M,n,m, MALLOC_X| MALLOC_F| MALLOC_F_HAT| FFTW_INIT, FFTW_MEASURE);
      copy_sparse_to_full_3d(&my_sparse_plan, &my_full_plan);

      t=second();
      nfft_trafo(&my_full_plan);
      t=second()-t;
   
      printf("%e\n",t);
      
      nfft_finalize(&my_full_plan);
    }
  else
    printf("nan\n");
  fflush(stdout);

  sparse_nfft_finalize_3d(&my_sparse_plan);
}

/** test memory snfft vs. nfft
 *  d=2
 *  with respect to J (and m)
 */
void test_memory(int J, int m, int FULL_NFFT)
{
  int M, my_N;

  int N[2];
  int n[2];
  int mem;

  nfft_plan my_full_plan;
  sparse_nfft_plan my_sparse_plan;

  my_N=int_2_pow(J+2);
  N[0]=my_N; n[0]=2*my_N;
  N[1]=my_N; n[1]=2*my_N;

  M=(J+4)*int_2_pow(J+1);
  
  mem=total_used_memory();
  sparse_nfft_init(&my_sparse_plan, J, M, m, 0);
  mem=total_used_memory()-mem;
  /* nfft without ndft indices */
  printf("%d\t",mem-my_sparse_plan.N_S*sizeof(int));
  /* ndft without nfft precomputed psi */
  printf("%d\t",mem-(my_sparse_plan.act_nfft_plan->K+1)*my_sparse_plan.act_nfft_plan->d*sizeof(double));
  fflush(stdout);
  sparse_nfft_finalize(&my_sparse_plan);

  if(FULL_NFFT)
    {
      mem=total_used_memory();
      nfft_init_specific(&my_full_plan,2,N,M,n,m,MALLOC_X| MALLOC_F| MALLOC_F_HAT| FFTW_INIT, FFTW_MEASURE);
      printf("%d\n",total_used_memory()-mem);      
      nfft_finalize(&my_full_plan); 
    }
  else
    printf("nan\n");
  fflush(stdout);
}

/** test memory snfft vs. nfft
 *  d=3
 *  with respect to J (and m)
 */
void test_memory_3d(int J, int m, int FULL_NFFT)
{
  int M, my_N;

  int N[3];
  int n[3];
  int mem;

  nfft_plan my_full_plan;
  sparse_nfft_plan_3d my_sparse_plan;

  my_N=int_2_pow(J+2);
  N[0]=my_N; n[0]=2*my_N;
  N[1]=my_N; n[1]=2*my_N;
  N[2]=my_N; n[2]=2*my_N;

  M=6*int_2_pow(J)*(int_2_pow((J+1)/2+1)-1)+int_2_pow(3*(J/2+1));
  
  mem=total_used_memory();
  sparse_nfft_init_3d(&my_sparse_plan, J, M, m, 0);
  mem=total_used_memory()-mem;
  /* nfft without ndft indices */
  printf("%d\t",mem-my_sparse_plan.N_S*sizeof(int));
  /* ndft without nfft precomputed psi */
  printf("%d\t",mem-(my_sparse_plan.act_nfft_plan->K+1)*my_sparse_plan.act_nfft_plan->d*sizeof(double));
  fflush(stdout);
  sparse_nfft_finalize_3d(&my_sparse_plan);

  if(FULL_NFFT)
    {
      mem=total_used_memory();
      nfft_init_specific(&my_full_plan,3,N,M,n,m,MALLOC_X| MALLOC_F| MALLOC_F_HAT| FFTW_INIT, FFTW_MEASURE);
      printf("%d\n",total_used_memory()-mem);      
      nfft_finalize(&my_full_plan); 
    }
  else
    printf("nan\n");
  fflush(stdout);
}

/**********************************************************/
/* main                                                   */
/**********************************************************/
int main(int argc, char *argv[])
{
  int J,default_J,max_J;
  int m,default_m,max_m;
  int max_J_NFFT,max_J_SNDFT;

  //test_snfft_3d(2,10);
  //fftw_cleanup();
  //exit(-1);

  /* test things d=2 */
  default_m=4;
  default_J=7;

  max_m=16;
  max_J=15;

  max_J_NFFT=10;
  max_J_SNDFT=12;

  printf("d=2; default_m_2=%d; default_J_2=%d; max_m_2=%d; max_J_2=%d;\n",
	 default_m,default_J,max_m,max_J);
  printf("Error_2d=[\n");
  for(m=2; m<=max_m; m++)
    {
      printf("%d\t",m);
      test_error(default_J,m);
    }
  printf("];\n");

  /*
  printf("Time_2d=[\n");
  for(J=2; J<=max_J; J++)
    {
      printf("%d\t",J); 
      if(J<=max_J_NFFT)
	test_time(J,default_m,1,1);
      else
	{
	  if(J<=max_J_SNDFT)
	    test_time(J,default_m,0,1);
	  else
	    test_time(J,default_m,0,0);
	}
    }
  printf("];\n");
  
  printf("Memory_2d=[\n");
  for(J=2; J<=max_J; J++)
    {
      printf("%d\t",J);
      if(J<=max_J_NFFT)
	test_memory(J,default_m,1);
      else
	test_memory(J,default_m,0);
    }
  printf("];\n");
  */

  /* test things d=3 */
  default_m=4;
  default_J=5;

  max_m=16;
  max_J=10;

  max_J_NFFT=5;
  max_J_SNDFT=8;

    printf("d=3; default_m_3=%d; default_J_3=%d; max_m_3=%d; max_J_3=%d;\n",
	 default_m,default_J,max_m,max_J);

  printf("Error_3d=[\n");
  for(m=1; m<=max_m; m++)
    {
      printf("%d\t",m);
      test_error_3d(default_J,m,0);
    }
  printf("];\n");
  /*
  printf("Time_3d=[\n");
  for(J=2; J<=max_J; J++)
    {
      printf("%d\t",J);
     
      if(J<=max_J_NFFT)
	test_time_3d(J,default_m,1,1);
      else
	{
	  if(J<=max_J_SNDFT)
	    test_time_3d(J,default_m,0,1);
	  else
	    test_time_3d(J,default_m,0,0);
	}
    }
  printf("];\n");
  
  printf("Memory_3d=[\n");
  for(J=2; J<=max_J; J++)
    {
      printf("%d\t",J);
      if(J<=max_J_NFFT)
	test_memory_3d(J,default_m,1);
      else
	test_memory_3d(J,default_m,0);
    }
  printf("];\n");
  */

  return(1);
}
