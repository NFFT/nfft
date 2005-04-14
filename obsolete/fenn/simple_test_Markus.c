#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <malloc.h>
#include <sys/resource.h>
#include <time.h>
#include <string.h>

#include "nfft.h"
#include "window_defines.h"
#include "utils.h"

/*========================================================*/
#ifndef sparse_nfft_plan_h
#define sparse_nfft_plan_h

typedef struct sparse_nfft_plan_ 
{
  int J;                      
  int sigma;
  int N_S;
  int M;

  int r_act_nfft_plan;
  nfft_plan* act_nfft_plan;
  nfft_plan* center_nfft_plan;

  fftw_plan* set_fftw_plan1;
  fftw_plan* set_fftw_plan2;

  fftw_complex *f_hat;
  fftw_complex *f;

  double **c_phi_inv;

  fftw_complex* dummy_nfft_f_hat;
} sparse_nfft_plan;

#endif

int int_2_pow(int a)
{
  return (1U<< a);
}

/*========================================================*/

/* copies this_sparse_plan->f_hat to this_plan->f_hat */
void copy_sparse_to_full(sparse_nfft_plan *this_sparse_plan, nfft_plan *this_full_plan)
{
 int r;
  int k1, k2;
  int a,b;
  const int J=this_sparse_plan->J;   /* N=2^J                  */
  const int N=this_full_plan->N[0];  /* size of full NFFT      */
  const int N_B=int_2_pow(J-2);      /* size of small blocks   */

  /* initialize f_hat with zero values */
  memset(this_full_plan->f_hat, 0, N*N*sizeof(fftw_complex));
  
   /* copy values at hyperbolic grid points */
  for (r=0; r<ceil(J/2.0); r++){
    a=int_2_pow(J-r-2); b=int_2_pow(r);
    for (k1=0; k1<a; k1++){
      for (k2=0; k2<b; k2++){
        /* top */
        this_full_plan->f_hat[( k1+a+N/2 )*N+(k2-b/2+N/2)][0]=this_sparse_plan->f_hat[(4*r+0)*N_B+ k1*b+k2 ][0];
        this_full_plan->f_hat[( k1+a+N/2 )*N+(k2-b/2+N/2)][1]=this_sparse_plan->f_hat[(4*r+0)*N_B+ k1*b+k2 ][1];
        /* bottom */
        this_full_plan->f_hat[(k1-2*a+N/2)*N+(k2-b/2+N/2)][0]=this_sparse_plan->f_hat[(4*r+1)*N_B+ k1*b+k2 ][0];
        this_full_plan->f_hat[(k1-2*a+N/2)*N+(k2-b/2+N/2)][1]=this_sparse_plan->f_hat[(4*r+1)*N_B+ k1*b+k2 ][1];
        /* right */
        this_full_plan->f_hat[(k2-b/2+N/2)*N+( k1+a+N/2 )][0]=this_sparse_plan->f_hat[(4*r+2)*N_B+ k2*a+k1 ][0];
        this_full_plan->f_hat[(k2-b/2+N/2)*N+( k1+a+N/2 )][1]=this_sparse_plan->f_hat[(4*r+2)*N_B+ k2*a+k1 ][1];
        /* left */
        this_full_plan->f_hat[(k2-b/2+N/2)*N+(k1-2*a+N/2)][0]=this_sparse_plan->f_hat[(4*r+3)*N_B+ k2*a+k1 ][0];
        this_full_plan->f_hat[(k2-b/2+N/2)*N+(k1-2*a+N/2)][1]=this_sparse_plan->f_hat[(4*r+3)*N_B+ k2*a+k1 ][1];
      }
    }
  }
  /* copy values at center points */
  a=pow(2,floor(J/2));
  for (k1=0; k1<a; k1++){
    for (k2=0; k2<a; k2++){
      this_full_plan->f_hat[(k1-a/2+N/2)*N+(k2-a/2+N/2)][0]=this_sparse_plan->f_hat[4*(int)ceil(J/2.0)*N_B + k1*a+k2][0];
      this_full_plan->f_hat[(k1-a/2+N/2)*N+(k2-a/2+N/2)][1]=this_sparse_plan->f_hat[4*(int)ceil(J/2.0)*N_B + k1*a+k2][1];
    }
  }

  /* copy nodes */
  memcpy(this_full_plan->x,this_sparse_plan->act_nfft_plan->x,this_sparse_plan->M*2*sizeof(double));
}

/* test copy_sparse_to_full */
void test_copy_sparse_to_full(sparse_nfft_plan *this_sparse_plan, nfft_plan *this_full_plan)
{
  int r;
  int k1, k2;
  int a,b;
  const int J=this_sparse_plan->J;   /* N=2^J                  */
  const int N=this_full_plan->N[0];  /* size of full NFFT      */
  const int N_B=int_2_pow(J-2);      /* size of small blocks   */

  /* copy sparse plan to full plan */
  copy_sparse_to_full(this_sparse_plan, this_full_plan);

  /* show blockwise f_hat */
  printf("f_hat blockwise\n");
  for (r=0; r<ceil(J/2.0); r++){
    a=int_2_pow(J-r-2); b=int_2_pow(r);

    printf("top\n");
    for (k1=0; k1<a; k1++){
      for (k2=0; k2<b; k2++){
        printf("(%1.1f,%1.1f) ", this_sparse_plan->f_hat[(4*r+0)*N_B+ k1*b+k2 ][0]
	                       , this_sparse_plan->f_hat[(4*r+0)*N_B+ k1*b+k2 ][1]);
      }
      printf("\n");
    }

    printf("bottom\n");
    for (k1=0; k1<a; k1++){
      for (k2=0; k2<b; k2++){
        printf("(%1.1f,%1.1f) ", this_sparse_plan->f_hat[(4*r+1)*N_B+ k1*b+k2 ][0]
                              , this_sparse_plan->f_hat[(4*r+1)*N_B+ k1*b+k2 ][1]);
      }
      printf("\n");
    }

    printf("right\n");
    for (k2=0; k2<b; k2++){
      for (k1=0; k1<a; k1++){
        printf("(%1.1f,%1.1f) ", this_sparse_plan->f_hat[(4*r+2)*N_B+ k2*a+k1 ][0]
                               , this_sparse_plan->f_hat[(4*r+2)*N_B+ k2*a+k1 ][1]);
      }
      printf("\n");
    }

    printf("left\n");
    for (k2=0; k2<b; k2++){
      for (k1=0; k1<a; k1++){
        printf("(%1.1f,%1.1f) ", this_sparse_plan->f_hat[(4*r+3)*N_B+ k2*a+k1 ][0]
	                       , this_sparse_plan->f_hat[(4*r+3)*N_B+ k2*a+k1 ][1]);
      }
      printf("\n");
    }
  }

  /* show full f_hat */
  printf("full f_hat\n");
  for (k1=0;k1<N;k1++){
    for (k2=0;k2<N;k2++){
      printf("(%1.1f,%1.1f) ", this_full_plan->f_hat[k1*N+k2][0],this_full_plan->f_hat[k1*N+k2][1]);
    }
    printf("\n");
  }
}

void sparse_init_random_nodes_coeffs(sparse_nfft_plan *this_sparse_plan)
{
  int k_S,j;

  /* init frequencies */
  for(k_S=0;k_S<this_sparse_plan->N_S;k_S++)
    {
      this_sparse_plan->f_hat[k_S][0]=0;//(double)rand()/RAND_MAX;
      this_sparse_plan->f_hat[k_S][1]=0;//(double)rand()/RAND_MAX;
    }

  this_sparse_plan->f_hat[4*((this_sparse_plan->J+1)/2)*int_2_pow(this_sparse_plan->J-2)+
			  int_2_pow(2*(this_sparse_plan->J/2)-1)+int_2_pow(this_sparse_plan->J/2-1)-1][0]=1.0;

  /* init nodes */
  //srand((unsigned)time(NULL));
  srand(1);
  for(j=0;j<this_sparse_plan->M;j++) 
    {
      this_sparse_plan->act_nfft_plan->x[2*j+0]=((double)rand())/RAND_MAX-0.5;
      this_sparse_plan->act_nfft_plan->x[2*j+1]=((double)rand())/RAND_MAX-0.5;
    }
}

/*========================================================*/
/* J >1, no precomputation at all!! */
void sparse_nfft_init(sparse_nfft_plan *this_plan, int J, int M, int m)
{
  int r;
  int N[2];
  int n[2];

  this_plan->sigma=2;
  this_plan->J=J;
  this_plan->M=M;
  this_plan->N_S=(J+2)*int_2_pow(J-1);
  
  /* memory allocation */
  this_plan->f = (fftw_complex *)fftw_malloc(M*sizeof(fftw_complex));
  this_plan->f_hat = (fftw_complex *)fftw_malloc(this_plan->N_S*sizeof(fftw_complex));

  this_plan->act_nfft_plan = (nfft_plan*)fftw_malloc(sizeof(nfft_plan));
  this_plan->center_nfft_plan = (nfft_plan*)fftw_malloc(sizeof(nfft_plan));

  this_plan->set_fftw_plan1=(fftw_plan*) fftw_malloc((J-1)*sizeof(fftw_plan));
  this_plan->set_fftw_plan2=(fftw_plan*) fftw_malloc((J-1)*sizeof(fftw_plan));

  /* planning the small nffts */
  /* r=0 */
  N[0]=int_2_pow(J-2); n[0]=this_plan->sigma*N[0];
  N[1]=1; n[1]=this_plan->sigma*N[1];
  
  nfft_init_specific(this_plan->act_nfft_plan,2,N,M,n,m, 0, FFTW_ESTIMATE);
  this_plan->dummy_nfft_f_hat=this_plan->act_nfft_plan->f_hat;

  this_plan->set_fftw_plan1[0]=this_plan->act_nfft_plan->my_fftw_plan1;
  this_plan->set_fftw_plan2[0]=this_plan->act_nfft_plan->my_fftw_plan2;

  for(r=1;r<J-1;r++)
    {
      N[0]=int_2_pow(J-r-2); n[0]=this_plan->sigma*N[0];
      N[1]=int_2_pow(r); n[1]=this_plan->sigma*N[1];
      this_plan->set_fftw_plan1[r] = 
	fftw_plan_dft(2, n, this_plan->act_nfft_plan->g1, this_plan->act_nfft_plan->g2,
		      FFTW_FORWARD, this_plan->act_nfft_plan->fftw_flags);
      this_plan->set_fftw_plan2[r] = 
	fftw_plan_dft(2, n, this_plan->act_nfft_plan->g2, this_plan->act_nfft_plan->g1,
		      FFTW_BACKWARD, this_plan->act_nfft_plan->fftw_flags);
    }

  /* center plan */
  /* J/2 == floor(((double)J) / 2.0) */
  N[0]=int_2_pow(J/2); n[0]=this_plan->sigma*N[0];
  N[1]=int_2_pow(J/2); n[1]=this_plan->sigma*N[1];
  nfft_init_specific(this_plan->center_nfft_plan,2,N,M,n, m, 0, FFTW_ESTIMATE);
}

/*========================================================*/
void sparse_nfft_finalize(sparse_nfft_plan *this_plan)
{
  int r;

  /* center plan */
  nfft_finalize(this_plan->center_nfft_plan);
  
  /* finalize the small nffts */
  this_plan->act_nfft_plan->my_fftw_plan1=this_plan->set_fftw_plan1[0];
  this_plan->act_nfft_plan->my_fftw_plan2=this_plan->set_fftw_plan2[0];

  for(r=1;r<this_plan->J-1;r++)
    {
      fftw_destroy_plan(this_plan->set_fftw_plan1[r]);
      fftw_destroy_plan(this_plan->set_fftw_plan2[r]);
    }

  /* r=0 */
  this_plan->act_nfft_plan->f_hat=this_plan->dummy_nfft_f_hat;
  nfft_finalize(this_plan->act_nfft_plan);

  fftw_free(this_plan->set_fftw_plan1);
  fftw_free(this_plan->set_fftw_plan2);

  fftw_free(this_plan->f_hat);
  fftw_free(this_plan->f);
}


void sparse_nfft_trafo(sparse_nfft_plan *this_plan)
{
  int r,j;
  double temp,ctx,cty,stx,sty;

  int M=this_plan->M;
  int J=this_plan->J;

  //memset(this_plan->f,0,M*sizeof(fftw_complex));
  for(j=0;j<M;j++)
    {
      this_plan->f[j][0]=0.0;
      this_plan->f[j][1]=0.0;
    }

  for(r=0;r<ceil(J/2.0);r++)
    {
      temp=-2.0*PI*3.0*ceil(pow(2.0,J-r-3));

      /* top */
      this_plan->act_nfft_plan->my_fftw_plan1 = this_plan->set_fftw_plan1[r];
      this_plan->act_nfft_plan->N[0]=int_2_pow(J-r-2); this_plan->act_nfft_plan->n[0]=this_plan->sigma*this_plan->act_nfft_plan->N[0];
      this_plan->act_nfft_plan->N[1]=int_2_pow(r); this_plan->act_nfft_plan->n[1]=this_plan->sigma*this_plan->act_nfft_plan->N[1];

      this_plan->act_nfft_plan->f_hat=this_plan->f_hat+4*r*int_2_pow(J-2)+0;

      if((this_plan->act_nfft_plan->N[0]<=this_plan->act_nfft_plan->m)||
	 (this_plan->act_nfft_plan->N[1]<=this_plan->act_nfft_plan->m))
	ndft_trafo(this_plan->act_nfft_plan);
      else
	nfft_trafo(this_plan->act_nfft_plan);

      for (j=0; j<M; j++)
	{
	  ctx=cos(temp*this_plan->act_nfft_plan->x[2*j+0]); stx=sin(temp*this_plan->act_nfft_plan->x[2*j+0]);
	  //cty=cos(temp*this_plan->act_nfft_plan->x[2*j+1]); sty=sin(temp*this_plan->act_nfft_plan->x[2*j+1]);

	  this_plan->f[j][0]+= (this_plan->act_nfft_plan->f[j][0]*ctx - this_plan->act_nfft_plan->f[j][1]*stx);
	  this_plan->f[j][1]+= (this_plan->act_nfft_plan->f[j][1]*ctx + this_plan->act_nfft_plan->f[j][0]*stx);
	}
      
      /* bottom */
      this_plan->act_nfft_plan->f_hat=this_plan->f_hat+4*r*int_2_pow(J-2)+1;
      
      if((this_plan->act_nfft_plan->N[0]<=this_plan->act_nfft_plan->m)||
	 (this_plan->act_nfft_plan->N[1]<=this_plan->act_nfft_plan->m))
	ndft_trafo(this_plan->act_nfft_plan);
      else
	nfft_trafo(this_plan->act_nfft_plan);

      for (j=0; j<M; j++)
	{
	  ctx=cos(temp*this_plan->act_nfft_plan->x[2*j+0]); stx=sin(temp*this_plan->act_nfft_plan->x[2*j+0]);
	  //cty=cos(temp*this_plan->act_nfft_plan->x[2*j+1]); sty=sin(temp*this_plan->act_nfft_plan->x[2*j+1]);

	  this_plan->f[j][0]+= (this_plan->act_nfft_plan->f[j][0]*ctx + this_plan->act_nfft_plan->f[j][1]*stx);
	  this_plan->f[j][1]+= (this_plan->act_nfft_plan->f[j][1]*ctx - this_plan->act_nfft_plan->f[j][0]*stx);
	}

      /* right */
      this_plan->act_nfft_plan->my_fftw_plan1 = this_plan->set_fftw_plan1[J-r-2];
      this_plan->act_nfft_plan->N[0]=int_2_pow(r); this_plan->act_nfft_plan->n[0]=this_plan->sigma*this_plan->act_nfft_plan->N[0];
      this_plan->act_nfft_plan->N[1]=int_2_pow(J-r-2); this_plan->act_nfft_plan->n[1]=this_plan->sigma*this_plan->act_nfft_plan->N[1];

      this_plan->act_nfft_plan->f_hat=this_plan->f_hat+4*r*int_2_pow(J-2)+2;
      
      if((this_plan->act_nfft_plan->N[0]<=this_plan->act_nfft_plan->m)||
	 (this_plan->act_nfft_plan->N[1]<=this_plan->act_nfft_plan->m))
	ndft_trafo(this_plan->act_nfft_plan);
      else
	nfft_trafo(this_plan->act_nfft_plan);

      for (j=0; j<M; j++)
	{
	  //ctx=cos(temp*this_plan->act_nfft_plan->x[2*j+0]); stx=sin(temp*this_plan->act_nfft_plan->x[2*j+0]);
	  cty=cos(temp*this_plan->act_nfft_plan->x[2*j+1]); sty=sin(temp*this_plan->act_nfft_plan->x[2*j+1]);

	  this_plan->f[j][0]+= (this_plan->act_nfft_plan->f[j][0]*cty - this_plan->act_nfft_plan->f[j][1]*sty);
	  this_plan->f[j][1]+= (this_plan->act_nfft_plan->f[j][1]*cty + this_plan->act_nfft_plan->f[j][0]*sty);
	}

      /* left */
      this_plan->act_nfft_plan->f_hat=this_plan->f_hat+4*r*int_2_pow(J-2)+3;
      
      if((this_plan->act_nfft_plan->N[0]<=this_plan->act_nfft_plan->m)||
	 (this_plan->act_nfft_plan->N[1]<=this_plan->act_nfft_plan->m))
	ndft_trafo(this_plan->act_nfft_plan);
      else
	nfft_trafo(this_plan->act_nfft_plan);

      for (j=0; j<M; j++)
	{
	  //ctx=cos(temp*this_plan->act_nfft_plan->x[2*j+0]); stx=sin(temp*this_plan->act_nfft_plan->x[2*j+0]);
	  cty=cos(temp*this_plan->act_nfft_plan->x[2*j+1]); sty=sin(temp*this_plan->act_nfft_plan->x[2*j+1]);

	  this_plan->f[j][0]+= (this_plan->act_nfft_plan->f[j][0]*cty + this_plan->act_nfft_plan->f[j][1]*sty);
	  this_plan->f[j][1]+= (this_plan->act_nfft_plan->f[j][1]*cty - this_plan->act_nfft_plan->f[j][0]*sty);
	}
    }


  /* Was ist mit den FFTW Plänen? */

  vpr_c(this_plan->f,this_plan->M,"vor center");


  /* center */
  r=4*((J+1)/2);
 
  //  this_plan->center_nfft_plan->f_hat=this_plan->f_hat+r*int_2_pow(J-2);
  memcpy(this_plan->center_nfft_plan->f_hat,this_plan->f_hat+r*int_2_pow(J-2),
	 this_plan->center_nfft_plan->N[0]*this_plan->center_nfft_plan->N[0]*sizeof(fftw_complex));

  /*  if (this_plan->center_nfft_plan->N[0]<=this_plan->center_nfft_plan->m) 
      ndft_trafo(this_plan->center_nfft_plan);
      else
      nfft_trafo(this_plan->center_nfft_plan);*/
  
  vpr_c(this_plan->center_nfft_plan->f_hat,this_plan->center_nfft_plan->N_L,"f_hat_center");
  ndft_trafo(this_plan->center_nfft_plan);
  vpr_c(this_plan->center_nfft_plan->f,this_plan->center_nfft_plan->M,"f_center");

  for (j=0; j<M; j++) 
    {
      this_plan->f[j][0]+= this_plan->center_nfft_plan->f[j][0];
      this_plan->f[j][1]+= this_plan->center_nfft_plan->f[j][1];
    }
}


/**********************************************************/
/* main                                                   */
/**********************************************************/
int main(int argc, char *argv[])
{
  int J, k, M;
  int my_N;
  int m=6;

  double t;                           /* time measurement */
  int N[2];
  int n[2];

  nfft_plan my_full_plan;
  sparse_nfft_plan my_sparse_plan;

  fftw_complex *slow_full, *slow_sparse;

  /* evaluate parameters */
  if (argc<2) {
    printf("simple_test J\n");
    exit(-1);
  }

  J=atoi(argv[1]);  
  my_N=int_2_pow(J);
  M=J*my_N;
  N[0]=my_N; N[1]=my_N; n[0]=2*my_N; n[1]=2*my_N;

  printf("N=[%d,%d], n=[%d,%d], M=%d\n",N[0],N[1],n[0],n[1],M);

  slow_full   = (fftw_complex*) fftw_malloc(my_full_plan.M*sizeof(fftw_complex));
  slow_sparse = (fftw_complex*) fftw_malloc(my_full_plan.M*sizeof(fftw_complex));

  /* initialise my_full_plan */
  sparse_nfft_init(&my_sparse_plan, J, M, m);

  nfft_init_specific(&my_full_plan,2,N,M,n,m,0, FFTW_ESTIMATE);

  /* init frequencies */
  sparse_init_random_nodes_coeffs(&my_sparse_plan);
  copy_sparse_to_full(&my_sparse_plan, &my_full_plan);
  //test_copy_sparse_to_full(&my_sparse_plan, &my_full_plan);

  ndft_trafo(&my_full_plan);
  vpr_c(my_full_plan.f,my_full_plan.M,"bla");

  sparse_nfft_trafo(&my_sparse_plan);
  vpr_c(my_sparse_plan.f,my_sparse_plan.M,"bla");

  printf("Error (fast sparse-fast full  ): %e\n",E_2_error_c(my_sparse_plan.f,my_full_plan.f,M));

  nfft_finalize(&my_full_plan);
  sparse_nfft_finalize(&my_sparse_plan);
  fftw_free(slow_full);
  fftw_free(slow_sparse);

  return 1;
}
