#include <malloc.h>

#include "nfft.h"
#include "utils.h"
#include "short_dim_2d.h"

#define MAX(a,b) ((a>b)?a:b)
#define MIN(a,b) ((a<b)?a:b)
#define SWAP(a,b) {void* c=(void*)a; a=b; b=c;}

#ifndef sparse_nfft_plan_h
#define sparse_nfft_plan_h

#define SNDFT		(1U<< 0)

typedef struct sparse_nfft_plan_ 
{
  int J;                      
  int sigma;
  int N_S;
  int M;

  unsigned snfft_flags;

  int* index_sparse_to_full;

  int r_act_nfft_plan;
  nfft_plan* act_nfft_plan;
  nfft_plan* center_nfft_plan;

  fftw_plan* set_fftw_plan1;
  //fftw_plan* set_fftw_plan2;

  nfft_plan* set_nfft_plan_1d;

  double* x_transposed;
  
  fftw_complex *f_hat;
  fftw_complex *f;

} sparse_nfft_plan;

typedef struct sparse_nfft_plan_3d_ 
{
  int J;                      
  int sigma;
  int N_S;
  int M;

  unsigned snfft_flags;

  int* index_sparse_to_full;         /**< overflow for J=9!        */

  int r_act_nfft_plan;
  nfft_plan* act_nfft_plan;
  nfft_plan* center_nfft_plan;

  fftw_plan* set_fftw_plan1;
  //fftw_plan* set_fftw_plan2;

  nfft_plan* set_nfft_plan_1d;
  nfft_plan* set_nfft_plan_2d;

  double *x_102,*x_201,*x_120,*x_021;	
  
  fftw_complex *f_hat;
  fftw_complex *f;

} sparse_nfft_plan_3d;

#endif

int int_2_pow(int a);
int total_used_memory();

void sparse_init_random_nodes_coeffs(sparse_nfft_plan *this_sparse_plan);
void copy_sparse_to_full(sparse_nfft_plan *this_sparse_plan, nfft_plan *this_full_plan);

void sparse_nfft_init(sparse_nfft_plan *this, int J, int M, int m, unsigned snfft_flags);
void sparse_nfft_finalize(sparse_nfft_plan *this);

void sparse_nfft_trafo(sparse_nfft_plan *this);
void sparse_ndft_trafo(sparse_nfft_plan *this);


void sparse_init_random_nodes_coeffs_3d(sparse_nfft_plan_3d *this);
void copy_sparse_to_full_3d(sparse_nfft_plan_3d *this, nfft_plan *this_full_plan);

void sparse_ndft_trafo_3d(sparse_nfft_plan_3d *this);
void sparse_nfft_trafo_3d(sparse_nfft_plan_3d *this);

void sparse_nfft_init_3d(sparse_nfft_plan_3d *this, int J, int M, int m, unsigned snfft_flags);
void sparse_nfft_finalize_3d(sparse_nfft_plan_3d *this);

