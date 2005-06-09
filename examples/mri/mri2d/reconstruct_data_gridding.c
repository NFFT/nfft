#include "util.h"
#include "nfft3.h"

/**
 * nfft makes a 2d-adjoint-nfft
 */
void nfft (char* filename, int N, int M, int weight)
{
  int j;                   /* some variables  */
  double weights;          /* store one weight temporary */
  double real,imag;        /* to read the real and imag part of a complex number */
  nfft_plan my_plan;       /* plan for the two dimensional nfft  */  
  FILE* fin;               /* input file  */
  FILE* fweight;           /* input file for the weights */
  FILE *fout_real;         /* output file  */
  FILE *fout_imag;         /* output file  */
  int my_N[2],my_n[2];
  int flags = PRE_PHI_HUT| PRE_PSI |MALLOC_X| MALLOC_F_HAT|
                      MALLOC_F| FFTW_INIT| FFT_OUT_OF_PLACE|
                      FFTW_MEASURE| FFTW_DESTROY_INPUT;

  /* initialise nfft */ 
  my_N[0]=N; my_n[0]=next_power_of_2(N);
  my_N[1]=N; my_n[1]=next_power_of_2(N);
  nfft_init_guru(&my_plan, 2, my_N, M, my_n, 6,flags,
                      FFTW_MEASURE| FFTW_DESTROY_INPUT);

  fin=fopen(filename,"r");

  fweight=fopen("weights.dat","r");
  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fweight,"%le ",&weights);
    fscanf(fin,"%le %le %le %le",&my_plan.x[2*j+0],&my_plan.x[2*j+1],&real,&imag);
    my_plan.f[j] = real + I*imag;
    if (weight)
      my_plan.f[j] = my_plan.f[j] * weights;
  } 
  fclose(fweight);

  /* precompute psi */
  if(my_plan.nfft_flags & PRE_PSI)
    nfft_precompute_psi(&my_plan);
  
  /* precompute full psi */
  if(my_plan.nfft_flags & PRE_FULL_PSI)
    nfft_precompute_full_psi(&my_plan);

  
  /* compute the adjoint nfft */
  nfft_adjoint(&my_plan);
    
  fout_real=fopen("output_real.dat","w");
  fout_imag=fopen("output_imag.dat","w");
  
  for (j=0;j<N*N;j++) {
    fprintf(fout_real,"%le ",creal(my_plan.f_hat[j]));
    fprintf(fout_imag,"%le ",cimag(my_plan.f_hat[j]));
  }

  fclose(fin);
  fclose(fout_real);
  fclose(fout_imag);
  
  nfft_finalize(&my_plan);
}


int main(int argc, char **argv)
{
  if (argc <= 5) {
    printf("usage: ./reconstruct_data_gridding FILENAME N M ITER WEIGHTS\n");
    return 1;
  }
  nfft(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[5]));

  return 1;
}

