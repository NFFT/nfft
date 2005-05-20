#include "utils.h"
#include "nfft.h"

/**
 * nfft makes an 2d-adjoint-nfft for every slice
 */
void nfft (char* filename,int N,int M,int Z, int weight ,fftw_complex *mem)
{
  int j,k,z;               /* some variables  */
  double weights;          /* store one weight temporary */
  double tmp;              /* tmp to read the obsolent z from the input file */
  nfft_plan my_plan;       /* plan for the two dimensional nfft  */
  int my_N[2],my_n[2];     /* to init the nfft */
  FILE* fin;               /* input file  */
  FILE* fweight;           /* input file for the weights */
  
  /* initialise my_plan */
  my_N[0]=N; my_n[0]=next_power_of_2(N);
  my_N[1]=N; my_n[1]=next_power_of_2(N);  
  nfft_init_specific(&my_plan, 2, my_N, M, my_n, 6, PRE_PHI_HUT| PRE_PSI|
                        MALLOC_X| MALLOC_F_HAT| MALLOC_F| 
                        FFTW_INIT| FFT_OUT_OF_PLACE| FFTW_MEASURE| FFTW_DESTROY_INPUT,
                        FFTW_MEASURE| FFTW_DESTROY_INPUT);
  
  /* precompute lin psi if set */
  if(my_plan.nfft_flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(&my_plan,10000);

  fin=fopen(filename,"r");

  for(z=0;z<Z;z++) {
    fweight=fopen("weights.dat","r");
    for(j=0;j<my_plan.M;j++)
    {
      fscanf(fweight,"%le ",&weights);
      fscanf(fin,"%le %le %le %le %le",
             &my_plan.x[2*j+0],&my_plan.x[2*j+1],&tmp,&my_plan.f[j][0],&my_plan.f[j][1]);
      if(weight) {
        my_plan.f[j][0] = my_plan.f[j][0] * weights;
        my_plan.f[j][1] = my_plan.f[j][1] * weights;
      }
    } 
    fclose(fweight);
    
    /* precompute psi if set just one time because the nodes equal each slice */
    if(z==0 && my_plan.nfft_flags & PRE_PSI)
      nfft_precompute_psi(&my_plan);
    
    /* precompute full psi if set just one time because the nodes equal each slice */    
    if(z==0 && my_plan.nfft_flags & PRE_FULL_PSI)
      nfft_full_psi(&my_plan,pow(10,-15));
    
    /* compute the adjoint nfft */
    nfft_adjoint(&my_plan);

    for(k=0;k<my_plan.N_L;k++) {
      /* write every slice in the memory.
      here we make an fftshift direct */
      mem[(Z*N*N/2+z*N*N+ k)%(Z*N*N)][0] = my_plan.f_hat[k][0];
      mem[(Z*N*N/2+z*N*N+ k)%(Z*N*N)][1] = my_plan.f_hat[k][1];
    }
  }
  fclose(fin);
  
  nfft_finalize(&my_plan);
}

 /**
 * print writes the memory back in a file
 * output_real.dat for the real part and output_imag.dat for the imaginary part
 */
void print(int N,int M,int Z, fftw_complex *mem)
{
  int i,j;
  FILE* fout_real;
  FILE* fout_imag;
  fout_real=fopen("output_real.dat","w");
  fout_imag=fopen("output_imag.dat","w");

  for(i=0;i<Z;i++) {
    for (j=0;j<N*N;j++) {
      fprintf(fout_real,"%le ",mem[i*N*N+j][0] /Z);
      fprintf(fout_imag,"%le ",mem[i*N*N+j][1] /Z);
    }
    fprintf(fout_real,"\n");
    fprintf(fout_imag,"\n");
  }

  fclose(fout_real);
  fclose(fout_imag);
}


int main(int argc, char **argv)
{
  fftw_complex *mem;
  fftw_plan plan;
  int N,M,Z;  

  if (argc <= 6) {
    printf("usage: ./reconstruct_data_gridding FILENAME N M Z ITER WEIGHTS\n");
    return 1;
  }

  N=atoi(argv[2]);
  M=atoi(argv[3]);
  Z=atoi(argv[4]);

  /* Allocate memory to hold every slice in memory after the
  2D-infft */
  mem = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * atoi(argv[2]) * atoi(argv[2]) * atoi(argv[4]));

  /* Create plan for the 1d-ifft */
  plan = fftw_plan_many_dft(1, &Z, N*N,
                                  mem, NULL,
                                  N*N, 1,
                                  mem, NULL,
                                  N*N,1 ,
                                  FFTW_BACKWARD, FFTW_MEASURE);
 
  /* execute the 2d-infft */
  nfft(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[6]),mem);

  /* execute the 1d-fft */
  fftw_execute(plan);
  
  /* write the memory back in files */
  print(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]), mem);

  /* free memory */
  fftw_free(mem);

  return 1;
}

