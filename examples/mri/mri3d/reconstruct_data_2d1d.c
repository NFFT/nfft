#include "util.h"
#include "nfft3.h"
#include "malloc.h"

/**
 * infft makes an inverse 2d-nfft for every slice
 */
void infft(char* filename,int N,int M,int Z,int iteration, int weight, fftw_complex *mem)
{
  int j,k,l,z;                  /* some variables  */
  double real,imag;             /* to read the real and imag part of a complex number */
  nfft_plan my_plan;            /* plan for the two dimensional nfft  */
  infft_plan my_iplan;          /* plan for the two dimensional infft */
  FILE* fin;                    /* input file */
  int my_N[2],my_n[2];          /* to init the nfft */
  double tmp, epsilon=0.0000003;/* tmp to read the obsolent z from the input file
                                   epsilon is the break criterium for
                                   the iteration */
  unsigned infft_flags = CGNR; /* flags for the infft */
                                   
  /* initialise my_plan */
  my_N[0]=N; my_n[0]=2*next_power_of_2(N);
  my_N[1]=N; my_n[1]=2*next_power_of_2(N);
  nfft_init_guru(&my_plan, 2, my_N, M/Z, my_n, 6, PRE_PHI_HUT| PRE_PSI|
                         MALLOC_X| MALLOC_F_HAT| MALLOC_F| 
                        FFTW_INIT| FFT_OUT_OF_PLACE| FFTW_MEASURE| FFTW_DESTROY_INPUT,
                        FFTW_MEASURE| FFTW_DESTROY_INPUT);
  
  /* precompute lin psi if set */
  if(my_plan.nfft_flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(&my_plan);
  
  /* set the flags for the infft*/
  if (weight)
    infft_flags = infft_flags | PRECOMPUTE_WEIGHT;
  
  /* initialise my_iplan, advanced */
  infft_init_advanced(&my_iplan,&my_plan, infft_flags );

  /* get the weights */
  if(my_iplan.flags & PRECOMPUTE_WEIGHT)
  {
    fin=fopen("weights.dat","r");
    for(j=0;j<my_plan.M_total;j++) 
    {
        fscanf(fin,"%le ",&my_iplan.w[j]);
    }
    fclose(fin);
  }
  
  /* open the input file */
  fin=fopen(filename,"r"); 

  /* For every Layer*/
  for(z=0;z<Z;z++) {

    /* read x,y,freal and fimag from the knots */
    for(j=0;j<my_plan.M_total;j++)
    {
      fscanf(fin,"%le %le %le %le %le ",&my_plan.x[2*j+0],&my_plan.x[2*j+1], &tmp,
      &real,&imag);
      my_iplan.y[j] = real + I*imag;
    }
    
    /* precompute psi if set just one time because the knots equal each plane */
    if(z==0 && my_plan.nfft_flags & PRE_PSI) 
      nfft_precompute_psi(&my_plan);
      
    /* precompute full psi if set just one time because the knots equal each plane */
    if(z==0 && my_plan.nfft_flags & PRE_FULL_PSI)
      nfft_precompute_full_psi(&my_plan);

    /* init some guess */
    for(k=0;k<my_plan.N_total;k++)
      my_iplan.f_hat_iter[k]=0.0;
 
    /* inverse trafo */
    infft_before_loop(&my_iplan);
    for(l=0;l<iteration;l++)
    {
      /* break if dot_r_iter is smaller than epsilon*/
      if(my_iplan.dot_r_iter<epsilon)
      break;
      fprintf(stderr,"%e,  %i of %i\n",sqrt(my_iplan.dot_r_iter),
      iteration*z+l+1,iteration*Z);
      infft_loop_one_step(&my_iplan);
    }
    for(k=0;k<my_plan.N_total;k++) {
      /* write every slice in the memory.
      here we make an fftshift direct */
      mem[(Z*N*N/2+z*N*N+ k)%(Z*N*N)] = my_iplan.f_hat_iter[k];
    }
  }
  /* finalize the infft */
  infft_finalize(&my_iplan);
  
  /* finalize the nfft */
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
      fprintf(fout_real,"%le ",creal(mem[(Z*N*N/2+i*N*N+ j)%(Z*N*N)]) /Z);
      fprintf(fout_imag,"%le ",cimag(mem[(Z*N*N/2+i*N*N+ j)%(Z*N*N)]) /Z);
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
    printf("usage: ./reconstruct FILENAME N M Z ITER WEIGHTS\n");
    return 1;
  }

  N=atoi(argv[2]);
  M=atoi(argv[3]);
  Z=atoi(argv[4]);

  /* Allocate memory to hold every layer in memory after the
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
  infft(argv[1],N,M,Z,atoi(argv[5]),atoi(argv[6]),mem);

  /* execute the 1d-fft */
  fftw_execute(plan);

  /* write the memory back in files */
  print(N,M,Z, mem);

  /* free memory */
  fftw_free(mem);
  fftw_destroy_plan(plan);
  return 1;
}
