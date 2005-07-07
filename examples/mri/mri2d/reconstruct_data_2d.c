#include "util.h"
#include "nfft3.h"
#include "math.h"

/**
 * infft makes an inverse 2d nfft
 */
void infft(char* filename,int N,int M,int iteration, int weight)
{
  int j,k,l;                    /* some variables  */
  double real,imag,t;             /* to read the real and imag part of a complex number */
  nfft_plan my_plan;            /* plan for the two dimensional nfft  */
  infft_plan my_iplan;          /* plan for the two dimensional infft */
  FILE* fin;                    /* input file                         */
  FILE* fout_real;              /* output file                        */
  FILE* fout_imag;              /* output file                        */
  int my_N[2],my_n[2];          /* to init the nfft */
  double epsilon=0.0000003;     /* epsilon is a the break criterium for
                                   the iteration */
  unsigned infft_flags = CGNR;  /* flags for the infft*/
                                    
  /* initialise my_plan */
  my_N[0]=N; my_n[0]=ceil(N*1.2);
  my_N[1]=N; my_n[1]=ceil(N*1.2);
  nfft_init_guru(&my_plan, 2, my_N, M, my_n, 6, PRE_PHI_HUT| PRE_PSI|
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

  /* read x,y,freal and fimag from the knots */
  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fin,"%le %le %le %le ",&my_plan.x[2*j+0],&my_plan.x[2*j+1],
    &real,&imag);
    my_iplan.y[j] = real + I*imag;
  }
  
  /* precompute psi */
  if(my_plan.nfft_flags & PRE_PSI)
    nfft_precompute_psi(&my_plan);

  /* precompute full psi */
  if(my_plan.nfft_flags & PRE_FULL_PSI)
      nfft_precompute_full_psi(&my_plan);

  /* init some guess */
  for(k=0;k<my_plan.N_total;k++)
    my_iplan.f_hat_iter[k]=0.0;

  t=second();
    
  /* inverse trafo */
  infft_before_loop(&my_iplan);
  for(l=0;l<iteration;l++)
  {
    /* break if dot_r_iter is smaller than epsilon*/
    if(my_iplan.dot_r_iter<epsilon)
      break;
    fprintf(stderr,"%e,  %i of %i\n",sqrt(my_iplan.dot_r_iter),
    l+1,iteration);
    infft_loop_one_step(&my_iplan);
  }

  
  t=second()-t;
  fprintf(stderr,"time: %e seconds mem: %i \n",t,total_used_memory());

  
  fout_real=fopen("output_real.dat","w");
  fout_imag=fopen("output_imag.dat","w");

  for(k=0;k<my_plan.N_total;k++) {
    fprintf(fout_real,"%le ", creal(my_iplan.f_hat_iter[k]));
    fprintf(fout_imag,"%le ", cimag(my_iplan.f_hat_iter[k]));
  }

  fclose(fout_real);
  fclose(fout_imag);

  /* finalize the infft */
  infft_finalize(&my_iplan);
  
  /* finalize the nfft */
  nfft_finalize(&my_plan);
}

int main(int argc, char **argv)
{
  if (argc <= 5) {
    printf("usage: ./reconstruct_data_2d FILENAME N M ITER WEIGHTS\n");
    return 1;
  }
  
  infft(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]));
   
  return 1;
}
