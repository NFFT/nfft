#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "util.h"
#include "nfft3.h"

/** 
 * \defgroup applications_mri2d_construct_data_2d construct_data_2d
 * \ingroup applications_mri2d
 * \{
 */

/**
 * construct makes an 2d-nfft
 */
void construct(char * file, int N, int M)
{
  int j,k;            /* some variables */
  double real;
  nfft_plan my_plan;  /* plan for the two dimensional nfft  */
  FILE* fp;
  FILE* fk;
  FILE* fi;
  
  /* initialise my_plan */
  nfft_init_2d(&my_plan,N,N,M);

  fp=fopen("knots.dat","r");
  
  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fp,"%le %le ",&my_plan.x[2*j+0],&my_plan.x[2*j+1]);
  }
  fclose(fp);

  fi=fopen("input_f.dat","r");
  fk=fopen(file,"w");

  for(j=0;j<N;j++)
  {
    for(k=0;k<N;k++) {
      fscanf(fi,"%le ",&real);
      my_plan.f_hat[(N*j+k)] = real;
    }
  }
    
  if(my_plan.nfft_flags & PRE_PSI)
    nfft_precompute_psi(&my_plan);

  nfft_trafo(&my_plan);

  for(j=0;j<my_plan.M_total;j++)
  {
    fprintf(fk,"%le %le %le %le\n",my_plan.x[2*j+0],my_plan.x[2*j+1],creal(my_plan.f[j]),cimag(my_plan.f[j]));
  }
  fclose(fk);
  fclose(fi);

  nfft_finalize(&my_plan);
}

int main(int argc, char **argv)
{ 
  if (argc <= 3) {
    printf("usage: ./construct_data FILENAME N M\n");
    return 1;
  }
  
  construct(argv[1],atoi(argv[2]),atoi(argv[3]));
  
  return 1;
}
/* \} */
