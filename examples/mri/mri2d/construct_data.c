#include "util.h"
#include "nfft3.h"

/**
 * nfft makes an 2d-nfft
 */
void nfft (char * file, int N, int M)
{
  int j,k;            /* some variables */
  nfft_plan my_plan;  /* plan for the two dimensional nfft  */
  FILE* fp;
  FILE* fk;
  FILE* fi;
  
  /* initialise my_plan */
  nfft_init_2d(&my_plan,N,N,M);

  fp=fopen("nodes.dat","r");
  
  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fp,"%le %le ",&my_plan.x[2*j+0],&my_plan.x[2*j+1]);
  }
  fclose(fp);

  fi=fopen("input_f.dat","r");
  fk=fopen(file,"w");

  for(j=0;j<N;j++)
  {
    for(k=0;k<N;k++)
      fscanf(fi,"%le ",&my_plan.f_hat[(N*j+k)]);
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
  
  nfft(argv[1],atoi(argv[2]),atoi(argv[3]));
  
  return 1;
}
