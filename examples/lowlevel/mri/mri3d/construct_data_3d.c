#include "util.h"
#include "nfft3.h"

void construct(char * file, int N, int M, int Z)
{
  int j,k,l;
  double real;
  nfft_plan my_plan;
  FILE* fp;
  FILE* fk;


  nfft_init_3d(&my_plan,Z,N,N,M);

  fp=fopen("knots.dat","r");
   
  for(j=0;j<M;j++)
    fscanf(fp,"%le %le %le",&my_plan.x[3*(j)+1],
      &my_plan.x[3*(j)+2],&my_plan.x[3*(j)+0]);

  fclose(fp);

  fp=fopen("input_f.dat","r");
  fk=fopen(file,"w");

  for(l=0;l<Z;l++) {
    for(j=0;j<N;j++)
    {
      for(k=0;k<N;k++)
      {
        //fscanf(fp,"%le ",&my_plan.f_hat[(N*N*(Z-l)+N*j+k+N*N*Z/2)%(N*N*Z)][0]);
        fscanf(fp,"%le ",&real);
        my_plan.f_hat[(N*N*l+N*j+k)] = real;
      }
    }
  }
    
    if(my_plan.nfft_flags & PRE_PSI)
      nfft_precompute_psi(&my_plan);

    nfft_trafo(&my_plan);

    
    for(j=0;j<my_plan.M_total;j++)
      fprintf(fk,"%le %le %le %le %le\n",my_plan.x[3*j+1],
      my_plan.x[3*j+2],my_plan.x[3*j+0],creal(my_plan.f[j]),cimag(my_plan.f[j]));

    
  
  fclose(fk);
  fclose(fp);

  nfft_finalize(&my_plan);
}

int main(int argc, char **argv)
{ 
  if (argc <= 4) {
    printf("usage: ./construct_data FILENAME N M Z\n");
    return 1;
  }

  construct(argv[1], atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));
  
	return 1;
}
