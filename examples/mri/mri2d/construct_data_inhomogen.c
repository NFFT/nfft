#include "options.h"
#include "util.h"
#include "nfft3.h"
#include "window_defines.h"

#define PHI_periodic(x) ((x>0.5)?(PHI(x-1.0,0)):((x<-0.5)?PHI(x+1.0,0):PHI(x,0)))

/**
 * nfft makes an 2d-nfft
 */
void nfft (char * file, int N, int M)
{
  int j,k,l;                  /* some variables */
  double real,w;
  nfft_plan my_plan,*ths;  /* plan for the two dimensional nfft  */
  FILE *fp,*fout,*fi,*finh;

  double Te = 2.0;
  double Aq = M * 0.00402542373;
  double Ts = (Aq + 2.0 * Te) / 2.0;
  double W = 1.3 ;
  int N3 = ceil(W*Aq);

  ths = (nfft_plan*) malloc(sizeof(nfft_plan));

  nfft_init_1d(ths,N3,1);
  N3=ths->n[0];

  
  /* initialise my_plan */
  nfft_init_3d(&my_plan,N3,N,N,M);

  fp=fopen("nodes.dat","r");
  
  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fp,"%le %le ",&my_plan.x[3*j+1],&my_plan.x[3*j+2]);
    //my_plan.x[3*j+0] =  W*( (((double)(j%(M/arms))*0.00402542373)) + Te  - Ts)  / ((double)N3);
  }
  fclose(fp);


  for(l=0;l<N3;l++)
  {
    fi=fopen("input_f.dat","r");
    finh=fopen("inh.dat","r");
    for(j=0;j<N;j++)
    {
      for(k=0;k<N;k++)
      {
        fscanf(fi,"%le ",&real);
        fscanf(finh,"%le ",&w);
        my_plan.f_hat[l*N*N+N*j+k] = real*PHI_periodic(w/W-((double)l)/((double)N3)+0.5);
      }
    }
    fclose(fi);
    fclose(finh);
  }
    
  if(my_plan.nfft_flags & PRE_PSI)
    nfft_precompute_psi(&my_plan);

  nfft_trafo(&my_plan);

  fout=fopen(file,"w");
  
  for(j=0;j<my_plan.M_total;j++)
  {
    my_plan.f[j]/=PHI_HUT(N3*my_plan.x[3*j+0],0);
    fprintf(fout,"%le %le %le %le\n",my_plan.x[2*j+0],my_plan.x[2*j+1],creal(my_plan.f[j]),cimag(my_plan.f[j]));
  }

  fclose(fout);

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
