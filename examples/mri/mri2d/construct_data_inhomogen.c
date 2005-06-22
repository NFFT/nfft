#include "options.h"
#include "util.h"
#include "nfft3.h"
#include "window_defines.h"
#include "math.h"

#define PHI_periodic(x) ((x>0.5)?(PHI(x-1.0,0)):((x<-0.5)?PHI(x+1.0,0):PHI(x,0)))

/**
 * nfft makes an 2d-nfft
 */
void nfft (char * file, int N, int M)
{
  int j,k,l;                  /* some variables */
  double real,imag;
  double *w;
  double time,min_time,max_time,min_inh,max_inh;
  nfft_plan my_plan,*ths;  /* plan for the two dimensional nfft  */
  FILE *fp,*fout,*fi,*finh,*ftime;
  int my_N[3],my_n[3];
  int flags = PRE_PHI_HUT| PRE_PSI |MALLOC_X| MALLOC_F_HAT|
                      MALLOC_F| FFTW_INIT| FFT_OUT_OF_PLACE|
                      FFTW_MEASURE| FFTW_DESTROY_INPUT;

  double Ts;
  double W;
  int N3;
  
  w = (double*) malloc(N*N*sizeof(double));

  ftime=fopen("readout_time.dat","r");
  finh=fopen("inh.dat","r");

  fprintf(stderr,"1\n");
  
  min_time=999999.0; max_time=-9999999.0;//Integer.maxValue!!!!
  for(j=0;j<M;j++)
  {
    fscanf(ftime,"%le ",&time);
    if(time<min_time)
      min_time = time;
    if(time>max_time)
      max_time = time;
  }

  fprintf(stderr,"2\n");
  
  fclose(ftime);
  
  Ts=min_time+max_time/2.0;


  min_inh=999999.0; max_inh=-9999999.0;//Integer.maxValue!!!!
  for(j=0;j<N*N;j++)
  {
    fscanf(finh,"%le ",&w[j]);
    if(w[j]<min_inh)
      min_inh = w[j];
    if(w[j]>max_inh)
      max_inh = w[j];
  }
  fclose(finh);

  W=2.0*MAX(fabs(min_inh),fabs(max_inh))*(1.2); //1.0+m/n!?!?!?!?!?
  N3=ceil(W*(max_time-min_time));

  fprintf(stderr,"3:  %i %e %e %e %e %e %e\n",N3,W,min_inh,max_inh,min_time,max_time,Ts);
  

  ths = (nfft_plan*) malloc(sizeof(nfft_plan));

  nfft_init_1d(ths,N3,1);
  N3=ths->n[0];

  my_N[0]=N3; my_n[0]=ths->n[0];
  my_N[1]=N; my_n[1]=next_power_of_2(N)+16;
  my_N[2]=N;my_n[2]=next_power_of_2(N)+16;
  
  /* initialise nfft */ 
  nfft_init_guru(&my_plan, 3, my_N, M, my_n, 6,flags,
                      FFTW_MEASURE| FFTW_DESTROY_INPUT);

  ftime=fopen("readout_time.dat","r");
  fp=fopen("knots.dat","r");
  
  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fp,"%le %le ",&my_plan.x[3*j+1],&my_plan.x[3*j+2]);
    fscanf(ftime,"%le ",&my_plan.x[3*j+0]);
    my_plan.x[3*j+0] = (my_plan.x[3*j+0]-Ts)*W/N3;
  }
  fclose(fp);
  fclose(ftime);


  for(l=0;l<N3;l++)
  {
    fi=fopen("input_f.dat","r");
    for(j=0;j<N*N;j++)
    {
      fscanf(fi,"%le ",&real);
      my_plan.f_hat[l*N*N+j] = real*PHI_periodic(w[j]/W-((double)l)/((double)N3)+0.5)/cexp(2.0*I*PI*Ts*w[j]);
    }
    fclose(fi);
  }
  
  if(my_plan.nfft_flags & PRE_PSI)
    nfft_precompute_psi(&my_plan);

  nfft_trafo(&my_plan);

  fout=fopen(file,"w");
  
  for(j=0;j<my_plan.M_total;j++)
  {
    my_plan.f[j]/=PHI_HUT(N3*my_plan.x[3*j+0],0);
    fprintf(fout,"%le %le %le %le\n",my_plan.x[3*j+1],my_plan.x[3*j+2],creal(my_plan.f[j]),cimag(my_plan.f[j]));
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
