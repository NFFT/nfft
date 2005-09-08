#include "nfft3.h"
#include "util.h"
#include "math.h"
#include "limits.h"

/**
 * nfft makes an 2d-nfft
 */
void nfft (char * file, int N, int M)
{
  int j,k,l;                  /* some variables */
  double real,imag;
  double w;
  double time,min_time,max_time,min_inh,max_inh;
  mri_inh_2d1d_plan my_plan;  /* plan for the two dimensional nfft  */
  FILE *fp,*fout,*fi,*finh,*ftime;
  int my_N[3],my_n[3];
  int flags = PRE_PHI_HUT| PRE_PSI |MALLOC_X| MALLOC_F_HAT|
                      MALLOC_F| FFTW_INIT| FFT_OUT_OF_PLACE|
                      FFTW_MEASURE| FFTW_DESTROY_INPUT;

  double Ts;
  double W,T;
  int N3;
  int m=2;
  double alpha = 1.25;

  ftime=fopen("readout_time.dat","r");
  finh=fopen("inh.dat","r");

  min_time=INT_MAX; max_time=INT_MIN;
  for(j=0;j<M;j++)
  {
    fscanf(ftime,"%le ",&time);
    if(time<min_time)
      min_time = time;
    if(time>max_time)
      max_time = time;
  }

  fclose(ftime);
  
  Ts=(min_time+max_time)/2.0;

  min_inh=INT_MAX; max_inh=INT_MIN;
  for(j=0;j<N*N;j++)
  {
    fscanf(finh,"%le ",&w);
    if(w<min_inh)
      min_inh = w;
    if(w>max_inh)
      max_inh = w;
  }
  fclose(finh);


  N3=ceil((MAX(fabs(min_inh),fabs(max_inh))*(max_time-min_time)/2.0+m/(2*alpha))*4*alpha);
  T=((max_time-min_time)/2.0)/(0.5-((double) m)/N3);
  W=N3/T;
  
  fprintf(stderr,"3:  %i %e %e %e %e %e %e\n",N3,W,min_inh,max_inh,min_time,max_time,Ts);


  my_N[0]=N; my_n[0]=ceil(N*alpha);
  my_N[1]=N; my_n[1]=ceil(N*alpha);
  my_N[2]=N3; my_n[2]=N3;
  
  /* initialise nfft */ 
  mri_inh_2d1d_init_guru(&my_plan, my_N, M, my_n, m,flags,
                      FFTW_MEASURE| FFTW_DESTROY_INPUT);

  ftime=fopen("readout_time.dat","r");
  fp=fopen("knots.dat","r");
  
  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fp,"%le %le ",&my_plan.x[2*j+0],&my_plan.x[2*j+1]);
    fscanf(ftime,"%le ",&my_plan.t[j]);
    my_plan.t[j] = (my_plan.t[j]-Ts)/T;
  }
  fclose(fp);
  fclose(ftime);

  finh=fopen("inh.dat","r");
  for(j=0;j<N*N;j++)
  {
    fscanf(finh,"%le ",&my_plan.w[j]);
    my_plan.w[j]/=W;
  }
  fclose(finh);


  fi=fopen("input_f.dat","r");
  for(j=0;j<N*N;j++)
  {
    fscanf(fi,"%le ",&real);
    my_plan.f_hat[j] = real*cexp(2.0*I*PI*Ts*my_plan.w[j]*W);
  }
  
  if(my_plan.nfft_flags & PRE_PSI)
    nfft_precompute_psi((nfft_plan*)&my_plan);

  mri_inh_2d1d_trafo(&my_plan);

  fout=fopen(file,"w");
  
  for(j=0;j<my_plan.M_total;j++)
  {
    fprintf(fout,"%le %le %le %le\n",my_plan.x[2*j+0],my_plan.x[2*j+1],creal(my_plan.f[j]),cimag(my_plan.f[j]));
  }

  fclose(fout);

  mri_inh_2d1d_finalize(&my_plan);
}

int main(int argc, char **argv)
{ 
  if (argc <= 3) {
    printf("usage: ./construct_data_inh_2d1d FILENAME N M\n");
    return 1;
  }
  
  nfft(argv[1],atoi(argv[2]),atoi(argv[3]));
  
  return 1;
}
