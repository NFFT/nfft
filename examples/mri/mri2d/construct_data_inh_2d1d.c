#include "nfft3.h"
#include "util.h"
#include "math.h"

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
  double W;
  int N3;

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
    fscanf(finh,"%le ",&w);
    if(w<min_inh)
      min_inh = w;
    if(w>max_inh)
      max_inh = w;
  }
  fclose(finh);


  W=2.0*MAX(fabs(min_inh),fabs(max_inh)); 
  //N3=ceil((MAX(fabs(min_inh),fabs(max_inh))*(max_time-min_time)+6/(2*2))*4*2);
  N3=ceil((max_time-min_time)*(MAX(fabs(min_inh),fabs(max_inh))+6/(2*2))*4*2);

  
  fprintf(stderr,"3:  %i %e %e %e %e %e %e\n",N3,W,min_inh,max_inh,min_time,max_time,Ts);
  


  my_N[0]=N; my_n[0]=ceil(N*2);
  my_N[1]=N; my_n[1]=ceil(N*2);
  my_N[2]=N3; my_n[2]=N3;
  
  /* initialise nfft */ 
  mri_inh_2d1d_init_guru(&my_plan, my_N, M, my_n, 6,flags,
                      FFTW_MEASURE| FFTW_DESTROY_INPUT);

  ftime=fopen("readout_time.dat","r");
  fp=fopen("knots.dat","r");
  
  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fp,"%le %le ",&my_plan.x[2*j+0],&my_plan.x[2*j+1]);
    fscanf(ftime,"%le ",&my_plan.t[j]);
    my_plan.t[j] = (my_plan.t[j]-Ts)*W/N3;
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
    printf("usage: ./construct_data FILENAME N M\n");
    return 1;
  }
  
  nfft(argv[1],atoi(argv[2]),atoi(argv[3]));
  
  return 1;
}
