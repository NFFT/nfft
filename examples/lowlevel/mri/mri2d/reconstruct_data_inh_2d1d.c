#include "nfft3.h"
#include "util.h"
#include "math.h"
#include "limits.h"

void reconstruct(char* filename,int N,int M,int iteration , int weight)
{
  int j,k,l;
  double weights;
  double time,min_time,max_time,min_inh,max_inh;
  double t,real,imag;
  double w,epsilon=0.0000003;     /* epsilon is a the break criterium for
                                   the iteration */;
  mri_inh_2d1d_plan my_plan;
  imri_inh_2d1d_plan my_iplan;
  FILE* fp,*fw,*fout_real,*fout_imag,*finh,*ftime;
  int my_N[3],my_n[3];
  int flags = PRE_PHI_HUT| PRE_PSI |MALLOC_X| MALLOC_F_HAT|
                      MALLOC_F| FFTW_INIT| FFT_OUT_OF_PLACE;
  unsigned infft_flags = CGNR | PRECOMPUTE_DAMP;

  double Ts;
  double W,T;
  int N3;
  int m=2;
  double sigma = 1.25;

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

  N3=ceil((MAX(fabs(min_inh),fabs(max_inh))*(max_time-min_time)/2.0+(m)/(2*sigma))*4*sigma);
	/* N3 has to be even */
	if(N3%2!=0)
	  N3++;
  
	T=((max_time-min_time)/2.0)/(0.5-((double) (m))/N3);
  W=N3/T;
  
  fprintf(stderr,"3:  %i %e %e %e %e %e %e %e\n",N3,W,T,min_inh,max_inh,min_time,max_time,Ts);
  

  my_N[0]=N; my_n[0]=ceil(N*sigma);
  my_N[1]=N; my_n[1]=ceil(N*sigma);
  my_N[2]=N3; my_n[2]=N3;
  
  /* initialise nfft */ 
  mri_inh_2d1d_init_guru(&my_plan, my_N, M, my_n, m, sigma, flags,
                      FFTW_MEASURE| FFTW_DESTROY_INPUT);


  /* precompute lin psi if set */
  if(my_plan.plan.nfft_flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(&my_plan.plan);
                      
  if (weight)
    infft_flags = infft_flags | PRECOMPUTE_WEIGHT;

  /* initialise my_iplan, advanced */
  imri_inh_2d1d_init_advanced(&my_iplan,&my_plan, infft_flags );

  /* get the weights */
  if(my_iplan.flags & PRECOMPUTE_WEIGHT)
  {
    fw=fopen("weights.dat","r");
    for(j=0;j<my_plan.M_total;j++)
    {
        fscanf(fw,"%le ",&my_iplan.w[j]);
    }
    fclose(fw);
  }
                      
  /* get the damping factors */
  if(my_iplan.flags & PRECOMPUTE_DAMP)
  {
    for(j=0;j<N;j++){
      for(k=0;k<N;k++) {
        int j2= j-N/2;
        int k2= k-N/2;
        double r=sqrt(j2*j2+k2*k2);
        if(r>(double) N/2) 
          my_iplan.w_hat[j*N+k]=0.0;
        else
          my_iplan.w_hat[j*N+k]=1.0;
      }   
    }
  }

  fp=fopen(filename,"r");
  ftime=fopen("readout_time.dat","r");

  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fp,"%le %le %le %le",&my_plan.plan.x[2*j+0],&my_plan.plan.x[2*j+1],&real,&imag);
    my_iplan.y[j]=real+I*imag;
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

    
  if(my_plan.plan.nfft_flags & PRE_PSI) {
    nfft_precompute_psi(&my_plan.plan);
  }
  if(my_plan.plan.nfft_flags & PRE_FULL_PSI) {
      nfft_precompute_full_psi(&my_plan.plan);
  } 

  /* init some guess */
  for(j=0;j<my_plan.N_total;j++)
  {
    my_iplan.f_hat_iter[j]=0.0;
  }
 
  t=second();
  
  /* inverse trafo */
  imri_inh_2d1d_before_loop(&my_iplan);
  for(l=0;l<iteration;l++)
  {
    /* break if dot_r_iter is smaller than epsilon*/
    if(my_iplan.dot_r_iter<epsilon)
    break;
    fprintf(stderr,"%e,  %i of %i\n",sqrt(my_iplan.dot_r_iter),
    l+1,iteration);
    imri_inh_2d1d_loop_one_step(&my_iplan);
  }

  t=second()-t;
#ifdef HAVE_MALLINFO
  fprintf(stderr,"time: %e seconds mem: %i \n",t,total_used_memory());
#else
  fprintf(stderr,"time: %e seconds mem: mallinfo not available\n",t);
#endif
  
  fout_real=fopen("output_real.dat","w");
  fout_imag=fopen("output_imag.dat","w");
  
  for (j=0;j<N*N;j++) {
    /* Verschiebung wieder herausrechnen */
    my_iplan.f_hat_iter[j]*=cexp(-2.0*I*PI*Ts*my_plan.w[j]*W);
    
    fprintf(fout_real,"%le ",creal(my_iplan.f_hat_iter[j]));
    fprintf(fout_imag,"%le ",cimag(my_iplan.f_hat_iter[j]));
  }

  fclose(fout_real);
  fclose(fout_imag);
  imri_inh_2d1d_finalize(&my_iplan);
  mri_inh_2d1d_finalize(&my_plan);
}


int main(int argc, char **argv)
{
  if (argc <= 5) {

    printf("usage: ./reconstruct_data_inh_2d1d FILENAME N M ITER WEIGHTS\n");
    return 1;
  }
  
  reconstruct(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]));

  return 1;
}

