#include "options.h"
#include "util.h"
#include "nfft3.h"
#include "window_defines.h"
#include "math.h"

/**
 * infft makes an inverse 2d nfft
 */
void infft(char* filename,int N,int M,int iteration, int weight)
{ int j,k,l;                    /* some variables  */
  nnfft_plan my_plan;            /* plan for the two dimensional nfft  */
  innfft_plan my_iplan;          /* plan for the two dimensional infft */
  FILE* fin;                    /* input file                         */
  FILE* finh;
  FILE* ftime;
  FILE* fout_real;              /* output file                        */
  FILE* fout_imag;              /* output file                        */
  int my_N[3],my_n[3];          /* to init the nfft */
  double t,tmp1,tmp2,tmp3,tmp4,epsilon=0.0000003;     /* epsilon is a the break criterium for
                                   the iteration */
  unsigned infft_flags = CGNR; /* flags for the infft*/
  double time,min_time,max_time,min_inh,max_inh;
  double real,imag;
  double *w;

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
  
  Ts=(min_time+max_time)/2.0;


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

  //W=2.0*MAX(fabs(min_inh),fabs(max_inh))*(1.2); //1.0+m/n!?!?!?!?!?
  //N3=2*ceil(W*(max_time-min_time));
  
  N3=ceil((MAX(fabs(min_inh),fabs(max_inh))*(max_time-min_time)/2.0)*4)+1;


  W=MAX(fabs(min_inh),fabs(max_inh))*2.0;
  //W= MAX(fabs(min_inh),fabs(max_inh))/(0.5-6.0/N3);

  
  fprintf(stderr,"3:  %i %e %e %e %e %e %e\n",N3,W,min_inh,max_inh,min_time,max_time,Ts);

  /* initialise my_plan */
  my_N[0]=N;my_n[0]=ceil(N*1.5);
  my_N[1]=N; my_n[1]=ceil(N*1.5);
  my_N[2]=N3; my_n[2]=ceil(N3*1.5);
  nnfft_init_guru(&my_plan, 3, N*N, M, my_N,my_n,6,
        PRE_PSI| PRE_PHI_HUT| MALLOC_X| MALLOC_V| MALLOC_F_HAT| MALLOC_F );
        
  /* precompute lin psi if set */
  if(my_plan.nnfft_flags & PRE_LIN_PSI)
    nnfft_precompute_lin_psi(&my_plan);
  
  /* set the flags for the infft*/
  if (weight)
    infft_flags = infft_flags | PRECOMPUTE_WEIGHT;

  /* initialise my_iplan, advanced */
  innfft_init_advanced(&my_iplan,&my_plan, infft_flags );

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
  ftime=fopen("readout_time.dat","r");

  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fin,"%le %le %le %le ",&my_plan.x[3*j+0],&my_plan.x[3*j+1],&real,&imag);
    my_iplan.y[j]=real+I*imag;
    fscanf(ftime,"%le ",&my_plan.x[3*j+2]);

    my_plan.x[3*j+2] = (my_plan.x[3*j+2]-Ts)*W/N3;
  }

  for(j=0;j<N;j++)
    { 
    for(l=0;l<N;l++)
      {
        my_plan.v[3*(N*j+l)+0]=(((double) j) -(((double) N)/2.0))/((double) N);
        my_plan.v[3*(N*j+l)+1]=(((double) l) -(((double) N)/2.0))/((double) N);
        my_plan.v[3*(N*j+l)+2] = w[N*j+l]/W ;
      }
    }
 
  /* precompute psi */
  if(my_plan.nnfft_flags & PRE_PSI) {
    nnfft_precompute_psi(&my_plan);
    if(my_plan.nnfft_flags & PRE_FULL_PSI)
      nnfft_precompute_full_psi(&my_plan);
  }
  
  if(my_plan.nnfft_flags & PRE_PHI_HUT)
    nnfft_precompute_phi_hut(&my_plan);
    
  /* init some guess */
  for(k=0;k<my_plan.N_total;k++)
  {
    my_iplan.f_hat_iter[k]=0.0;
  }
 
  t=second();
  
  /* inverse trafo */
  innfft_before_loop(&my_iplan);
  for(l=0;l<iteration;l++)
  {
    /* break if dot_r_iter is smaller than epsilon*/
    if(my_iplan.dot_r_iter<epsilon)
    break;
    fprintf(stderr,"%e,  %i of %i\n",sqrt(my_iplan.dot_r_iter),
    l+1,iteration);
    innfft_loop_one_step(&my_iplan);
  }

  t=second()-t;
  fprintf(stderr,"time: %e seconds mem: %i \n",t,total_used_memory());

  fout_real=fopen("output_real.dat","w");
  fout_imag=fopen("output_imag.dat","w");

  for(k=0;k<my_plan.N_total;k++) {

    my_iplan.f_hat_iter[k]*=cexp(2.0*I*PI*Ts*w[k]);

    fprintf(fout_real,"%le ", creal(my_iplan.f_hat_iter[k]));
    fprintf(fout_imag,"%le ", cimag(my_iplan.f_hat_iter[k]));
  }
  

  fclose(fout_real);
  fclose(fout_imag);


  /* finalize the infft */
  innfft_finalize(&my_iplan);
  
  /* finalize the nfft */
  nnfft_finalize(&my_plan);

  free(w);
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
