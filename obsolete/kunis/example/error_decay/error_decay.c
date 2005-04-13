#include "utils.h"
#include "nfft.h"

#define DIRICHLET 0
#define FEJER 1
#define JACKSON2 2
#define JACKSON4 3
#define SOBOLEV 4
#define MULTIQ 5

double native_error(fftw_complex *x0, fftw_complex *x, double *w,int n)
{
  double xd=0;
  fftw_complex xda;
  int k;
  
  for (k=0; k<n; k++)
    {
      xda[0]=(double)fabs(x0[k][0]-x[k][0]);
      xda[1]=(double)fabs(x0[k][1]-x[k][1]);
      xd+=w[k]*(xda[0]*xda[0]+xda[1]*xda[1]);
    }
  
  return (double)sqrt(xd);
}

void error_decay(int N,int M,int kernel)
{
  int j,k,l;                            /* index nodes, freqencies, iterations*/
  nfft_plan my_plan;                    /* plan for the one dimensional nfft  */
  infft_plan my_iplan;                  /* plan for the one dimensional infft */
  FILE* fp;                             /* file input_data                    */
  double temp1,temp2,temp3;             /* for reading input                  */
  fftw_complex* f_hat_solution;         /* known solution                     */
   
  /* allocate memory for the known solution */
  f_hat_solution=(fftw_complex*)fftw_malloc(N*sizeof(fftw_complex));
  
  /* initialise my_plan */
  nfft_init_1d(&my_plan,N,M);

  /* initialise my_iplan */
  infft_init_specific(&my_iplan,&my_plan, CGNE_R | PRECOMPUTE_DAMP);

  /* init nodes */
  fp=fopen("input_data.dat","r");
  for(j=0;j<my_plan.M;j++)
    {
      fscanf(fp,"%le %le %le",&temp1,&temp2,&temp3);
      my_plan.x[j]=temp1;
      my_plan.f[j][0]=temp2;
      my_plan.f[j][1]=temp3;
    }
  fclose(fp);
  
  /* precompute psi */
  if(my_plan.nfft_flags & PRE_PSI)
    nfft_precompute_psi(&my_plan);

   /* initialise damping factors */
  if(my_iplan.infft_flags & PRECOMPUTE_DAMP)
    for(k=0;k<my_plan.N_L;k++)
      {
	if(kernel==DIRICHLET)
	  my_iplan.w_hat[k]=1.0;
	if(kernel==FEJER)
	  my_iplan.w_hat[k]=modified_fejer(N,k-N/2);
	if(kernel==JACKSON2)
	  my_iplan.w_hat[k]=modified_jackson2(N,k-N/2);
	if(kernel==JACKSON4)
	  my_iplan.w_hat[k]=modified_jackson4(N,k-N/2);
	if(kernel==SOBOLEV)
	  my_iplan.w_hat[k]=modified_sobolev(2,k-N/2);
	if(kernel==MULTIQ)
	  my_iplan.w_hat[k]=modified_multiquadric(2,1,k-N/2);
      }

  /* compute the exact solution f_hat_solution from f_tilde (my_plan.f) */
  nfft_adjoint(&my_plan);
  if(my_iplan.infft_flags & PRECOMPUTE_DAMP)
    copyc_w(my_plan.f_hat,my_iplan.w_hat,my_plan.f_hat,my_plan.N_L);
  
  /* compute the right hand side */
  nfft_trafo(&my_plan);

  /* save the exact solution */
  SWAPC(f_hat_solution,my_plan.f_hat);

  /* initialise the right hand side */
  SWAPC(my_plan.f,my_iplan.given_f);

  /* init some guess */
  for(k=0;k<my_plan.N_L;k++)
    {
      my_iplan.f_hat_iter[k][0]=0.0;
      my_iplan.f_hat_iter[k][1]=0.0;
    }

  /* inverse trafo */  
  infft_before_loop(&my_iplan);
  for(l=0;l<100;l++)
    {
      infft_loop_one_step(&my_iplan);

      if(my_iplan.infft_flags & PRECOMPUTE_DAMP)
	printf("%e\n",native_error(f_hat_solution,my_iplan.f_hat_iter,my_iplan.w_hat,my_plan.N_L));
    }

  /* finalize */
  infft_finalize(&my_iplan);  
  nfft_finalize(&my_plan);  
  fftw_free(f_hat_solution);
}

int main(int argc, char **argv)
{
  error_decay(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));

  return 1;
}
