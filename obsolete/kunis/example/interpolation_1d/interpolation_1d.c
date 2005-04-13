#include "utils.h"
#include "nfft.h"

#define DIRICHLET 0
#define FEJER 1
#define JACKSON2 2
#define JACKSON4 3
#define SOBOLEV 4
#define MULTIQ 5

void interpolation_1d(int N, int M, int kernel)
{
  int j,k,l;                            /**< nodes, freqencies, iterations    */
  nfft_plan my_plan;                    /**< plan for nfft                    */
  infft_plan my_iplan;                  /**< plan for infft                   */
  FILE* fp;                             /**< file input_data                  */
   
  /** initialise my_plan */
  nfft_init_1d(&my_plan,N,M);

  /** initialise my_iplan */
  infft_init_specific(&my_iplan,&my_plan, CGNE_R | PRECOMPUTE_DAMP |
		      PRECOMPUTE_WEIGHT);

  /** init nodes and given data */
  fp=fopen("input_data.dat","r");
  for(j=0; j<my_plan.M; j++)
    {
      fscanf(fp,"%le %le",&my_plan.x[j],&my_iplan.given_f[j][0]);
      my_iplan.given_f[j][1]=0.0;
    }
  fclose(fp);
  
  /** precompute psi */
  if(my_plan.nfft_flags & PRE_PSI)
    nfft_precompute_psi(&my_plan);
  
  /** initialise damping factors */
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
 
  /** initialise weights, assumes ordered x_j */
  if(my_iplan.infft_flags & PRECOMPUTE_WEIGHT)
    voronoi_weights_1d(my_iplan.w, my_plan.x, my_plan.M);

  /** init some guess */
  for(k=0; k<my_plan.N_L; k++)
    {
      my_iplan.f_hat_iter[k][0]=0.0;
      my_iplan.f_hat_iter[k][1]=0.0;
    }

  /** inverse trafo */  
  infft_before_loop(&my_iplan);
  for(l=0; l<my_plan.M; l++)
    {
      infft_loop_one_step(&my_iplan);
    }

  /** write output to stdout */
  for(k=0;k<my_plan.N_L;k++)
    printf("%le\t%le\n",my_iplan.f_hat_iter[k][0],my_iplan.f_hat_iter[k][1]);

  /** finalise */
  infft_finalize(&my_iplan);  
  nfft_finalize(&my_plan);  
}

int main(int argc, char **argv)
{
  interpolation_1d(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));

  return 1;
}
