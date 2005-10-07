#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "util.h"
#include "nfft3.h"



double norm2( int k0, int k1)
{
  return sqrt(k0*k0+k1*k1);
}

void func_rec_2d( int my_N,
                  int M,
                  int iter,
                  int my_M,
                  double border_eps,
                  double mu, double c, char *fout_name) 
{

  int j,j0,j1,k,k0,k1,l;                      /* nodes, freqencies iterations       */
  nfst_plan  my_cplan, my_other_cplan;        /* plans for the two dimensional nfst */
  infst_plan my_icplan;                       /* plan for the two dimensional infst */

  FILE *fp_in, *fp_out;                       /* input/output file                  */
  double min_x, max_x, min_y, max_y;


  printf( "N=%d  M=%d  iter=%d  rec_M=%d  b_eps=%f  mu=%f  c=%f  fout=%s\n", 
      my_N, M, iter, my_M, border_eps, mu, c, fout_name);

  nfst_init_2d(&my_cplan, my_N, my_N, M);

  /* initialise my_icplan, specific */
  infst_init_advanced( &my_icplan, &my_cplan, CGNE | PRECOMPUTE_DAMP);


  /* init nodes */
  if( (fp_in = fopen( "input.dat", "r")) == NULL)
  {
    fprintf( stderr, "%s\n", "Can't open inputfile");
    exit( 1);
  }
      
  
  fscanf( fp_in, "%le %le %le", &my_cplan.x[0], &my_cplan.x[1], &my_icplan.y[0]);
  min_x = max_x = my_cplan.x[0];
  min_y = max_y = my_cplan.x[1];

  j=0;
  for( j = 1; j < my_cplan.M_total; j++)
  {
      fscanf(fp_in, "%le %le %le", &my_cplan.x[2*j+0],
                                   &my_cplan.x[2*j+1],
                                   &my_icplan.y[j]);

      min_x = (my_cplan.x[2*j+0] < min_x) ? my_cplan.x[2*j+0] : min_x;
      max_x = (my_cplan.x[2*j+0] > max_x) ? my_cplan.x[2*j+0] : max_x;

      min_y = (my_cplan.x[2*j+1] < min_y) ? my_cplan.x[2*j+1] : min_y;
      max_y = (my_cplan.x[2*j+1] > max_y) ? my_cplan.x[2*j+1] : max_y;
  }

  fclose(fp_in);

  // scale input
  double diff_x = 1.0/(0.5-2*border_eps) * (max_x - min_x);
  double diff_y = 1.0/(0.5-2*border_eps) * (max_y - min_y);
  for( j = 0; j < my_cplan.M_total; j++)
  {
    my_cplan.x[2*j+0] = border_eps + (my_cplan.x[2*j+0] - min_x) / diff_x;
    my_cplan.x[2*j+1] = border_eps + (my_cplan.x[2*j+1] - min_y) / diff_y;
  }

  /* precompute psi */
  if(my_cplan.nfst_flags & PRE_PSI)
    nfst_precompute_psi( &my_cplan);


  /** initialise damping factors */
  if(my_icplan.flags & PRECOMPUTE_DAMP)
    {
      for(k0=0;k0<my_cplan.N[0]-1;k0++)
        for(k1=0;k1<my_cplan.N[1]-1;k1++)
        {
          my_icplan.w_hat[k0*(my_cplan.N[1]-1)+k1]=
            modified_multiquadric( mu, c, norm2(k0, k1));
        }
      
//      for(k0=0;k0<my_cplan.N[0]-1;k0++)
//        my_icplan.w_hat[k0*(my_cplan.N[1]-1)+(my_cplan.N[1]-1)-1]=0;
//      for(k1=0;k1<my_cplan.N[1]-1;k1++)
//        my_icplan.w_hat[((my_cplan.N[0]-1)-1)*(my_cplan.N[1]-1)+k1]=0;    
    }


  /* init some guess */
  for( k = 0; k < my_cplan.N_total; k++)
    my_icplan.f_hat_iter[k] = 0.0;

  /* inverse trafo */  
  infst_before_loop( &my_icplan);
  for(l=0; l < iter; l++)
  { 
      infst_loop_one_step( &my_icplan);
      //fprintf( stderr, "%e\n", my_icplan.dot_r_iter);
  }



  nfst_init_2d( &my_other_cplan, my_N, my_N, my_M * my_M);
  
  // nodes on grid
  for( j0 = 0; j0 < my_M; j0++) 
    for( j1 = 0; j1 < my_M; j1++)
    {
      my_other_cplan.x[2*(j0 * my_M + j1) + 0] = 
	border_eps + ((double)((double)j0*((0.5-2*border_eps)/my_M)));
      my_other_cplan.x[2*(j0 * my_M + j1) + 1] = 
	border_eps + ((double)((double)j1*((0.5-2*border_eps)/my_M)));
    }


  /* precompute psi */
  if( my_other_cplan.nfst_flags & PRE_PSI)
    nfst_precompute_psi( &my_other_cplan);

  SWAP_double( my_icplan.f_hat_iter, my_other_cplan.f_hat);
  nfst_trafo( &my_other_cplan);
  SWAP_double( my_icplan.f_hat_iter, my_other_cplan.f_hat);

  
  if( (fp_out = fopen( fout_name, "w")) == NULL)
  {
    fprintf( stderr, "%s\n", "Can't open outputfile");
    exit( 2);
  }
    
  for( j = 0; j < my_other_cplan.M_total; j++)
  {
    fprintf( fp_out, "%f  %f  %f\n", my_other_cplan.x[2*j], 
                                     my_other_cplan.x[2*j+1],
                                     my_other_cplan.f[j]);
  }

  fclose( fp_out);

  infst_finalize( &my_icplan);  
  nfst_finalize( &my_cplan);
  nfst_finalize( &my_other_cplan);
}



int main(int argc, char **argv)
{
  func_rec_2d( atoi( argv[1]), /* N_0 */
               atoi( argv[2]), /* M   */
               atoi( argv[3]), /* iterations */
               atoi( argv[4]), /* rec_M */
               atof( argv[5]), /* border_eps */
               atof( argv[6]), /* mu (multiquadric damping factors) */
               atof( argv[7]), /*  c (multiquadric damping factors) */
                     argv[8]   /* output filename */
               );

  return 1;
}
