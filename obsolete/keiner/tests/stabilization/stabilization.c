/**
 * Program template for using the NFSFT library.
 */

/* Include NFSFT header */
#include<nfsft.h>
#include<stdlib.h>
#include<stdio.h>

#define M_MIN 1
#define M_STRIDE 1
#define M_MAX 128

#define T_MIN 1000
#define T_MAX 1000
#define T_STRIDE 1000

/**
 * The main program
 */
int main(int argc, char **argv)
{
  int m_min; 
  int m_max; 
  int m_stride; 
  double t_min; 
  double t_max; 
  double t_stride; 
  int N; /** Next greater power of two with respect to M */

  /** 
   * For the Fourier coefficients. 
   * Layout: f_hat[M+n][k] = a_k^n for n=-M,...,M; k=0,...,M
   */
  complex **f_hat;
  /* NFSFT transform plan */
  nfsft_plan plan;
  
  /** Loop counters */
  int m,n,nstab,ntotal;
  double t;
  char filename[100];
  FILE *file;
  
  if (argc == 1)
  {
    m_min = M_MIN;
    m_max = M_MAX; 
    m_stride = M_STRIDE; 
    t_min = T_MIN; 
    t_max = T_MAX; 
    t_stride = T_STRIDE; 
  }
  else if (argc == 7)
  {
    sscanf(argv[1],"%d",&m_min);
    sscanf(argv[2],"%d",&m_max);
    sscanf(argv[3],"%d",&m_stride);
    sscanf(argv[4],"%lf",&t_min);
    sscanf(argv[5],"%lf",&t_max);
    sscanf(argv[6],"%lf",&t_stride);
  }
  else
  {
    fprintf(stderr,"Stabilization - Stabilization test for NFSFT\n");
    fprintf(stderr,"Usage: Stabilization M_MIN M_MAX M_STRIDE T_MIN T_MAX T_STRIDE\n");
    return -1;
  }  
  
  N = 1<<ngpt(m_max);
  
  /* Allocate memory */
  f_hat = (complex**) malloc((2*m_max+1)*sizeof(complex*));
  for (n = -m_max; n <= m_max; n++)
  {
    /* Length must be a power of two + 1 >= M */
    f_hat[n+m_max] = (complex*) calloc((N+1),sizeof(complex));
  }  
   
  sprintf(filename,"stabilization.dat");   
  file = fopen(filename,"w");
  fclose(file);
  
  for (t = t_min; t <= t_max; t = t * t_stride)
  {
    nfsft_precompute_stab(m_max,t,0U);
    for (m = m_min; m <= m_max; m = m + m_stride)
    {
      /* Compute NFSFT */
      plan = nfsft_init_stab(1, m, 0U, f_hat, 0U, 0U);
      nstab = nfsft_trafo_stab(plan);
      nfsft_get_stat(&nstab,&ntotal);
      nfsft_finalize_stab(plan);
      file = fopen(filename,"a");
      if (file != NULL)
      {
        fprintf(file,"%10.0f %5d %10d %10d %15.10f\n",t,m,nstab,ntotal,ntotal==0?0:((double)nstab)/ntotal);
        fclose(file);
      }  
      fprintf(stdout,"%10.0f %5d %10d %10d %15.10f\n",t,m,nstab,ntotal,ntotal==0?0:((double)nstab)/ntotal);
    }
    nfsft_forget_stab();
  }
  
  for (n = -m_max; n <= m_max; n++)
  {
    free(f_hat[n+m_max]);
  }  
  free(f_hat);
  
  return EXIT_SUCCESS;
}