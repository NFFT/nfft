#include "nfft3.h"
#include "dpt.h"
#include "util.h"

#include <stdlib.h>
#include <stdio.h>

/* Include CUnit header. */
#include <CUnit/CUnit.h>

/** Maximum filename length */
#define FILENAME_LENGTH_MAX 50

/** Name of the file containing the test data filenames for NDSFT. */
const char TESTFILES_DPT[] = "dpt.txt\0";

/* Default threshold. Not important here since NDSFT-algorithm is used. */
#define THRESHOLD 1000.0

void test_dpt_trafo(void)
{
  /** The set of DPTs */
  dpt_set set;
  /** DPT mode */
  int function_values;
  /** The transform length (must be a power of two) */
  int N;
  /** Start index */
  int k_start;
  /** End index */
  int k_end;
  /** The exponent of N */
  int t;
  /** The Legendre coefficients */
  complex* x;
  /** The reference Chebyshev coefficients */
  complex* y;
  /** The degree \f$k\f$ */
  int k;
  /** Auxilliary variables used to read in complex numbers. */
  double d1,d2;
  /** Recursion coefficients \f$\alpha\f$. */
  double *alpha;
  /** Recursion coefficients \f$\alpha\f$. */
  double *beta;
  /** Recursion coefficients \f$\alpha\f$. */
  double *gamma;
  /** The file containing the names of the testdata files. */
  FILE *testfiles;
  /** The file containg the testcase data */
  FILE *file;
  /** Name of file containing test data. */
  char filename[FILENAME_LENGTH_MAX+1];

  /* Tell what we're doing. */
  fprintf(stdout,"dpt_trafo: Testing dpt_trafo ...\n");

  /* Try to open file containing the names of the test data files. */
  testfiles = fopen(TESTFILES_DPT,"r");

  fprintf(stdout,"%d\n",testfiles);
  fflush(stdout);

  /* Test if successful. */
  if (testfiles == NULL)
  {
    CU_FAIL("Couldn't open %s to read test data filenames!\n");      
    return;
  }

  while (fscanf(testfiles,"%s",filename) == 1)
  { 
    fprintf(stdout,"filename = %s,",filename);
    /* Open input file. */ 
    file = fopen(filename,"r");
    /* Check if file was opened successfully. */
    if (file != NULL)
    {
      /* Read in DPT mode */
      fscanf(file,"%d",&function_values);
      fprintf(stdout," function_values = %d,",function_values);      
      
      /* Read in transfrom length. */
      fscanf(file,"%d",&t);
      N = 1<<t;
      fprintf(stdout," t = %d,",t);
      fprintf(stdout," N = %d,",N);
      
      /* Read in start index. */
      fscanf(file,"%d",&k_start);
      fprintf(stdout," k_start = %d,",k_start);
      
      /* Read in end index. */
      fscanf(file,"%d",&k_end);
      fprintf(stdout," k_end = %d,",k_end);

      /* Allocate memory for recursion coefficients. */
      alpha = (double*) malloc((N+2)*sizeof(double));
      beta = (double*) malloc((N+2)*sizeof(double));
      gamma = (double*) malloc((N+2)*sizeof(double));
      
      /* Read in recursion coeffcients. */
      for (k = 0; k < N+2; k++)
      {
        fscanf(file,"%le\n",&alpha[k]);
      }
      for (k = 0; k < N+2; k++)
      {
        fscanf(file,"%le\n",&beta[k]);
      }
      for (k = 0; k < N+2; k++)
      {
        fscanf(file,"%le\n",&gamma[k]);
      }
      
      /* Print out recursion coefficients. */
      /*for (k = 0; k < N+2; k++)
      {
        fprintf(stdout,"alpha_%d^%d = %.16le\n",k-1,k_start,alpha[k]);
      }
      for (k = 0; k < N+2; k++)
      {
        fprintf(stdout,"beta_%d^%d = %.16le\n",k-1,k_start,beta[k]);
      }
      for (k = 0; k < N+2; k++)
      {
        fprintf(stdout,"gamma_%d^%d = %.16le\n",k-1,k_start,gamma[k]);
      }*/
      
      /* Allocate memory for Legendre coefficients. */
      x = (complex*) calloc((k_end+1),sizeof(complex));
      
      /* Read in Legendre coefficients. */
      for (k = k_start; k <= k_end; k++)
      {
        fscanf(file,"%le",&d1);
        fscanf(file,"%le",&d2);
        x[k] = d1 + I*d2;
      }
      
      /* Print out Legendre coefficients. */
      /*for (k = k_start; k <= k_end; k++)
      {
        fprintf(stdout,"x[%d] = %le + I*%le\n",k,creal(x[k]),cimag(x[k]));
      }*/
            
      /* Allocate memory for Chebyshev coefficients. */
      y = (complex*) calloc((k_end+1),sizeof(complex));
      
      /* Read in Chebyshev coefficients. */
      for (k = 0; k <= k_end; k++)
      {
        fscanf(file,"%le",&d1);
        fscanf(file,"%le",&d2);
        y[k] = d1 + I*d2;
      }
      
      /* Print out Chebyshev coefficients. */
      /*for (k = 0; k <= k_end; k++)
      {
        fprintf(stdout,"y[%d] = %le + I*%le\n",k,creal(y[k]),cimag(y[k]));
      }*/            
            
      /* Initialize DPT. */
      set = dpt_init(0,t,0U);
      
      /* Precompute DPT. */
      dpt_precompute(set,0,alpha,beta,gamma,k_start,THRESHOLD);
      
      /* Execute DPT. */
      fpt_trafo(set,0,&x[k_start],x,k_end,0U | (function_values?DPT_FUNCTION_VALUES:0U));   
      
      /* Print out computed and reference coefficients. */
      /*fprintf(stdout,"\n");
      for (k = 0; k <= k_end; k++)
      {
        fprintf(stdout,"x[%d] = %+1.16le + I*%+1.16le, \t y[%d] = %+1.16le + I*%+1.16le\n",k,
          creal(x[k]),cimag(x[k]),k,creal(y[k]),cimag(y[k]));
      }*/
      
      /* Print out the infinity-norm error. */
      fprintf(stdout," e_infty = %le,",error_l_infty_complex(y,x,k_end+1));
      fprintf(stdout," e_2 = %le",error_l_2_complex(y,x,k_end+1));
 
      /* CLose the file. */
      fclose(file);
      file = NULL;
      
      /* Forget precomputed data. */
      dpt_finalize(set);
      set = NULL;
      
      /* Free memory. */
      free(alpha);
      free(beta);
      free(gamma);
      free(x);
      free(y);     
      alpha = NULL;
      beta = NULL;
      gamma = NULL;
      x = NULL;
      y = NULL;
      //free(f_orig);
      /* Test passed. */
      fprintf(stdout," ok\n");
      CU_PASS("ok");
    }
    else
    {
      fprintf(stdout," failed: Couldn't open file %s.\n",filename);
      CU_FAIL("Couldn't open file!\n");
    }
  }
  close(testfiles);
  testfiles = NULL;  
}

/**
 * The main program.
 *
 * \param argc The number of arguments
 * \param argv An array containing the arguments as C-strings
 *
 * \return Exit code 
 */
int main (int argc, char **argv)
{  
  /* Initialise registry. */
  CU_initialize_registry();
  /* Create test suite. */
  CU_pSuite suite = CU_add_suite("NDSFT", NULL, NULL);
  /* Add test for ndsft_trafo. */
  CU_add_test(suite,"NDSFT",test_dpt_trafo);
  /* Run the tests. */
  CU_automated_run_tests();
  /* Cleanup registry. */
  CU_cleanup_registry();   
  return EXIT_SUCCESS;
}
