#include "nfft3.h"
#include "util.h"

#include <stdlib.h>
#include <stdio.h>

#define FILENAME_LENGTH_MAX 50
#define REPEAT 1
const char TESTFILES_DPT[] = "dpt.txt\0";
const char TESTFILES_DPT_TRANSPOSED[] = "dpt_transposed.txt\0";
#define THRESHOLD 1000.0

void test_dpt_trafo(void)
{
  /** The set of DPTs */
  fpt_set set;
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
  complex* y_ref;
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
  /** For time measurements. */
  double time;

  /* Tell what we're doing. */
  fprintf(stdout,"dpt_trafo: Testing dpt_trafo ...\n");

  /* Try to open file containing the names of the test data files. */
  testfiles = fopen(TESTFILES_DPT,"r");

  fprintf(stdout,"%d\n",testfiles);
  fflush(stdout);

  /* Test if successful. */
  if (testfiles == NULL)
  {
    fprintf(stderr,"Couldn't open %s to read test data filenames!\n");
    return;
  }

  while (fscanf(testfiles,"%s",filename) == 1)
  {
    fprintf(stdout,"filename = %s\t,",filename);
    /* Open input file. */
    file = fopen(filename,"r");
    /* Check if file was opened successfully. */
    if (file != NULL)
    {
      /* Read in DPT mode */
      fscanf(file,"%d",&function_values);
      fprintf(stdout," function_values = %1d,",function_values);

      /* Read in transfrom length. */
      fscanf(file,"%d",&t);
      N = 1<<t;
      fprintf(stdout," t = %2d,",t);
      fprintf(stdout," N = %4d,",N);

      /* Read in start index. */
      fscanf(file,"%d",&k_start);
      fprintf(stdout," k_start = %4d,",k_start);

      /* Read in end index. */
      fscanf(file,"%d",&k_end);
      fprintf(stdout," k_end = %4d,",k_end);

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
      for (k = 0; k < N+2; k++)
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
      }

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
     fprintf(stdout,"\n Legendre coeffs \n");
      for (k = k_start; k <= k_end; k++)
      {
        fprintf(stdout,"x[%d] = %le + I*%le\n",k,creal(x[k]),cimag(x[k]));
      }

      /* Allocate memory for Chebyshev coefficients. */
      y = (complex*) calloc((k_end+1),sizeof(complex));
      y_ref = (complex*) calloc((k_end+1),sizeof(complex));

      /* Read in Chebyshev coefficients. */
      for (k = 0; k <= k_end; k++)
      {
        fscanf(file,"%le",&d1);
        fscanf(file,"%le",&d2);
        y_ref[k] = d1 + I*d2;
      }

      /* Print out Chebyshev coefficients. */
     fprintf(stdout,"\n Chebychev coeffs \n");
      for (k = 0; k <= k_end; k++)
      {
        fprintf(stdout,"y_ref[%d] = %le + I*%le\n",k,creal(y_ref[k]),cimag(y_ref[k]));
      }

      /* Initialize DPT. */
      //fprintf(stderr,"t = %d -> N = %d\n",t,1<<t);
      //fflush(stderr);
      set = fpt_init(0,t,0U/*FPT_AL_SYMMETRY*/);

      /* Precompute DPT. */
      fpt_precompute(set,0,alpha,beta,gamma,k_start,THRESHOLD);

      /* Execute DPT. */
      time = nfft_second();
      for (k = 0; k < REPEAT; k++)
      {
        dpt_trafo(set,0,&x[k_start],y,k_end,0U | (function_values?FPT_FUNCTION_VALUES:0U));
      }
      time = (nfft_second() - time)/((double)REPEAT);

      /* Print out computed and reference coefficients. */
      fprintf(stdout,"\ncomputed and reference:\n");
      for (k = 0; k <= k_end; k++)
      {
        fprintf(stdout,"y_ref[%d] = %+1.10le + I*%+1.10le, \t y[%d] = %+1.10le + I*%+1.10le\n",k,
          creal(y_ref[k]),cimag(y_ref[k]),k,creal(y[k]),cimag(y[k]));
      }

      /* Print out the infinity-norm error. */
      fprintf(stdout," e_infty = %11le,",error_l_infty_complex(y_ref,y,k_end+1));
      fprintf(stdout," e_2 = %11le",error_l_2_complex(y_ref,y,k_end+1));
      fprintf(stdout," time = %11le",time);

      /* CLose the file. */
      fclose(file);
      file = NULL;

      /* Forget precomputed data. */
      fpt_finalize(set);
      set = NULL;

      /* Free memory. */
      free(alpha);
      free(beta);
      free(gamma);
      free(x);
      free(y);
      free(y_ref);
      alpha = NULL;
      beta = NULL;
      gamma = NULL;
      x = NULL;
      y = NULL;
      y_ref = NULL;
      //free(f_orig);
      /* Test passed. */
      fprintf(stdout,"\n");
    }
    else
    {
      fprintf(stdout," failed: Couldn't open file %s.\n",filename);
    }
  }
  close(testfiles);
  testfiles = NULL;
}

void test_dpt_transposed(void)
{
  /** The set of DPTs */
  fpt_set set;
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
  /** The reference Legendre coefficients */
  complex* x;
  complex* x_ref;
  /** The Chebyshev coefficients */
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
  /** For time measurements. */
  double time;

  /* Tell what we're doing. */
  fprintf(stdout,"dpt_transposed: Testing dpt_transposed ...\n");

  /* Try to open file containing the names of the test data files. */
  testfiles = fopen(TESTFILES_DPT_TRANSPOSED,"r");

  fprintf(stdout,"%d\n",testfiles);
  fflush(stdout);

  /* Test if successful. */
  if (testfiles == NULL)
  {
    fprintf(stderr,"Couldn't open %s to read test data filenames!\n");
    return;
  }

  while (fscanf(testfiles,"%s",filename) == 1)
  {
    fprintf(stdout,"filename = %s\t,",filename);
    /* Open input file. */
    file = fopen(filename,"r");
    /* Check if file was opened successfully. */
    if (file != NULL)
    {
      /* Read in DPT mode */
      fscanf(file,"%d",&function_values);
      fprintf(stdout," function_values = %1d,",function_values);

      /* Read in transfrom length. */
      fscanf(file,"%d",&t);
      N = 1<<t;
      fprintf(stdout," t = %2d,",t);
      fprintf(stdout," N = %4d,",N);

      /* Read in start index. */
      fscanf(file,"%d",&k_start);
      fprintf(stdout," k_start = %4d,",k_start);

      /* Read in end index. */
      fscanf(file,"%d",&k_end);
      fprintf(stdout," k_end = %4d,",k_end);

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

      /* Allocate memory for reference Legendre coefficients. */
      x = (complex*) calloc((k_end+1),sizeof(complex));
      x_ref = (complex*) calloc((k_end+1),sizeof(complex));

      /* Read in Legendre coefficients. */
      for (k = k_start; k <= k_end; k++)
      {
        fscanf(file,"%le",&d1);
        fscanf(file,"%le",&d2);
        x_ref[k] = d1 + I*d2;
      }

      /* Print out reference Legendre coefficients. */
      /*fprintf(stdout,"\n Legendre coefficients \n");
      for (k = k_start; k <= k_end; k++)
      {
        fprintf(stdout,"x_ref[%d] = %le + I*%le\n",k,creal(x_ref[k]),cimag(x_ref[k]));
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
      /*fprintf(stdout,"\n");
      for (k = 0; k <= k_end; k++)
      {
        fprintf(stdout,"y[%d] = %le + I*%le\n",k,creal(y[k]),cimag(y[k]));
      }*/

      /* Initialize DPT. */
      set = fpt_init(0,t,0U/*FPT_AL_SYMMETRY*/);

      /* Precompute DPT. */
      fpt_precompute(set,0,alpha,beta,gamma,k_start,THRESHOLD);

      /* Execute DPT. */
      time = nfft_second();
      for (k = 0; k < REPEAT; k++)
      {
        dpt_transposed(set,0,&x[k_start],y,k_end, 0U | (function_values?FPT_FUNCTION_VALUES:0U));
      }
      time = (nfft_second() - time)/((double)REPEAT);

      /* Print out computed and reference coefficients. */
      /*fprintf(stdout,"\n");
      for (k = k_start; k <= k_end; k++)
      {
        fprintf(stdout,"x[%d] = %+1.16le + I*%+1.16le, \t x_ref[%d] = %+1.16le + I*%+1.16le\n",k,
          creal(x[k]),cimag(x[k]),k,creal(x_ref[k]),cimag(x_ref[k]));
      }*/

      /* Print out the infinity-norm error. */
      fprintf(stdout," e_infty = %11le,",error_l_infty_complex(&x_ref[k_start],&x[k_start],k_end-k_start+1));
      fprintf(stdout," e_2 = %11le",error_l_2_complex(&x_ref[k_start],&x[k_start],k_end-k_start+1));
      fprintf(stdout," time = %11le",time);

      /* CLose the file. */
      fclose(file);
      file = NULL;

      /* Forget precomputed data. */
      fpt_finalize(set);
      set = NULL;

      /* Free memory. */
      free(alpha);
      free(beta);
      free(gamma);
      free(x);
      free(x_ref);
      free(y);
      alpha = NULL;
      beta = NULL;
      gamma = NULL;
      x = NULL;
      x_ref = NULL;
      y = NULL;
      //free(f_orig);
      /* Test passed. */
      fprintf(stdout,"\n");
    }
    else
    {
      fprintf(stdout," failed: Couldn't open file %s.\n",filename);
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
  /* Add test for ndsft_trafo. */
  test_dpt_trafo();
  /* Add test for ndsft_trafo. */
  test_dpt_transposed();
  /* Exit the program. */
  return EXIT_SUCCESS;
}
