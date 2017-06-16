/*
 * Copyright (c) 2002, 2016 Jens Keiner, Stefan Kunis, Daniel Potts
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
#include "config.h"

#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include "nfft3.h"
#include "infft.h"
#include "imex.h"
#include "fastsum.c"
#include "kernels.c"

#ifdef HAVE_MEXVERSION_C
  #include "mexversion.c"
#endif

#define PLANS_MAX 100 /* maximum number of plans */
#define CMD_LEN_MAX 40 /* maximum length of command argument */

/* global flags */
#define FASTSUM_MEX_FIRST_CALL (1U << 0)
unsigned short gflags = FASTSUM_MEX_FIRST_CALL;

fastsum_plan* plans[PLANS_MAX]; /* plans */

static inline void check_nargs(const int nrhs, const int n, const char* errmsg)
{
  DM(if (nrhs != n)
    mexErrMsgTxt(errmsg);)
}

static inline void check_plan(int i)
{
  DM(if (i < 0 || i >= PLANS_MAX)
    mexErrMsgTxt("Invalid plan");)
  DM(if (plans[i] == 0)
    mexErrMsgTxt("Plan was not initialized or has already been finalized");)
}

static inline void check_plan_source_nodes(int i)
{
  DM(
	if (!(plans[i]->x))
		mexErrMsgTxt("Required to set source nodes first");
	)
}

static inline void check_plan_target_nodes(int i)
{
  DM(
	if (!(plans[i]->y))
		mexErrMsgTxt("Required to set target nodes first");
	)
}

static inline void check_plan_nodes(int i)
{
  DM(
    if (!(plans[i]->x) && !(plans[i]->y))
		mexErrMsgTxt("Required to set source and target nodes first");
	else if (!(plans[i]->x))
		mexErrMsgTxt("Required to set source nodes first");
	else if (!(plans[i]->y))
		mexErrMsgTxt("Required to set target nodes first");
	)
}

static inline int mkplan()
{
  int i = 0;
  while (i < PLANS_MAX && plans[i] != 0) i++;
  if (i == PLANS_MAX)
    mexErrMsgTxt("fastsum: Too many plans already allocated.");
  plans[i] = nfft_malloc(sizeof(fastsum_plan));
  return i;
}

static kernel get_kernel(const mxArray *p)
{
	kernel ker;
	char s[CMD_LEN_MAX+1]; /**< name of kernel          */
	if (mxGetString(p, s, CMD_LEN_MAX))
		mexErrMsgTxt("Could not get kernel string.");
    if (strcmp(s, "gaussian") == 0)
      ker = gaussian;
    else if (strcmp(s, "multiquadric") == 0)
      ker = multiquadric;
    else if (strcmp(s, "inverse_multiquadric") == 0)
      ker = inverse_multiquadric;
    else if (strcmp(s, "logarithm") == 0)
      ker = logarithm;
    else if (strcmp(s, "thinplate_spline") == 0)
      ker = thinplate_spline;
    else if (strcmp(s, "one_over_square") == 0)
      ker = one_over_square;
    else if (strcmp(s, "one_over_modulus") == 0)
      ker = one_over_modulus;
    else if (strcmp(s, "one_over_x") == 0)
      ker = one_over_x;
    else if (strcmp(s, "inverse_multiquadric3") == 0)
      ker = inverse_multiquadric3;
    else if (strcmp(s, "sinc_kernel") == 0)
      ker = sinc_kernel;
    else if (strcmp(s, "cosc") == 0)
      ker = cosc;
    else if (strcmp(s, "cot") == 0)
      ker = kcot;
    else if (strcmp(s, "one_over_cube") == 0)
      ker = one_over_cube;
    else if (strcmp(s, "log_sin") == 0)
      ker = log_sin;
    else
    {
      mexErrMsgTxt("fastsum: Unknown kernel function.");
    }
	return ker;
}

static inline void zero_nodes_pointer(int i)
{
	// Initialize pointers that are set in init_nodes
	plans[i]->x = 0;
	plans[i]->y = 0;
}

/* cleanup on mex function unload */
static void cleanup(void)
{
  int i;

  if (!(gflags & FASTSUM_MEX_FIRST_CALL))
  {
    for (i = 0; i < PLANS_MAX; i++)
      if (plans[i])
      {
		nfft_free(plans[i]->kernel_param);
	    if(plans[i]->x)
		  fastsum_finalize_source_nodes(plans[i]);
	    if(plans[i]->y)
		  fastsum_finalize_target_nodes(plans[i]);
        fastsum_finalize_kernel(plans[i]);
        nfft_free(plans[i]);
        plans[i] = 0;
      }
    gflags |= FASTSUM_MEX_FIRST_CALL;
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char cmd[CMD_LEN_MAX];
  if (gflags & FASTSUM_MEX_FIRST_CALL)
  {
    /* Force Matlab to load libfftw3. There is at least one version of Matlab
     * which otherwise crashes upon invocation of this mex function. */
    mexEvalString("fft([1,2,3,4]);");

    nfft_mex_install_mem_hooks();

    /* plan pointers to zeros */
    {
      int i;
      for (i = 0; i < PLANS_MAX; i++)
        plans[i] = 0;
    }

    mexAtExit(cleanup);
    gflags &= ~FASTSUM_MEX_FIRST_CALL;
  }

  /* command string */
  DM(if (nrhs == 0)
    mexErrMsgTxt("At least one input required.");)

  DM(if (!mxIsChar(prhs[0]))
    mexErrMsgTxt("First argument must be a string.");)

  if (mxGetString(prhs[0], cmd, CMD_LEN_MAX))
    mexErrMsgTxt("Could not get command string.");

  if(strcmp(cmd,"get_num_threads") == 0)
  {
    int32_t nthreads = X(get_num_threads)();
    plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    *((int32_t *)mxGetData(plhs[0])) = nthreads;

    return;
  }
  
  else if (strcmp(cmd,"init") == 0)
  {
	check_nargs(nrhs,9,"Wrong number of arguments for init.");
	{
    int i;
	
	int d; /**< number of dimensions    */
	int n; /**< expansion degree        */
	int p; /**< degree of smoothness    */
	kernel ker; /**< kernel function         */
	double *param; /**< parameter for kernel    */
	double eps_I; /**< inner boundary          */
	double eps_B; /**< outer boundary          */
	param = nfft_malloc(sizeof(double));
	
    d = nfft_mex_get_int(prhs[1],"fastsum init: Input argument d must be a scalar.");
	DM(if (d < 1)
		mexErrMsgTxt("fastsum init: Input argument d must be positive.");)
	
	ker = get_kernel(prhs[2]);
	
    *param = nfft_mex_get_double(prhs[3],"fastsum init: Input argument c must be a scalar.");
	
	int flags_int = nfft_mex_get_int(prhs[4],"fastsum init: Input argument flags must be a scalar.");
	DM( if (flags_int < 0)
		mexErrMsgTxt("fastsum init: Input argument flags must be non-negative.");)
	unsigned flags = (unsigned) flags_int;
	
    n = nfft_mex_get_int(prhs[5],"fastsum init: Input argument n must be a scalar.");
	DM(if (n < 1)
		mexErrMsgTxt("fastsum init: Input argument n must be positive.");)
	
    p = nfft_mex_get_int(prhs[6],"fastsum init: Input argument p must be a scalar.");
	DM(if (p < 1)
		mexErrMsgTxt("fastsum init: Input argument p must be positive.");)
	
    eps_I = nfft_mex_get_double(prhs[7],"fastsum init: Input argument eps_I must be a scalar.");
    eps_B = nfft_mex_get_double(prhs[8],"fastsum init: Input argument eps_B must be a scalar.");

    i = mkplan();
    
	fastsum_init_guru_kernel(plans[i], d, ker, param, flags, n, p, eps_I, eps_B);
	
	zero_nodes_pointer(i);

    plhs[0] = mxCreateDoubleScalar((double)i);
	}
    return;
  }
  
  else if (strcmp(cmd,"set_x_alpha") == 0)
  {
    check_nargs(nrhs,6,"Wrong number of arguments for set_x_alpha.");
    {
	  int i = nfft_mex_get_int(prhs[1],"fastsum set_x_alpha: Input argument plan must be a scalar.");
	  check_plan(i);
	  
      int d = plans[i]->d;
      DM(if (!mxIsDouble(prhs[2]) || (mxGetNumberOfDimensions(prhs[2]) > 2) || (mxGetN(prhs[2]) != (unsigned)d))
        mexErrMsgTxt("Input argument x must be a matrix with d columns");
	  if (mxGetM(prhs[2]) > INT_MAX)
		mexErrMsgTxt("fastsum set_x_alpha: Input argument x is too large");)
		
      int N = (int) mxGetM(prhs[2]);
		
	  DM(if (!(mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) || (mxGetNumberOfDimensions(prhs[3]) > 2) || (mxGetN(prhs[3]) != 1))
        mexErrMsgTxt("Input argument alpha must be an N x 1 array");
	  if (mxGetM(prhs[3]) > INT_MAX)
		mexErrMsgTxt("fastsum set_x_alpha: Input argument alpha is too large");)
      int N_alpha = (int) mxGetM(prhs[3]);
	  
	  int nn_oversampled = nfft_mex_get_int(prhs[4],"fastsum set_x_alpha: Input argument nn_oversampled must be a scalar.");
	  DM(if (nn_oversampled < 1)
		mexErrMsgTxt("fastsum set_x_alpha: Input argument nn_oversampled must be positive.");)
	  
	  int m = nfft_mex_get_int(prhs[5],"fastsum set_x_alpha: Input argument m must be a scalar.");
	  DM(if (m < 1)
		mexErrMsgTxt("fastsum set_x_alpha: Input argument m must be positive.");)
	  
	  DM(if(N != N_alpha)
	    mexErrMsgTxt("fastsum set_x_alpha: Input arguments x and alpha must have the same number of rows");)
	  
	  if(!(plans[i]->x)) 
		fastsum_init_guru_source_nodes(plans[i], N, nn_oversampled, m);
	  else if( (N != plans[i]->N_total) || (m != plans[i]->mv1.m) || (nn_oversampled != plans[i]->mv1.n[0]) )
		{
		fastsum_finalize_source_nodes(plans[i]);
		fastsum_init_guru_source_nodes(plans[i], N, nn_oversampled, m);
		}
		
      {
        double *x = mxGetPr(prhs[2]);
        double *ar = mxGetPr(prhs[3]), *ai=0;
		if(mxIsComplex(prhs[3])) 
			ai = mxGetPi(prhs[3]);
		DM(double norm_max = (.25-(plans[i]->eps_B)*.5)*(.25-(plans[i]->eps_B)*.5);
		short warn=0;)
        for (int k = 0; k < N; k++)
        {
			DM(double norm = 0;)
			for (int t = 0; t < d; t++)
			{
			plans[i]->x[k*d+t] = x[k+t*N];
			if(ai)
				plans[i]->alpha[k] = ar[k] + _Complex_I*ai[k];
			else
				plans[i]->alpha[k] = ar[k];
			DM( if((plans[i]->x[k*d+t] < -0.5) || (plans[i]->x[k*d+t] >= 0.5))
				mexErrMsgTxt("x must be in interval [-0.5,0.5)");
			norm += plans[i]->x[k*d+t] * plans[i]->x[k*d+t];)
			}
			DM(if( norm > norm_max)
				warn = 1;)
        }
		DM(if(warn)
			mexWarnMsgTxt("x must be in ball with radius 1/4-eps_B/2.\nThis may cause wrong results or crashes!!");)
      }
	  
	  fastsum_precompute_source_nodes(plans[i]);
    }
    return;
  }
  
  else if (strcmp(cmd,"set_y") == 0)
  {
    check_nargs(nrhs,5,"Wrong number of arguments for set_y.");
    {
      int i = nfft_mex_get_int(prhs[1],"fastsum set_y: Input argument plan must be a scalar.");
	  check_plan(i);
	
      int d = plans[i]->d;
      DM(if (!mxIsDouble(prhs[2]) || mxGetNumberOfDimensions(prhs[2]) > 2 || mxGetN(prhs[2]) != (unsigned)d)
        mexErrMsgTxt("fastsum set_y: Input argument y must be an M x d double array");
	  if(mxGetM(prhs[2]) > INT_MAX)
	    mexErrMsgTxt("fastsum set_y: Input argument y is too large");)
      int M = (int) mxGetM(prhs[2]);
	  
	  int nn_oversampled = nfft_mex_get_int(prhs[3],"fastsum set_y: Input argument nn_oversampled must be a scalar.");
	  DM(if (nn_oversampled < 1)
		mexErrMsgTxt("fastsum set_y: Input argument nn_oversampled must be positive.");)
	  
	  int m = nfft_mex_get_int(prhs[4],"fastsum set_y: Input argument m must be a scalar.");
	  DM(if (m < 1)
		mexErrMsgTxt("fastsum set_y: Input argument m must be positive.");)
	  
	  if(!(plans[i]->y))
		fastsum_init_guru_target_nodes(plans[i], M, nn_oversampled, m);
	  else if( (M != plans[i]->M_total) || (m != plans[i]->mv2.m) || (nn_oversampled != plans[i]->mv2.n[0]) )
		{
		fastsum_finalize_target_nodes(plans[i]);
		fastsum_init_guru_target_nodes(plans[i], M, nn_oversampled, m);
		}
	  
      {
        double *y = mxGetPr(prhs[2]);
		DM(double norm_max = (.25-(plans[i]->eps_B)*.5)*(.25-(plans[i]->eps_B)*.5);
		short warn=0;)
        for (int j = 0; j < M; j++)
        {
			DM(double norm = 0;)
			for (int t = 0; t < d; t++)
			{
			plans[i]->y[d*j+t] = y[j+t*M];
			DM( if((plans[i]->y[d*j+t] < -0.5) || (plans[i]->y[d*j+t] >= 0.5))
				mexErrMsgTxt("y must be in interval [-0.5,0.5)");
			norm += plans[i]->y[d*j+t] * plans[i]->y[d*j+t];)
			}
			DM(if( norm > norm_max)
				warn = 1;)
        }
		DM(if(warn)
			mexWarnMsgTxt("y must be in ball with radius 1/4-eps_B/2.\nThis may cause wrong results or crashes!!");)
	  
	  fastsum_precompute_target_nodes(plans[i]);
      }
    }
    return;
  }
  
  else if (strcmp(cmd,"trafo_direct") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for trafo_direct.");
    {
      const int i = nfft_mex_get_int(prhs[1],"fastsum trafo_direct: Input argument plan must be a scalar.");
      check_plan(i);
	  check_plan_nodes(i);
      fastsum_exact(plans[i]);
    }
    return;
  }
  
  else if (strcmp(cmd,"trafo") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for trafo.");
    {
      const int i = nfft_mex_get_int(prhs[1],"fastsum trafo: Input argument plan must be a scalar.");
      check_plan(i);
	  check_plan_nodes(i);
      fastsum_trafo(plans[i]);
    }
    return;
  }
  
  else if (strcmp(cmd,"get_f") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for get_f.");
    {
      const int i = nfft_mex_get_int(prhs[1],"fastsum get_f: Input argument plan must be a scalar.");
      check_plan(i);
	  check_plan_nodes(i);
      const int M = plans[i]->M_total;
      plhs[0] = mxCreateDoubleMatrix((unsigned int)M, 1, mxCOMPLEX);
      {
        double *fr = mxGetPr(plhs[0]), *fi = mxGetPi(plhs[0]);
        for (int j = 0; j < M; j++)
        {
          fr[j] = creal(plans[i]->f[j]);
          fi[j] = cimag(plans[i]->f[j]);
        }
      }
    }
    return;
  }
  
  else if (strcmp(cmd,"finalize") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for finalize.");
    {
      const int i = nfft_mex_get_int(prhs[1],"fastsum finalize: Input argument plan must be a scalar.");
      check_plan(i);
	  nfft_free(plans[i]->kernel_param);
	  if(plans[i]->x || plans[i]->alpha)
		fastsum_finalize_source_nodes(plans[i]);
	  if(plans[i]->y)
		fastsum_finalize_target_nodes(plans[i]);
      fastsum_finalize_kernel(plans[i]);
      nfft_free(plans[i]);
      plans[i] = 0;
    }
    return;
  }
  
  
  
  
  // Auxiliary functions
  else if (strcmp(cmd,"init_guru") == 0)
  {
	check_nargs(nrhs,11,"Wrong number of arguments for init_guru.");
	{
    int i;
	
	int d; /**< number of dimensions    */
	int N; /**< number of source nodes  */
	int M; /**< number of target nodes  */
	int n; /**< expansion degree        */
	int m; /**< cut-off parameter       */
	int p; /**< degree of smoothness    */
	kernel ker; /**< kernel function         */
	double *param; /**< parameter for kernel    */
	double eps_I; /**< inner boundary          */
	double eps_B; /**< outer boundary          */
	param = nfft_malloc(sizeof(double));
	
    d = nfft_mex_get_int(prhs[1],"fastsum init_guru: Input argument d must be a scalar.");
	DM(if (d < 1)
		mexErrMsgTxt("fastsum init_guru: Input argument d must be positive.");)
	N = nfft_mex_get_int(prhs[2],"fastsum init_guru: Input argument N must be a scalar.");
	DM(if (N < 1)
		mexErrMsgTxt("fastsum init_guru: Input argument N must be positive.");)
	M = nfft_mex_get_int(prhs[3],"fastsum init_guru: Input argument M must be a scalar.");
	DM(if (M < 1)
		mexErrMsgTxt("fastsum init_guru: Input argument M must be positive.");)
    n = nfft_mex_get_int(prhs[4],"fastsum init_guru: Input argument n must be a scalar.");
	DM(if (n < 1)
		mexErrMsgTxt("fastsum init_guru: Input argument n must be positive.");)
    m = nfft_mex_get_int(prhs[5],"fastsum init_guru: Input argument m must be a scalar.");
	DM(if (m < 1)
		mexErrMsgTxt("fastsum init_guru: Input argument m must be positive.");)
    p = nfft_mex_get_int(prhs[6],"fastsum init_guru: Input argument p must be a scalar.");
	DM(if (p < 1)
		mexErrMsgTxt("fastsum init_guru: Input argument p must be positive.");)
	ker = get_kernel(prhs[7]);
    *param = nfft_mex_get_double(prhs[8],"fastsum init_guru: Input argument c must be a scalar.");
    eps_I = nfft_mex_get_double(prhs[9],"fastsum init_guru: Input argument eps_I must be a scalar.");
    eps_B = nfft_mex_get_double(prhs[10],"fastsum init_guru: Input argument eps_B must be a scalar.");

    i = mkplan();
    
	fastsum_init_guru(plans[i], d, N, M, ker, param, 0, n, m, p, eps_I, eps_B);

    plhs[0] = mxCreateDoubleScalar((double)i);
	}
    return;
  }
  
  else if (strcmp(cmd,"get_x") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for get_x.");
    {
      const int i = nfft_mex_get_int(prhs[1],"fastsum: Input argument plan must be a scalar.");
      check_plan(i);
	  check_plan_source_nodes(i);
      const int d = plans[i]->d;
      const int N = plans[i]->N_total;
      plhs[0] = mxCreateDoubleMatrix((unsigned int)N, (unsigned int)d, mxREAL);
      {
		double *x = mxGetPr(plhs[0]);
        for (int j = 0; j < N; j++)
        {
			for (int t = 0; t < d; t++)
			{
			x[j+t*N] = plans[i]->x[d*j+t];
			}
        }
      }
    }
    return;
  }
  
  else if (strcmp(cmd,"get_alpha") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for get_alpha.");
    {
      const int i = nfft_mex_get_int(prhs[1],"fastsum: Input argument plan must be a scalar.");
      check_plan(i);
	  check_plan_source_nodes(i);
      const int N = plans[i]->N_total;
      plhs[0] = mxCreateDoubleMatrix((unsigned int)N, 1, mxCOMPLEX);
      {
		double *ar = mxGetPr(plhs[0]), *ai = mxGetPi(plhs[0]);
        for (int j = 0; j < N; j++)
        {
			ar[j] = creal(plans[i]->alpha[j]);
			ai[j] = cimag(plans[i]->alpha[j]);
        }
      }
    }
    return;
  }
  
  else if (strcmp(cmd,"get_y") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for get_y.");
    {
      const int i = nfft_mex_get_int(prhs[1],"fastsum: Input argument plan must be a scalar.");
      check_plan(i);
	  check_plan_target_nodes(i);
      const int d = plans[i]->d;
      const int M = plans[i]->M_total;
      plhs[0] = mxCreateDoubleMatrix((unsigned int)M, (unsigned int)d, mxREAL);
      {
		double *y = mxGetPr(plhs[0]);
        for (int j = 0; j < M; j++)
        {
			for (int t = 0; t < d; t++)
			{
			y[j+t*M] = plans[i]->y[d*j+t];
			}
        }
      }
    }
    return;
  }
  
  else if (strcmp(cmd,"get_b") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for get_b.");
    {
      const int i = nfft_mex_get_int(prhs[1],"fastsum: Input argument plan must be a scalar.");
      check_plan(i);
      const int d = plans[i]->d;
	  size_t dims[d];
	  int n_total = 1;
	  for(int j=0; j<d; j++)
	  {
		  dims[j] = (unsigned int) plans[i]->n;
		  n_total *= plans[i]->n;
	  }
      plhs[0] = mxCreateNumericArray((unsigned int)d, dims, mxDOUBLE_CLASS, mxCOMPLEX);
      {
		double *br = mxGetPr(plhs[0]), *bi = mxGetPi(plhs[0]);
        for (int j = 0; j < n_total; j++)
        {
			br[j] = creal(plans[i]->b[j]);
			bi[j] = cimag(plans[i]->b[j]);
        }
      }
    }
    return;
  }
  
  else if (strcmp(cmd,"display") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for display.");
    {
      const int i = nfft_mex_get_int(prhs[1],"fastsum: Input argument plan must be a scalar.");
	  check_plan(i);
      mexPrintf("Plan %d\n",i);
      mexPrintf("  pointer: %p\n",plans[i]);
      mexPrintf("        d: %d\n",plans[i]->d);
      mexPrintf("  N_total: %d\n",plans[i]->N_total);
      mexPrintf("  M_total: %d\n",plans[i]->M_total);
      mexPrintf("        n: %d\n",plans[i]->n);
      mexPrintf("        p: %d\n",plans[i]->p);
      mexPrintf("    eps_I: %f\n",plans[i]->eps_I);
      mexPrintf("    eps_B: %f\n",plans[i]->eps_B);
      mexPrintf("        x: %p\n",plans[i]->x);
      mexPrintf("        y: %p\n",plans[i]->y);
      mexPrintf("    alpha: %p\n",plans[i]->alpha);
      mexPrintf("        f: %p\n",plans[i]->f);
      mexPrintf("   kernel: %p\n",plans[i]->k);
      mexPrintf("   _param: %p\n",plans[i]->kernel_param);
      mexPrintf("  *_param: %f\n",*(plans[i]->kernel_param));
      mexPrintf("    flags: %d\n",plans[i]->flags);
    }
    return;
  }
  
  else
    mexErrMsgTxt("fastsum: Unknown command.\n");
}
