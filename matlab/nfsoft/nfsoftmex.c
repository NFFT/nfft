/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
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

#define PLANS_START 10 /* initial number of plans */
#define CMD_LEN_MAX 20 /* maximum length of command argument */

/* global flags */
#define NFSOFT_MEX_FIRST_CALL (1U << 0)
static unsigned short gflags = NFSOFT_MEX_FIRST_CALL;

static nfsoft_plan** plans = NULL; /* plans */
static unsigned int plans_num_allocated = 0;
static int n_max = -1; /* maximum degree precomputed */
static char cmd[CMD_LEN_MAX];

static inline int get_plan(const mxArray *pm)
{
  int i = nfft_mex_get_int(pm,"Input argument plan must be a scalar.");
  DM(if (i < 0 || i >= plans_num_allocated || plans[i] == 0)
    mexErrMsgTxt("Plan was not initialized or has already been finalized");)
  return i;
}

static inline int mkplan()
{
  int i = 0;
  while (i < plans_num_allocated && plans[i] != 0) i++;
  if (i == plans_num_allocated)
  {
    int l;

    if (plans_num_allocated >= INT_MAX - PLANS_START - 1)
      mexErrMsgTxt("nfsoft: Too many plans already allocated.");

    nfsoft_plan** plans_old = plans;
    plans = nfft_malloc((plans_num_allocated+PLANS_START)*sizeof(nfsoft_plan*));
    for (l = 0; l < plans_num_allocated; l++)
      plans[l] = plans_old[l];
    for (l = plans_num_allocated; l < plans_num_allocated+PLANS_START; l++)
      plans[l] = 0;
    if (plans_num_allocated > 0)
      nfft_free(plans_old);
    plans_num_allocated += PLANS_START;
  }
  plans[i] = nfft_malloc(sizeof(nfsoft_plan));
  return i;
}

/* cleanup on mex function unload */
static void cleanup(void)
{
  int i;

  if (!(gflags & NFSOFT_MEX_FIRST_CALL))
  {
    for (i = 0; i < plans_num_allocated; i++)
      if (plans[i])
      {
        nfsoft_finalize(plans[i]);
        nfft_free(plans[i]);
        plans[i] = 0;
      }

    if (plans_num_allocated > 0)
    {
      nfft_free(plans);
      plans = NULL;
      plans_num_allocated = 0;
    }
    gflags |= NFSOFT_MEX_FIRST_CALL;
    n_max = -1;
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (gflags & NFSOFT_MEX_FIRST_CALL)
  {
    /* Force Matlab to load libfftw3. There is at least one version of Matlab
     * which otherwise crashes upon invocation of this mex function. */
    mexEvalString("fft([1,2,3,4]);");

    nfft_mex_install_mem_hooks();


    mexAtExit(cleanup);
    gflags &= ~NFSOFT_MEX_FIRST_CALL;
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
    INT nthreads = X(get_num_threads)();
    plhs[0] = mxCreateDoubleScalar((double) nthreads);
    return;
  }
  else if(strcmp(cmd,"set_num_threads") == 0)
  {
    int nthreads_new = nfft_mex_set_num_threads_check(nrhs, prhs, (void **) plans, plans_num_allocated);
    X(set_num_threads)(nthreads_new);

    return;
  }
  else if(strcmp(cmd,"has_threads_enabled") == 0)
  {
    INT threads_enabled = X(has_threads_enabled)();
    plhs[0] = mxCreateDoubleScalar((double) threads_enabled);
    return;
  }
  else if(strcmp(cmd,"init") == 0)
  {
    nfft_mex_check_nargs(nrhs,8,"Wrong number of arguments for init.");
	int N = nfft_mex_get_int(prhs[1],"N must be scalar");
 	DM( if (N <= 0)
		mexErrMsgTxt("Input argument N must be positive.");)
	int M = nfft_mex_get_int(prhs[2],"M must be scalar");
	DM( if (M <= 0)
		mexErrMsgTxt("Input argument M must be positive.");)
	int flags_int = nfft_mex_get_int(prhs[3],"Input argument flags must be a scalar.");
	DM( if (flags_int < 0)
		mexErrMsgTxt("Input argument flags must be non-negative.");)
	unsigned flags = (unsigned) flags_int;
	int nfft_flags_int = nfft_mex_get_int(prhs[4],"Input argument nfft_flags must be a scalar.");
	DM( if (nfft_flags_int < 0)
		mexErrMsgTxt("Input argument nfft_flags must be non-negative.");
	  if (nfft_flags_int & (MALLOC_X | MALLOC_F))
	    mexErrMsgTxt("nfft_flags must not contain MALLOC_X or MALLOC_F.");)
	unsigned nfft_flags = (unsigned) nfft_flags_int;
	int nfft_cutoff = nfft_mex_get_int(prhs[5],"Input argument nfft_cutoff must be a scalar.");
	DM( if (nfft_cutoff <= 0)
		mexErrMsgTxt("Input argument nfft_cutoff must be positive.");)
	int fpt_kappa = nfft_mex_get_int(prhs[6],"Input argument fpt_kappa must be a scalar.");
	DM( if (fpt_kappa <= 0)
		mexErrMsgTxt("Input argument fpt_kappa must be positive.");)
	int fftw_size = nfft_mex_get_int(prhs[7],"Input argument nn_oversampled must be a scalar.");
	DM( if ((fftw_size %2) || (fftw_size < 2*N+2))
		mexErrMsgTxt("Input argument fftw_size must be even and at least 2*N+2.");)
	
	int i = mkplan();
	nfsoft_init_guru_advanced(plans[i], N, M, flags | NFSOFT_MALLOC_X | NFSOFT_MALLOC_F | NFSOFT_MALLOC_F_HAT,
	nfft_flags | PRE_PHI_HUT | PRE_PSI | MALLOC_F_HAT | FFTW_INIT, nfft_cutoff, fpt_kappa, fftw_size);
    plans[i]->p_nfft.f = plans[i]->f;
    plans[i]->p_nfft.x = plans[i]->x;
    plhs[0] = mxCreateDoubleScalar((double)i);
    return;
  }

  else if(strcmp(cmd,"set_x") == 0)
  {
    nfft_mex_check_nargs(nrhs,3,"Wrong number of arguments for set_x.");
    int i = get_plan(prhs[1]);

    DM(if (!(mxIsDouble(prhs[2])) || (mxGetNumberOfDimensions(prhs[2]) > 2) || (mxGetN(prhs[2]) != plans[i]->M_total) || mxGetM(prhs[2]) != 3))
    mexErrMsgTxt("Input argument x must be a real 3 x M array");
    double *x = mxGetPr(prhs[2]);

    for (int j = 0; j < plans[i]->M_total; j++)
    {
    plans[i]->p_nfft.x[3*j]   = x[3*j+2] / K2PI;  // gamma
    plans[i]->p_nfft.x[3*j+1] = x[3*j]   / K2PI;  // alpha
    plans[i]->p_nfft.x[3*j+2] = x[3*j+1] / K2PI;  // beta
    }
    nfsoft_precompute(plans[i]);
    return;
  }

  else if(strcmp(cmd,"set_f_hat") == 0)
  {
    nfft_mex_check_nargs(nrhs,3,"Wrong number of arguments for set_f_hat.");
	int i = get_plan(prhs[1]);
	
	int N = plans[i]->N_total;
	int fh_size = NFSOFT_F_HAT_SIZE(N);
	
	DM(if (!(mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) || (mxGetNumberOfDimensions(prhs[2]) > 2) || (mxGetN(prhs[2]) != 1) || mxGetM(prhs[2]) != fh_size))
        mexErrMsgTxt("Input argument f_hat must be a NFSOFT_F_HAT_SIZE(N) x 1 array");
	double *fh_real = mxGetPr(prhs[2]), *fh_imag = mxGetPi(prhs[2]);
	int glo1 = 0;
	for (int k = -N; k <= N; k++)
	{
	for (int m = -N; m <= N; m++)
	{
		int max = (ABS(m) > ABS(k) ? ABS(m) : ABS(k));
		for (int j = max; j <= N; j++)		// j = polynomial degree
		{
		int my_ind = NFSOFT_F_HAT_SIZE(j-1) + (m+j)*(2*j+1) + (k+j);
		if(fh_imag)
			plans[i]->f_hat[glo1] = fh_real[my_ind] + I*fh_imag[my_ind];
		else 
			plans[i]->f_hat[glo1] = fh_real[my_ind];
		glo1++;
		}
	}
	}
    return;
  }

  else if(strcmp(cmd,"set_f") == 0)
  {
    nfft_mex_check_nargs(nrhs,3,"Wrong number of arguments for set_f.");
	int i = get_plan(prhs[1]);
	
	DM(if (!(mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) || (mxGetNumberOfDimensions(prhs[2]) > 2) || (mxGetN(prhs[2]) != 1) || mxGetM(prhs[2]) != plans[i]->M_total))
        mexErrMsgTxt("Input argument f must be an M x 1 array");
	double *f_real = mxGetPr(prhs[2]), *f_imag=0;
	if(mxIsComplex(prhs[2])) 
		f_imag = mxGetPi(prhs[2]);
	for (int j = 0; j < plans[i]->M_total;j++)
		{
		if(f_imag)
			plans[i]->f[j] = f_real[j] + I*f_imag[j];
		else plans[i]->f[j] = f_real[j];
		}
    return;
  }

  else if(strcmp(cmd,"precompute") == 0)
  {
      /* Do nothing here.
       * nfsft_precompute_x has been moved to the set_x routine. */
      return;
    }

  else if(strcmp(cmd,"trafo") == 0)
  {
    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for trafo.");
	int i = get_plan(prhs[1]);
	nfsoft_trafo(plans[i]);
    return;
  }

  else if(strcmp(cmd,"adjoint") == 0)
  {
    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for adjoint.");
	int i = get_plan(prhs[1]);
	nfsoft_adjoint(plans[i]);
    return;
  }

  else if(strcmp(cmd,"get_x") == 0)
  {
    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for get_x.");
  int i = get_plan(prhs[1]);

  plhs[0] = mxCreateDoubleMatrix(3, (unsigned int) plans[i]->M_total, mxREAL);
  double *x = mxGetPr(plhs[0]);

  for (int j = 0; j < plans[i]->M_total; j++)
    {
    x[3*j+2] = plans[i]->p_nfft.x[3*j]   * K2PI;  // gamma
    x[3*j]   = plans[i]->p_nfft.x[3*j+1] * K2PI;  // alpha
    x[3*j+1] = plans[i]->p_nfft.x[3*j+2] * K2PI;  // beta
    }
  return;
  }

  else if(strcmp(cmd,"get_f") == 0)
  {
    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for get_f.");
	int i = get_plan(prhs[1]);
	
      plhs[0] = mxCreateDoubleMatrix((unsigned int)plans[i]->M_total, 1, mxCOMPLEX);
      {
        double *fr = mxGetPr(plhs[0]), *fi = mxGetPi(plhs[0]);
        for (int j = 0; j < plans[i]->M_total; j++)
        {
          fr[j] = creal(plans[i]->f[j]);
          fi[j] = cimag(plans[i]->f[j]);
        }
      }
    return;
  }

  else if(strcmp(cmd,"get_f_hat") == 0)
  {
    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for get_f_hat.");
	int i = get_plan(prhs[1]);
	int N = plans[i]->N_total;
	int fh_size = NFSOFT_F_HAT_SIZE(N);
	
    plhs[0] = mxCreateDoubleMatrix((unsigned int)fh_size, 1, mxCOMPLEX);
    double *fr = mxGetPr(plhs[0]), *fi = mxGetPi(plhs[0]);
	
	int glo1 = 0;
	for (int k = -N; k <= N; k++)
	{
	  for (int m = -N; m <= N; m++)
	  {
		int max = (ABS(m) > ABS(k) ? ABS(m) : ABS(k));
		for (int j = max; j <= N; j++)		// j = polynomial degree
		{
		int my_ind = NFSOFT_F_HAT_SIZE(j-1) + (m+j)*(2*j+1) + (k+j);
		fr[my_ind] = creal(plans[i]->f_hat[glo1]);
		fi[my_ind] = cimag(plans[i]->f_hat[glo1]);
		glo1++;
		}
	  }
	}
    return;
  }

  else if(strcmp(cmd,"finalize") == 0)
  {
    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments finalize.");
	int i = get_plan(prhs[1]);
	
    nfsoft_finalize(plans[i]);
    nfft_free(plans[i]);
    plans[i] = 0;
    return;
  }
  
  else if (strcmp(cmd,"display") == 0)
  {
    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for display.");
    {
	  int i = get_plan(prhs[1]);
      mexPrintf("Plan %d\n",i);
      mexPrintf("  pointer: %p\n",plans[i]);
      mexPrintf("  N_total: %d\n",plans[i]->N_total);
      mexPrintf("  M_total: %d\n",plans[i]->M_total);
      mexPrintf("        x: %p\n",plans[i]->x);
      mexPrintf("        f: %p\n",plans[i]->f);
      mexPrintf("    f_hat: %p\n",plans[i]->f_hat);
      mexPrintf("    flags: %d\n",plans[i]->flags);
    }
    return;
  }
  
  else
    mexErrMsgTxt("nfsoft: Unknown command.\n");
}
