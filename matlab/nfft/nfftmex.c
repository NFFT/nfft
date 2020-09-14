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
#define CMD_LEN_MAX 40 /* maximum length of command argument */

/* global flags */
#define NFFT_MEX_FIRST_CALL (1U << 0)
static unsigned short gflags = NFFT_MEX_FIRST_CALL;

static X(plan)** plans = NULL; /* plans */
static unsigned int plans_num_allocated = 0;
static char cmd[CMD_LEN_MAX];

static inline void get_guru(const mxArray *prhs[], int d, int *N, int *M, int *n, int *m, unsigned int *f1, unsigned int *f2)
{
  DM(if (mxGetNumberOfElements(prhs[1]) != 2*d+5)
    mexErrMsgTxt("init_guru: Cell array must have 2*d+5 many entries: {d, N(1),...,N(d), M, n(1),...,n(d), m, nfft_flag, fftw_flags}");)

  int k;

  for(k=0;k<d;k++)
  {
    N[k] = mxGetScalar(mxGetCell(prhs[1], (unsigned int)(1+k)));
    DM(if (N[k] < 2 || (N[k]%2 != 0))
      mexErrMsgTxt("init_guru: input vector N must consist of even natural numbers");)
  }

  *M = mxGetScalar(mxGetCell(prhs[1], (unsigned int)(d+1)));
  DM(if (*M < 1)
    mexErrMsgTxt("init_guru: input argument M must be a natural number");)

  for(k=0;k<d;k++)
    n[k] = mxGetScalar(mxGetCell(prhs[1], (unsigned int)(d+2+k)));

  *m = mxGetScalar(mxGetCell(prhs[1], (unsigned int)(2*d+2)));
  DM(if (*m < 1)
    mexErrMsgTxt("init_guru: input argument m must be a natural number");)

  *f1 = mxGetScalar(mxGetCell(prhs[1], (unsigned int)(2*d+3)));

  *f2 = mxGetScalar(mxGetCell(prhs[1], (unsigned int)(2*d+4)));
}

static inline void check_plan(int i)
{
  DM(if (i < 0 || i >= plans_num_allocated || plans[i] == 0)
    mexErrMsgTxt("Plan was not initialized or has already been finalized");)
}

static inline int mkplan(void)
{
  int i = 0;
  while (i < plans_num_allocated && plans[i] != 0) i++;
  if (i == plans_num_allocated)
  {
    int l;

    if (plans_num_allocated >= INT_MAX - PLANS_START - 1)
      mexErrMsgTxt("nfft: Too many plans already allocated.");

    X(plan)** plans_old = plans;
    plans = X(malloc)((plans_num_allocated+PLANS_START)*sizeof(nfft_plan*));
    for (l = 0; l < plans_num_allocated; l++)
      plans[l] = plans_old[l];
    for (l = plans_num_allocated; l < plans_num_allocated+PLANS_START; l++)
      plans[l] = 0;
    if (plans_num_allocated > 0)
      X(free)(plans_old);
    plans_num_allocated += PLANS_START;
  }
  plans[i] = X(malloc)(sizeof(nfft_plan));
  return i;
}

/* cleanup on mex function unload */
static void cleanup(void)
{
  int i;

  if (!(gflags & NFFT_MEX_FIRST_CALL))
  {
    for (i = 0; i < plans_num_allocated; i++)
      if (plans[i])
      {
        X(finalize)(plans[i]);
	X(free)(plans[i]);
	plans[i] = 0;
      }

    if (plans_num_allocated > 0)
    {
      X(free)(plans);
      plans = NULL;
      plans_num_allocated = 0;
    }
    gflags |= NFFT_MEX_FIRST_CALL;
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (gflags & NFFT_MEX_FIRST_CALL)
  {
    /* Force Matlab to load libfftw3. There is at least one version of Matlab
     * which otherwise crashes upon invocation of this mex function. */
    mexEvalString("fft([1,2,3,4]);");

    nfft_mex_install_mem_hooks();

    mexAtExit(cleanup);
    gflags &= ~NFFT_MEX_FIRST_CALL;
  }

  /* command string */
  DM(if (nrhs == 0)
    mexErrMsgTxt("At least one input required.");)

  DM(if (!mxIsChar(prhs[0]))
    mexErrMsgTxt("First argument must be a string.");)

  if (mxGetString(prhs[0], cmd, CMD_LEN_MAX))
    mexErrMsgTxt("Could not get command string.");

  if (strcmp(cmd,"init_1d") == 0)
  {
    nfft_mex_check_nargs(nrhs,3,"Wrong number of arguments for init_1d.");
    {
      int i;
      int N, M;
      nfft_mex_get_nm(prhs,&N,&M);
      i = mkplan();
      X(init_1d)(plans[i],N,M);
      plhs[0] = mxCreateDoubleScalar((double)i);
    }
    return;
  }
  else if (strcmp(cmd,"init_2d") == 0)
  {
    nfft_mex_check_nargs(nrhs,4,"Wrong number of arguments for init_2d.");
    {
      int i;
      int N1, N2, M;
      nfft_mex_get_n1n2m(prhs,&N1,&N2,&M);
      i = mkplan();
      X(init_2d)(plans[i],N1,N2,M);
      plhs[0] = mxCreateDoubleScalar((double)i);
    }
    return;
  }
  else if (strcmp(cmd,"init_3d") == 0)
  {
    nfft_mex_check_nargs(nrhs,5,"Wrong number of arguments for init_3d.");
    {
      int i;
      int N1, N2, N3, M;
      nfft_mex_get_n1n2n3m(prhs,&N1,&N2,&N3,&M);
      i = mkplan();
      X(init_3d)(plans[i],N1,N2,N3,M);
      plhs[0] = mxCreateDoubleScalar((double)i);
    }
    return;
  }
  else if (strcmp(cmd,"init") == 0)
  {
    nfft_mex_check_nargs(nrhs,3,"Wrong number of arguments for init.");
    {
      int i;
      int d, M;
      
      DM(if (!mxIsDouble(prhs[1]) || mxGetNumberOfDimensions(prhs[1]) > 2)
        mexErrMsgTxt("init: input argument N must be a double array");)

      d = mxGetNumberOfElements(prhs[1]);
      if (d < 1)
        mexErrMsgTxt("init: input argument N must be a double array of length >= 1");
      else
      {
        int t, N[d];
        double *N_vec = mxGetPr(prhs[1]);
        for (t = 0; t < d; t++)
        {
          N[t] = (int) N_vec[t];
          DM(if (N[t] < 2 || (N[t]%2 != 0))
            mexErrMsgTxt("init: input vector N must consist of even natural numbers");)
        }
        M = nfft_mex_get_int(prhs[2],"init: input argument M must be a scalar.");
        DM(if (M < 1)
          mexErrMsgTxt("init: input argument NM must be a natural number");)
        i = mkplan();
        X(init)(plans[i],d,N,M);
        plhs[0] = mxCreateDoubleScalar((double)i);
      }
    }
    return;
  }
  else if (strcmp(cmd,"init_guru") == 0)
  {
    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for init_guru");

    DM(if (!mxIsCell(prhs[1]) || mxGetNumberOfElements(prhs[1]) < 3)
      mexErrMsgTxt("init_guru: second argument must be a cell array of length 2*d+5");)

    DM(if (!mxIsScalar(mxGetCell(prhs[1], 0)))
      mexErrMsgTxt("init_guru: cell array entry at index 0 (dimension d) must be scalar");)

    const int d = mxGetScalar(mxGetCell(prhs[1], 0));

    DM(if (d < 1)
      mexErrMsgTxt("init_guru: cell array entry at index 0 (dimension d) must be a natural number");)
    
    int N[d],n[d],m,M;
    unsigned int f1,f2;

    get_guru(prhs,d,N,&M,n,&m,&f1,&f2);
    int i = mkplan();
    X(init_guru)(plans[i],d,N,M,n,m,
		   f1 | MALLOC_X | MALLOC_F | MALLOC_F_HAT | FFTW_INIT,
		   f2);

    plhs[0] = mxCreateDoubleScalar((double)i);

    return;
  }
  else if (strcmp(cmd,"precompute_psi") == 0)
  {
    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for precompute_one_psi.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      check_plan(i);
      X(precompute_one_psi)(plans[i]);
    }
    return;
  }
  else if (strcmp(cmd,"trafo") == 0)
  {
    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for trafo.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      check_plan(i);
      X(trafo)(plans[i]);
    }
    return;
  }
  else if (strcmp(cmd,"adjoint") == 0)
  {
    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for adjoint.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      check_plan(i);
      X(adjoint)(plans[i]);
    }
    return;
  }
  else if (strcmp(cmd,"finalize") == 0)
  {
    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for finalize.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      check_plan(i);
      X(finalize)(plans[i]);
      X(free)(plans[i]);
      plans[i] = 0;
    }
    return;
  }
  else if (strcmp(cmd,"trafo_direct") == 0)
  {
    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for trafo direct.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      check_plan(i);
      X(trafo_direct)(plans[i]);
    }
    return;
  }
  else if (strcmp(cmd,"adjoint_direct") == 0)
  {
    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for adjoint direct.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      check_plan(i);
      X(adjoint_direct)(plans[i]);
    }
    return;
  }
  else if (strcmp(cmd,"solver") == 0)
  {
    nfft_mex_check_nargs(nrhs,4,"Wrong number of arguments for solver.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      check_plan(i);
      unsigned iterations = nfft_mex_get_int(prhs[2],"nfft: Input argument iterations must be a scalar.");
      unsigned flags = nfft_mex_get_int(prhs[3],"nfft: Input argument flags must be a scalar.");
      if (flags == 0U)
        flags = CGNR;

      SOLVER(plan_complex) ip; /**< plan for the inverse nfft       */
      SOLVER(init_advanced_complex)(&ip, (NFFT(mv_plan_complex)*) plans[i], flags);

      for (int k=0; k < plans[i]->M_total; k++)
        ip.y[k] = plans[i]->f[k];

      for (int k = 0; k < plans[i]->N_total; k++) // copy initial guess
        ip.f_hat_iter[k] = plans[i]->f_hat[k];

      SOLVER(before_loop_complex)(&ip);
      for (int l=0; l < iterations; l++)
        SOLVER(loop_one_step_complex)(&ip);
      for (int k=0; k < plans[i]->N_total; k++)
        plans[i]->f_hat[k] = ip.f_hat_iter[k];
      SOLVER(finalize_complex)(&ip);
    }
    return;
  }
  else if (strcmp(cmd,"get_x") == 0)
  {
    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for get_x.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      int m, d;
      check_plan(i);
      m = plans[i]->M_total;
      d = plans[i]->d;
      const mwSize dims[2] = {d, m};
#if defined(NFFT_SINGLE)
      plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
      float *x = (float*)mxGetData(plhs[0]);
#else
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      double *x = (double*)mxGetData(plhs[0]);
#endif
      {
        int j,t;
        for (j = 0; j < m; j++)
          for (t = 0; t < d; t++)
            x[d*j+t] = plans[i]->x[d*j+t];
      }
    }
    return;
  }
  else if (strcmp(cmd,"get_f") == 0)
  {
    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for get_f.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      int m;
      check_plan(i);
      m = plans[i]->M_total;
      const mwSize dims[2] = {m, 1};
#if defined(NFFT_SINGLE)
      plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxCOMPLEX);
      float *fr = (float*)mxGetData(plhs[0]), *fi = (float*)mxGetImagData(plhs[0]);
#else
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);
      double *fr = (double*)mxGetData(plhs[0]), *fi = (double*)mxGetImagData(plhs[0]);
#endif
      {
        int j;
        for (j = 0; j < m; j++)
        {
          fr[j] = CREAL(plans[i]->f[j]);
          fi[j] = CIMAG(plans[i]->f[j]);
        }
      }
    }
    return;
  }
  else if (strcmp(cmd,"get_f_hat") == 0)
  {
    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for get_f_hat.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      int n;
      check_plan(i);
      n = plans[i]->N_total;
      const mwSize dims[2] = {n, 1};
#if defined(NFFT_SINGLE)
      plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxCOMPLEX);
      float *f_hatr = (float*)mxGetData(plhs[0]), *f_hati = (float*)mxGetImagData(plhs[0]);
#else
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);
      double *f_hatr = (double*)mxGetData(plhs[0]), *f_hati = (double*)mxGetImagData(plhs[0]);
#endif
      {
        int k;
        for (k = 0; k < n; k++)
        {
          f_hatr[k] = CREAL(plans[i]->f_hat[k]);
          f_hati[k] = CIMAG(plans[i]->f_hat[k]);
        }
      }
    }
    return;
  }
  else if (strcmp(cmd,"set_x") == 0)
  {
    nfft_mex_check_nargs(nrhs,3,"Wrong number of arguments for set_x.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      int m, d;
      check_plan(i);
      m = plans[i]->M_total;
      d = plans[i]->d;
      DM(if (mxGetNumberOfDimensions(prhs[2]) > 2)
        mexErrMsgTxt("Input argument x must be a d x M double array");)
      DM(if (mxGetM(prhs[2]) != (unsigned int)d || mxGetN(prhs[2]) != (unsigned int)m)
        mexErrMsgTxt("Input argument x must have correct size (d x M).");)
      if (mxIsDouble(prhs[2]))
      {
        double *x = mxGetPr(prhs[2]);
        for (int j = 0; j < m; j++)
          for (int t = 0; t < d; t++)
            plans[i]->x[d*j+t] = x[d*j+t];
      }
      else if (mxIsSingle(prhs[2]))
      {
        float *x = (float*)mxGetData(prhs[2]);
        for (int j = 0; j < m; j++)
          for (int t = 0; t < d; t++)
            plans[i]->x[d*j+t] = x[d*j+t];
      }
      else
        DM(mexErrMsgTxt("Input argument x must be a d x M double array");)
    }
    return;
  }
  else if (strcmp(cmd,"set_f") == 0)
  {
    nfft_mex_check_nargs(nrhs,3,"Wrong number of arguments for set_f.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      int m;
      check_plan(i);
      m = plans[i]->M_total;
      DM(if (mxGetM(prhs[2]) * mxGetN(prhs[2]) != (unsigned int)m)
        mexErrMsgTxt("Input argument f must have correct size.");)
      if (mxIsDouble(prhs[2]))
      {
        double *fr = mxGetPr(prhs[2]), *fi = mxGetPi(prhs[2]);
        if (fi)
          for (int j = 0; j < m; j++)
            plans[i]->f[j] = fr[j] + I*fi[j];
        else
          for (int j = 0; j < m; j++)
            plans[i]->f[j] = fr[j];
      }
      else if (mxIsSingle(prhs[2]))
      {
        float *fr = (float*)mxGetData(prhs[2]), *fi = (float*)mxGetImagData(prhs[2]);
        if (fi)
          for (int j = 0; j < m; j++)
            plans[i]->f[j] = fr[j] + I*fi[j];
        else
          for (int j = 0; j < m; j++)
            plans[i]->f[j] = fr[j];
      }
      else
        DM(mexErrMsgTxt("Input argument f must be a double array");)
    }
    return;
  }
  else if (strcmp(cmd,"set_f_hat") == 0)
  {
    nfft_mex_check_nargs(nrhs,3,"Wrong number of arguments for set_f_hat.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      int n;
      check_plan(i);
      n = plans[i]->N_total;
      DM(if (   mxGetM(prhs[2]) != (unsigned int)n || mxGetN(prhs[2]) != 1)
        mexErrMsgTxt("Input argument f_hat must have correct size.");)
      if (mxIsDouble(prhs[2]))
      {
        double *f_hatr = mxGetPr(prhs[2]), *f_hati = mxGetPi(prhs[2]);
        if (f_hati)
          for (int k = 0; k < n; k++)
            plans[i]->f_hat[k] = f_hatr[k] + I*f_hati[k];
        else
          for (int k = 0; k < n; k++)
            plans[i]->f_hat[k] = f_hatr[k];
      }
      else if (mxIsSingle(prhs[2]))
      {
        float *f_hatr = (float*)mxGetData(prhs[2]), *f_hati = (float*)mxGetImagData(prhs[2]);
        if (f_hati)
          for (int k = 0; k < n; k++)
            plans[i]->f_hat[k] = f_hatr[k] + I*f_hati[k];
        else
          for (int k = 0; k < n; k++)
            plans[i]->f_hat[k] = f_hatr[k];
      }
      else
        DM(mexErrMsgTxt("Input argument f_hat must be a double array");)
    }
    return;
  }
  else if (strcmp(cmd,"display") == 0)
  {
    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for set_f_hat_linear.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      mexPrintf("Plan %d\n",i);
      check_plan(i);
      mexPrintf("   pointer: %p\n",plans[i]);
      mexPrintf("         d: %d\n",plans[i]->d);
      mexPrintf("   N_total: %d\n",plans[i]->N_total);
      mexPrintf("   M_total: %d\n",plans[i]->M_total);
      mexPrintf("      n[1]: %d\n",plans[i]->n[1]);
      mexPrintf("         m: %d\n",plans[i]->m);
      mexPrintf("         x: %p\n",plans[i]->x);
      mexPrintf("         f: %p\n",plans[i]->f);
      mexPrintf("     f_hat: %p\n",plans[i]->f_hat);
      mexPrintf("     flags: %d\n",plans[i]->flags);
      mexPrintf("fftw_flags: %d\n",plans[i]->fftw_flags);
    }
    return;
  }
  else if(strcmp(cmd,"get_num_threads") == 0)
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
  else if(strcmp(cmd,"get_default_window_cut_off_m") == 0)
  {
    plhs[0] = mxCreateDoubleScalar((double) WINDOW_HELP_ESTIMATE_m);
    return;
  }
  else if(strcmp(cmd,"get_epsilon") == 0)
  {
    plhs[0] = mxCreateDoubleScalar((double) X(float_property)(NFFT_EPSILON));
    return;
  }
  else
    mexErrMsgTxt("nfft: Unknown command.\n");
}
