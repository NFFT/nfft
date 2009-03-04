/*
 * Copyright (c) 2002, 2009 Jens Keiner, Daniel Potts, Stefan Kunis
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

/* $Id: nfftmex.c 3058 2009-03-03 14:01:44Z keiner $ */

#include <complex.h>
#include <string.h>
#include <stdio.h>
#include "nfft3.h"
#include "infft.h"
#include "util.h"
#include "imex.h"

#ifndef HAVE_MEXVERSION_C
  #include "mexversion.c"
#endif

#define PLANS_MAX 100 /* maximum number of plans */
#define CMD_LEN_MAX 20 /* maximum length of command argument */

/* global flags */
#define NFFT_MEX_FIRST_CALL (1U << 0)
#define NFFT_MEX_PRECOMPUTED (1U << 1)
unsigned short gflags = NFFT_MEX_FIRST_CALL;

nfft_plan* plans[PLANS_MAX]; /* plans */
int n_max = -1; /* maximum degree precomputed */
char cmd[CMD_LEN_MAX];

static inline void get_nm(const mxArray *prhs[], int *n, int *m)
{
  int t = nfft_mex_get_int(prhs[1],"nfft: Input argument N must be a scalar.");
  DM(if (t < 0)
    mexErrMsgTxt("nfft: Input argument N must be non-negative.");)
  *n = t;
  t = nfft_mex_get_int(prhs[2],"nfft: Input argument M must be a scalar.");
  DM(if (t < 1)
    mexErrMsgTxt("nfft: Input argument M must be positive.");)
  *m = t;
}

static inline void get_nmf(const mxArray *prhs[], int *n, int *m,
  unsigned int *f)
{
  get_nm(prhs,n,m);
  *f = nfft_mex_get_int(prhs[3],"nfft: Input argument flags must be a scalar.");
}

static inline void get_nmffc(const mxArray *prhs[], int *n, int *m,
  unsigned int *f, unsigned int *f2, int *c)
{
  get_nmf(prhs,n,m,f);
  *f2 = nfft_mex_get_int(prhs[4],"nfft: Input argument flags2 must be a scalar.");
  {
    int t = nfft_mex_get_int(prhs[5],"nfft: Input argument c must be a scalar.");
    DM(if (t < 1)
      mexErrMsgTxt("nfft: Input argument c must be positive.");)
    *c = t;
  }
}

static inline void check_nargs(const int nrhs, const int n, const char* errmsg)
{
  DM(if (nrhs != n)
    mexErrMsgTxt(errmsg);)
}

static inline int mkplan(void)
{
  int i = 0;
  while (i < PLANS_MAX && plans[i] != 0) i++;
  if (i == PLANS_MAX)
    mexErrMsgTxt("nfft: Too many plans already allocated.");
  plans[i] = nfft_malloc(sizeof(nfft_plan));
  return i;
}

/* cleanup on mex function unload */
static void cleanup(void)
{
  int i;

  if (!(gflags & NFFT_MEX_FIRST_CALL))
  {
    for (i = 0; i < PLANS_MAX; i++)
      if (plans[i])
      {
        nfft_finalize(plans[i]);
        plans[i] = 0;
      }
    gflags |= NFFT_MEX_FIRST_CALL;
    gflags &= ~NFFT_MEX_PRECOMPUTED;
    n_max = -1;
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

    /* plan pointers to zeros */
    {
      int i;
      for (i = 0; i < PLANS_MAX; i++)
        plans[i] = 0;
    }

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
    check_nargs(nrhs,3,"Wrong number of arguments for init.");
    {
      const int i = mkplan();
      int n, m;
      get_nm(prhs,&n,&m);
      nfft_init_1d(plans[i],n,m);
      plhs[0] = mxCreateDoubleScalar((double)i);
    }
    return;
  }
  else if (strcmp(cmd,"precompute_psi") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for precompute_one_psi.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      nfft_precompute_one_psi(plans[i]);
    }
    return;
  }
  else if (strcmp(cmd,"trafo") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for trafo.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      nfft_trafo(plans[i]);
    }
    return;
  }
  else if (strcmp(cmd,"adjoint") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for adjoint.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      nfft_adjoint(plans[i]);
    }
    return;
  }
  else if (strcmp(cmd,"finalize") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for finalize.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      nfft_finalize(plans[i]);
      nfft_free(plans[i]);
      plans[i] = 0;
    }
    return;
  }
  else if (strcmp(cmd,"trafo_direct") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for trafo direct.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      ndft_trafo(plans[i]);
    }
    return;
  }
  else if (strcmp(cmd,"adjoint_direct") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for adjoint direct.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      ndft_adjoint(plans[i]);
    }
    return;
  }
  else if (strcmp(cmd,"get_x") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for get_x.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      const int m = plans[i]->M_total;
      const int d = plans[i]->d;
      plhs[0] = mxCreateDoubleMatrix(d, m, mxREAL);
      {
        double *x = mxGetPr(plhs[0]);
        int j,t;
        for (j = 0; j < d*m; j++)
	  for (t = 0; t < d; t++)
	    x[d*j+t] = plans[i]->x[d*j+t];
      }
    }
    return;
  }
  else if (strcmp(cmd,"get_f") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for get_f.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      const int m = plans[i]->M_total;
      plhs[0] = mxCreateDoubleMatrix(m, 1, mxCOMPLEX);
      {
        double *fr = mxGetPr(plhs[0]), *fi = mxGetPi(plhs[0]);
        int j;
        for (j = 0; j < m; j++)
        {
          fr[j] = creal(plans[i]->f[j]);
          fi[j] = cimag(plans[i]->f[j]);
        }
      }
    }
    return;
  }
  else if (strcmp(cmd,"get_f_hat") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for get_f_hat.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      const int n = plans[i]->N_total;
      plhs[0] = mxCreateDoubleMatrix(n, 1, mxCOMPLEX);
      {
        double *f_hatr = mxGetPr(plhs[0]), *f_hati = mxGetPi(plhs[0]);
        int k;
        for (k = 0; k < n; k++)
        {
          f_hatr[k] = creal(plans[i]->f_hat[k]);
          f_hati[k] = cimag(plans[i]->f_hat[k]);
        }
      }
    }
    return;
  }
  else if (strcmp(cmd,"set_x") == 0)
  {
    check_nargs(nrhs,3,"Wrong number of arguments for set_x.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      const int m = plans[i]->M_total;
      const int d = plans[i]->d;
      DM(if (!mxIsDouble(prhs[2]) || mxGetNumberOfDimensions(prhs[2]) > 2)
        mexErrMsgTxt("Input argument x must be a d x M double array");)
      DM(if (mxGetM(prhs[2]) != d || mxGetN(prhs[2]) != m)
        mexErrMsgTxt("Input argument x must have correct size.");)
      {
        double *x = mxGetPr(prhs[2]);
        int j,t;
        for (j = 0; j < m; j++)
	  for (t = 0; t < d; t++)
	    plans[i]->x[d*j+t] = x[d*j+t];
      }
    }
    return;
  }
  else if (strcmp(cmd,"set_f") == 0)
  {
    check_nargs(nrhs,3,"Wrong number of arguments for set_f.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      const int m = plans[i]->M_total;
      DM(if (mxGetM(prhs[2]) * mxGetN(prhs[2]) != m)
        mexErrMsgTxt("Input argument f must have correct size.");)
      {
        double *fr = mxGetPr(prhs[2]), *fi = mxGetPi(prhs[2]);
        int j;
        if (fi)
          for (j = 0; j < m; j++)
            plans[i]->f[j] = fr[j] + _Complex_I*fi[j];
        else
          for (j = 0; j < m; j++)
            plans[i]->f[j] = fr[j];
      }
    }
    return;
  }
  else if (strcmp(cmd,"set_f_hat") == 0)
  {
    check_nargs(nrhs,3,"Wrong number of arguments for set_f_hat.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      const int n = plans[i]->N_total;
      DM(if (!mxIsDouble(prhs[2]))
        mexErrMsgTxt("Input argument f must be a double array");)
      DM(if (   mxGetM(prhs[2]) != n || mxGetN(prhs[2]) != 1)
        mexErrMsgTxt("Input argument f must have correct size.");)
      {
        double *f_hatr = mxGetPr(prhs[2]), *f_hati = mxGetPi(prhs[2]);
        int k;
        if (f_hati)
          for (k = 0; k < n; k++)
            plans[i]->f_hat[k] = f_hatr[k] + _Complex_I*f_hati[k];
        else
          for (k = 0; k < n; k++)
            plans[i]->f_hat[k] = f_hatr[k];
      }
    }
    return;
  }
  else if (strcmp(cmd,"display") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for set_f_hat_linear.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      mexPrintf("Plan %d\n",i);
      mexPrintf("  pointer: %p\n",plans[i]);
      mexPrintf("        d: %d\n",plans[i]->d);
      mexPrintf("  N_total: %d\n",plans[i]->N_total);
      mexPrintf("  M_total: %d\n",plans[i]->M_total);
      mexPrintf("        x: %p\n",plans[i]->x);
      mexPrintf("        f: %p\n",plans[i]->f);
      mexPrintf("    f_hat: %p\n",plans[i]->f_hat);
      mexPrintf("    flags: %d\n",plans[i]->nfft_flags);
    }
    return;
  }
  else
    mexErrMsgTxt("nfft: Unknown command.\n");
}
