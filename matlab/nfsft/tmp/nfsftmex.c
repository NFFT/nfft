/* $Id: nfsftmex.c 2636 2008-12-01 08:31:10Z keiner $
 *
 * Copyright (c) 2007, 2008 Jens Keiner, Stefan Kunis, Daniel Potts
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

#include <complex.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "nfft3.h"
#include "infft.h"
#include "util.h"
#include "../imex.h"

#ifndef HAVE_MEXVERSION_C
  #include "mexversion.c"
#endif

#define PLANS_MAX 100 /* maximum number of plans */
#define CMD_LEN_MAX 20 /* maximum length of command argument */

/* global flags */
#define NFSFT_MEX_FIRST_CALL (1U << 0)
#define NFSFT_MEX_PRECOMPUTED (1U << 1)
static unsigned short gflags = NFSFT_MEX_FIRST_CALL;

static nfsft_plan* plans[PLANS_MAX]; /* plans */
static int n_max = -1; /* maximum degree precomputed */
static char cmd[CMD_LEN_MAX];

static inline int get_int(const mxArray *p, const char *errmsg)
{
  if (!mxIsDouble(p) || mxIsComplex(p) || mxGetM(p) != 1 || mxGetN(p) != 1)
    mexErrMsgTxt(errmsg);
  return mxGetScalar(p);
}

static inline double get_double(const mxArray *p, const char *errmsg)
{
  if (!mxIsDouble(p) || mxIsComplex(p) || mxGetM(p) != 1 || mxGetN(p) != 1)
    mexErrMsgTxt(errmsg);
  return mxGetScalar(p);
}

static inline void get_nm(const mxArray *prhs[], int *n, int *m)
{
  int t = get_int(prhs[1],"nfsft: Input argument N must be a scalar.");
  if (t < 0)
    mexErrMsgTxt("nfsft: Input argument N must be non-negative.");
  *n = t;
  t = get_int(prhs[2],"nfsft: Input argument M must be a scalar.");
  if (t < 1)
    mexErrMsgTxt("nfsft: Input argument M must be positive.");
  *m = t;
}

static inline void get_nmf(const mxArray *prhs[], int *n, int *m,
  unsigned int *f)
{
  get_nm(prhs,n,m);
  *f = get_int(prhs[3],"nfsft: Input argument flags must be a scalar.");
}

static inline void get_nmffc(const mxArray *prhs[], int *n, int *m,
  unsigned int *f, unsigned int *f2, int *c)
{
  get_nmf(prhs,n,m,f);
  *f2 = get_int(prhs[4],"nfsft: Input argument flags2 must be a scalar.");
  {
    int t = get_int(prhs[5],"nfsft: Input argument c must be a scalar.");
    if (t < 1)
      mexErrMsgTxt("nfsft: Input argument c must be positive.");
    *c = t;
  }
}

static inline void check_nargs(const int nrhs, const int n, const char* errmsg)
{
  if (nrhs != n)
    mexErrMsgTxt(errmsg);
}

static inline int mkplan()
{
  int i = 0;
  while (i < PLANS_MAX && plans[i] != 0) i++;
  if (i == PLANS_MAX)
    mexErrMsgTxt("nfsft: Too many plans already allocated.");
  plans[i] = nfft_malloc(sizeof(nfsft_plan));
  return i;
}

/* cleanup on mex function unload */
static void cleanup(void)
{
  int i;

  if (!(gflags & NFSFT_MEX_FIRST_CALL))
  {
    for (i = 0; i < PLANS_MAX; i++)
      if (plans[i])
      {
        nfsft_finalize(plans[i]);
        nfft_free(plans[i]);
        plans[i] = 0;
      }
    nfsft_forget();
    gflags |= NFSFT_MEX_FIRST_CALL;
    gflags &= ~NFSFT_MEX_PRECOMPUTED;
    n_max = -1;
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (gflags & NFSFT_MEX_FIRST_CALL)
  {
    /* Force Matlab to load libfftw3. There is at least one version of Matlab
     * which otherwise crashes upon invocation of this mex function. */
    mexEvalString("fft([1,2,3,4]);");

    install_mem_hooks();

    /* plan pointers to zeros */
    {
      int i;
      for (i = 0; i < PLANS_MAX; i++)
        plans[i] = 0;
    }

    mexAtExit(cleanup);
    gflags &= ~NFSFT_MEX_FIRST_CALL;
  }

  /* command string */
  if (nrhs == 0)
    mexErrMsgTxt("At least one input required.");

  if (!mxIsChar(prhs[0]))
    mexErrMsgTxt("First argument must be a string.");

  if (mxGetString(prhs[0], cmd, CMD_LEN_MAX))
    mexErrMsgTxt("Could not get command string.");

  if (strcmp(cmd,"init") == 0)
  {
    check_nargs(nrhs,3,"Wrong number of arguments for init.");
    {
      const int i = mkplan();
      int n, m;
      get_nm(prhs,&n,&m);
      nfsft_init(plans[i],n,m);
      plhs[0] = mxCreateDoubleScalar(i);
    }
    return;
  }
  else if (strcmp(cmd,"init_advanced") == 0)
  {
    check_nargs(nrhs,4,"Wrong number of arguments for init_advanced.");
    {
      const int i = mkplan();
      int n, m;
      unsigned int f;
      get_nmf(prhs,&n,&m,&f);
      nfsft_init_advanced(plans[i],n,m,f | NFSFT_MALLOC_X | NFSFT_MALLOC_F |
        NFSFT_MALLOC_F_HAT);
      plhs[0] = mxCreateDoubleScalar(i);
    }
//    mexPrintf("init_advanced: N = %d, M = %d, flags = %d, index = %d, "
//        "plan = %p, M_total = %d, N_total = %d\n",dp1[0],dp2[0],dp3[0],i,
//        plans[i],plans[i]->M_total,plans[i]->N_total);
    return;
  }
  else if (strcmp(cmd,"init_guru") == 0)
  {
    check_nargs(nrhs,6,"Wrong number of arguments for init_guru.");
    {
      const int i = mkplan();
      int n, m, c;
      unsigned int f, f2;
      get_nmffc(prhs,&n,&m,&f,&f2,&c);
      nfsft_init_guru(plans[i],n,m,f | NFSFT_MALLOC_X | NFSFT_MALLOC_F |
        NFSFT_MALLOC_F_HAT, PRE_PHI_HUT | PRE_PSI | FFTW_INIT
        | FFT_OUT_OF_PLACE, c);
      plhs[0] = mxCreateDoubleScalar(i);
    }
    return;
  }
  else if (strcmp(cmd,"precompute") == 0)
  {
    check_nargs(nrhs,5,"Wrong number of arguments for precompute.");
    {
      int n = get_int(prhs[1],"nfsft: Input argument n must be a scalar.");
      double k = get_double(prhs[2],"nfsft: Input argument kappa must be a scalar.");
      unsigned int f = get_int(prhs[3],"nfsft: Input argument flags must be a scalar.");
      unsigned int f2 = get_int(prhs[4],"nfsft: Input argument flags2 must be a scalar.");
      if (n_max < n)
      {
        if (gflags & NFSFT_MEX_PRECOMPUTED)
          nfsft_forget();

        nfsft_precompute(n,k,f,f2);
        n_max = n;
        gflags |= NFSFT_MEX_PRECOMPUTED;
      }
    }
    return;
  }
  else if (strcmp(cmd,"forget") == 0)
  {
    check_nargs(nrhs,1,"Wrong number of arguments for forget.");
    nfsft_forget();
    gflags &= ~NFSFT_MEX_PRECOMPUTED;
    n_max = -1;
    return;
  }
  else if (strcmp(cmd,"trafo") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for trafo.");
    {
      int i = get_int(prhs[1],"nfsft: Input argument plan must be a scalar.");
      nfsft_trafo(plans[i]);
    }
    return;
  }
  else if (strcmp(cmd,"adjoint") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for adjoint.");
    {
      int i = get_int(prhs[1],"nfsft: Input argument plan must be a scalar.");
      nfsft_adjoint(plans[i]);
    }
    return;
  }
  else if (strcmp(cmd,"finalize") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for finalize.");
    {
      int i = get_int(prhs[1],"nfsft: Input argument plan must be a scalar.");
      nfsft_finalize(plans[i]);
      nfft_free(plans[i]);
      plans[i] = 0;
    }
    return;
  }
  else if (strcmp(cmd,"precompute_x") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for precompute_x.");
    {
      int i = get_int(prhs[1],"nfsft: Input argument plan must be a scalar.");
      nfsft_precompute_x(plans[i]);
    }
    return;
  }
  else if (strcmp(cmd,"trafo_direct") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for trafo direct.");
    {
      int i = get_int(prhs[1],"nfsft: Input argument plan must be a scalar.");
      ndsft_trafo(plans[i]);
    }
    return;
  }
  else if (strcmp(cmd,"adjoint_direct") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for adjoint direct.");
    {
      int i = get_int(prhs[1],"nfsft: Input argument plan must be a scalar.");
      ndsft_adjoint(plans[i]);
    }
    return;
  }
  else if (strcmp(cmd,"get_x") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for get_x.");
    {
      int i = get_int(prhs[1],"nfsft: Input argument plan must be a scalar.");
      const int m = plans[i]->M_total;
      plhs[0] = mxCreateDoubleMatrix(2, m, mxREAL);
      {
        double *x = mxGetPr(plhs[0]);
        int j;
        for (j = 0; j < m; j++)
        {
          x[2*j] = 2*PI*plans[i]->x[2*j];
          if (x[2*j] < 0.0)
            x[2*j] += 2*PI;
          x[2*j+1] = 2*PI*plans[i]->x[2*j+1];
        }
      }
    }
    return;
  }
  else if (strcmp(cmd,"get_f") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for get_f.");
    {
      int i = get_int(prhs[1],"nfsft: Input argument plan must be a scalar.");
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
      int i = get_int(prhs[1],"nfsft: Input argument plan must be a scalar.");
      const int n = plans[i]->N;
      plhs[0] = mxCreateDoubleMatrix(2*n+1, n+1, mxCOMPLEX);
      {
        double *f_hatr = mxGetPr(plhs[0]), *f_hati = mxGetPi(plhs[0]);
        int k,j;
        for (k = 0; k <= n; k++)
          for (j = -k; j <= k; j++)
          {
            f_hatr[k*(2*n+1)+n+j] = creal(plans[i]->f_hat[NFSFT_INDEX(k,j,plans[i])]);
            f_hati[k*(2*n+1)+n+j] = cimag(plans[i]->f_hat[NFSFT_INDEX(k,j,plans[i])]);
          }
      }
    }
    return;
  }
  else if (strcmp(cmd,"get_f_hat_linear") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for get_f_hat_linear.");
    {
      int i = get_int(prhs[1],"nfsft: Input argument plan must be a scalar.");
      const int n = plans[i]->N;
      plhs[0] = mxCreateDoubleMatrix((n+1)*(n+1), 1, mxCOMPLEX);
      {
        double *f_hatr = mxGetPr(plhs[0]), *f_hati = mxGetPi(plhs[0]);
        int j,k,idx = 0;
        for (k = 0; k <= n; k++)
          for (j = -k; j <= k; j++, idx++)
          {
            f_hatr[idx] = creal(plans[i]->f_hat[NFSFT_INDEX(k,j,plans[i])]);
            f_hati[idx] = cimag(plans[i]->f_hat[NFSFT_INDEX(k,j,plans[i])]);
          }
      }
    }
    return;
  }
  else if (strcmp(cmd,"set_x") == 0)
  {
    check_nargs(nrhs,3,"Wrong number of arguments for set_x.");
    {
      int i = get_int(prhs[1],"nfsft: Input argument plan must be a scalar.");
      const int m = plans[i]->M_total;
      if (!mxIsDouble(prhs[2]) || mxGetNumberOfDimensions(prhs[2]) > 2)
        mexErrMsgTxt("Input argument x must be a 2 x M double array");
      if (mxGetM(prhs[2]) != 2 || mxGetN(prhs[2]) != m)
        mexErrMsgTxt("Input argument x must have correct size.");
      {
        double *x = mxGetPr(prhs[2]);
        int j;
        for (j = 0; j < m; j++)
        {
          plans[i]->x[2*j] = ((x[2*j] > PI)?(x[2*j] - 2*PI):(x[2*j]))/(2*PI);
          plans[i]->x[2*j+1] = x[2*j+1]/(2*PI);
        }
      }
    }
//    mexPrintf("set_x: nrhs = %d, index = %d, plan = %p, m_total = %d\n",nrhs,dp1[0],
//      plans[(int)(dp1[0])],plans[(int)(dp1[0])]->M_total);
    return;
  }
  else if (strcmp(cmd,"set_f") == 0)
  {
    check_nargs(nrhs,3,"Wrong number of arguments for set_f.");
    {
      int i = get_int(prhs[1],"nfsft: Input argument plan must be a scalar.");
      const int m = plans[i]->M_total;
      if (mxGetM(prhs[2]) * mxGetN(prhs[2]) != m)
        mexErrMsgTxt("Input argument f must have correct size.");
      {
        double *fr = mxGetPr(prhs[2]), *fi = mxGetPi(prhs[2]);
        int j;
        for (j = 0; j < m; j++)
          plans[i]->f[j] = fr[j] + _Complex_I*((fi)?(fi[j]):0.0);
      }
    }
    return;
  }
  else if (strcmp(cmd,"set_f_hat") == 0)
  {
    check_nargs(nrhs,3,"Wrong number of arguments for set_f_hat.");
    {
      int i = get_int(prhs[1],"nfsft: Input argument plan must be a scalar.");
      const int n = plans[i]->N;
      if (!mxIsDouble(prhs[2]))
        mexErrMsgTxt("Input argument f must be a double array");
      if (mxGetM(prhs[2]) != 2*n+1 || mxGetN(prhs[2]) != n+1)
        mexErrMsgTxt("Input argument f must have correct size.");
      {
        double *f_hatr = mxGetPr(prhs[2]), *f_hati = mxGetPi(prhs[2]);
        int j,k;
        for (k = 0; k <= n; k++)
          for (j = -k; j <= k; j++)
            plans[i]->f_hat[NFSFT_INDEX(k,j,plans[i])] =
              f_hatr[k*(2*n+1)+n+j] + _Complex_I*((f_hati)?(f_hati[k*(2*n+1)+n+j]):0.0);
      }
    }
    return;
  }
  else if (strcmp(cmd,"set_f_hat_linear") == 0)
  {
    check_nargs(nrhs,3,"Wrong number of arguments for set_f_hat_linear.");
    {
      int i = get_int(prhs[1],"nfsft: Input argument plan must be a scalar.");
      const int n = plans[i]->N;
      if (!mxIsDouble(prhs[2]))
        mexErrMsgTxt("Input argument f_hat must be a double array");
      if (mxGetM(prhs[2]) != (n+1)*(n+1) || mxGetN(prhs[2]) != 1)
        mexErrMsgTxt("Input argument f_hat must have correct size.");
      {
        double *f_hatr = mxGetPr(prhs[2]), *f_hati = mxGetPi(prhs[2]);
        int j,k,idx=0;
        for (k = 0; k <= n; k++)
          for (j = -k; j <= k; j++, idx++)
            plans[i]->f_hat[NFSFT_INDEX(k,j,plans[i])] =
              f_hatr[idx] + _Complex_I*((f_hati)?(f_hati[idx]):0.0);
      }
    }
    return;
  }
  else if (strcmp(cmd,"display") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for set_f_hat_linear.");
    {
      int i = get_int(prhs[1],"nfsft: Input argument plan must be a scalar.");
      mexPrintf("Plan %d\n",i);
      mexPrintf("  pointer: %p\n",plans[i]);
      mexPrintf("        N: %d\n",plans[i]->N);
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
    mexErrMsgTxt("nfsft: Unknown command.\n");
}

