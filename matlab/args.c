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
#include "nfft3.h"
#include "infft.h"
#include "imex.h"

int nfft_mex_get_int(const mxArray *p, const char *errmsg)
{
  DM(if (!mxIsDouble(p) || mxIsComplex(p) || mxGetM(p) != 1 || mxGetN(p) != 1)
    mexErrMsgTxt(errmsg);)
  return mxGetScalar(p);
}

double nfft_mex_get_double(const mxArray *p, const char *errmsg)
{
  DM(if (!mxIsDouble(p) || mxIsComplex(p) || mxGetM(p) != 1 || mxGetN(p) != 1)
    mexErrMsgTxt(errmsg);)
  return mxGetScalar(p);
}

void nfft_mex_get_nm(const mxArray *prhs[], int *n, int *m)
{
  int t = nfft_mex_get_int(prhs[1],"Input argument N must be a scalar.");
  DM(if ((t < 0) || (t%2!=0))
    mexErrMsgTxt("Input argument N must be non-negative and multiple of two.");)
  *n = t;
  t = nfft_mex_get_int(prhs[2],"Input argument M must be a scalar.");
  DM(if (t < 1)
    mexErrMsgTxt("Input argument M must be positive.");)
  *m = t;
}

void nfft_mex_get_nm_odd(const mxArray *prhs[], int *n, int *m)
{
  int t = nfft_mex_get_int(prhs[1],"Input argument N must be a scalar.");
  DM(if (t < 0)
    mexErrMsgTxt("Input argument N must be non-negative.");)
  *n = t;
  t = nfft_mex_get_int(prhs[2],"Input argument M must be a scalar.");
  DM(if (t < 1)
    mexErrMsgTxt("Input argument M must be positive.");)
  *m = t;
}

void nfft_mex_get_n1n2m(const mxArray *prhs[], int *n1, int *n2, int *m)
{
  int t = nfft_mex_get_int(prhs[1],"Input argument N1 must be a scalar.");
  DM(if ((t < 0) || (t%2!=0))
    mexErrMsgTxt("Input argument N1 must be non-negative and even.");)
  *n1 = t;

  t = nfft_mex_get_int(prhs[2],"Input argument N2 must be a scalar.");
  DM(if ((t < 0) || (t%2!=0))
    mexErrMsgTxt("Input argument N2 must be non-negative and even.");)
  *n2 = t;

  t = nfft_mex_get_int(prhs[3],"Input argument M must be a scalar.");
  DM(if (t < 1)
    mexErrMsgTxt("Input argument M must be positive.");)
  *m = t;
}

void nfft_mex_get_n1n2m_odd(const mxArray *prhs[], int *n1, int *n2, int *m)
{
  int t = nfft_mex_get_int(prhs[1],"Input argument N1 must be a scalar.");
  DM(if (t < 0)
    mexErrMsgTxt("Input argument N1 must be non-negative.");)
  *n1 = t;

  t = nfft_mex_get_int(prhs[2],"Input argument N2 must be a scalar.");
  DM(if (t < 0)
    mexErrMsgTxt("Input argument N2 must be non-negative.");)
  *n2 = t;

  t = nfft_mex_get_int(prhs[3],"Input argument M must be a scalar.");
  DM(if (t < 1)
    mexErrMsgTxt("Input argument M must be positive.");)
  *m = t;
}

void nfft_mex_get_n1n2n3m(const mxArray *prhs[], int *n1, int *n2, int *n3, int *m)
{
  int t = nfft_mex_get_int(prhs[1],"Input argument N1 must be a scalar.");
  DM(if ((t < 0) || (t%2!=0))
    mexErrMsgTxt("Input argument N1 must be non-negative and even.");)
  *n1 = t;

  t = nfft_mex_get_int(prhs[2],"Input argument N2 must be a scalar.");
  DM(if ((t < 0) || (t%2!=0))
    mexErrMsgTxt("Input argument N2 must be non-negative and even.");)
  *n2 = t;

  t = nfft_mex_get_int(prhs[3],"Input argument N3 must be a scalar.");
  DM(if ((t < 0) || (t%2!=0))
    mexErrMsgTxt("Input argument N3 must be non-negative and even.");)
  *n3 = t;

  t = nfft_mex_get_int(prhs[4],"Input argument M must be a scalar.");
  DM(if (t < 1)
    mexErrMsgTxt("Input argument M must be positive.");)
  *m = t;
}

void nfft_mex_get_n1n2n3m_odd(const mxArray *prhs[], int *n1, int *n2, int *n3, int *m)
{
  int t = nfft_mex_get_int(prhs[1],"Input argument N1 must be a scalar.");
  DM(if (t < 0)
    mexErrMsgTxt("Input argument N1 must be non-negative.");)
  *n1 = t;

  t = nfft_mex_get_int(prhs[2],"Input argument N2 must be a scalar.");
  DM(if (t < 0)
    mexErrMsgTxt("Input argument N2 must be non-negative.");)
  *n2 = t;

  t = nfft_mex_get_int(prhs[3],"Input argument N3 must be a scalar.");
  DM(if (t < 0)
    mexErrMsgTxt("Input argument N3 must be non-negative.");)
  *n3 = t;

  t = nfft_mex_get_int(prhs[4],"Input argument M must be a scalar.");
  DM(if (t < 1)
    mexErrMsgTxt("Input argument M must be positive.");)
  *m = t;
}

void nfft_mex_check_nargs(const int nrhs, const int n, const char* errmsg)
{
  DM(if (nrhs != n)
    mexErrMsgTxt(errmsg);)
}

int nfft_mex_set_num_threads_check(const int nrhs, const mxArray *prhs[], void **plans, const int plans_num_allocated)
{
  nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for set_num_threads.");

  int nthreads_new = nfft_mex_get_int(prhs[1],"Input argument nthreads must be a scalar.");

  if (nthreads_new < 1)
    mexErrMsgTxt("Number of threads must be at least 1.");

  if (nthreads_new > 1 && !X(has_threads_enabled)())
    mexErrMsgTxt("Threads are not enabled.");

  int nthreads_old = X(get_num_threads)();
  if (nthreads_new != nthreads_old)
  {
    int i;
    int is_plan_allocated = 0;
    for (i = 0; i < plans_num_allocated; i++)
      if (plans[i] != 0)
      {
        is_plan_allocated = 1;
        break;
      }
    if (is_plan_allocated)
      mexWarnMsgIdAndTxt("nfft:set_num_threads:plansAllocated","At least one plan is allocated. New number of threads may not affect the FFT step of any allocated plans.");
  }

  return nthreads_new;
}

